#!/usr/bin/env python3
"""
Fetch PubMed metadata for a list of PMIDs and write a fully‑flattened CSV.

Changes (2025‑06‑08):
  • Populate DOI, PMCID, Abstract, Affiliations, Grants columns.
  • Removed Measurement Type column.
"""

import csv, sys, os, re, time, argparse
from datetime import datetime
from collections import OrderedDict
try:
    from Bio import Entrez
except ImportError:
    sys.exit("Biopython is required.  Install with `pip install biopython`.")

# ───────────────────────── CLI & helpers ─────────────────────────
def parse_arguments():
    p = argparse.ArgumentParser(description="Fetch PubMed metadata for a list of PMIDs.")
    p.add_argument("--pmid_study", required=True,
                   help="CSV/TSV/TXT file with rows: PMID, StudyIdentifier")
    p.add_argument("--email",    help="NCBI email")
    p.add_argument("--api_key",  help="NCBI API key")
    return p.parse_args()

def prompt_if_missing(a):
    if not a.email:    a.email    = input("Enter your NCBI email: ")
    if not a.api_key:  a.api_key  = input("Enter your NCBI API key (Enter to skip): ")
    return a

def detect_delimiter(path):
    with open(path, encoding="utf-8") as fh:
        sample = fh.read(2048)
    try:
        return csv.Sniffer().sniff(sample).delimiter
    except csv.Error:
        return ","

def load_pairs(path):
    if not os.path.exists(path):
        sys.exit(f"File not found: {path}")
    delim = detect_delimiter(path)
    with open(path, newline="", encoding="utf-8") as fh:
        return [(r[0].strip(), r[1].strip())
                for r in csv.reader(fh, delimiter=delim) if len(r) >= 2]

def fetch_pubmed(pmid):
    while True:
        try:
            h = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
            rec = Entrez.read(h);  h.close();  return rec
        except Exception:
            time.sleep(30)

def flatten(prefix, obj, out):
    if isinstance(obj, dict):
        for k, v in obj.items():
            flatten(f"{prefix}_{k}", v, out)
    elif isinstance(obj, list):
        if obj and isinstance(obj[0], dict):
            for i, v in enumerate(obj, 1):
                flatten(f"{prefix}_{i}", v, out)
        else:
            out[prefix] = "; ".join(str(x) for x in obj)
    else:
        out[prefix] = str(obj)

def detect_sequencing_type(keywords, mesh):
    txt = "; ".join(keywords + mesh).lower()
    if re.search(r"\b16s\b|\bamplicon\b", txt):   return "16S/amplicon"
    if "shotgun"      in txt:                    return "shotgun"
    if "metagenomic"  in txt:                    return "metagenomics"
    return ""

def write_with_retry(writer, row, outfile, header):
    for _ in range(3):
        try:
            writer.writerow({k: row.get(k, "") for k in header})
            return header
        except ValueError as e:
            if "dict contains fields not in fieldnames" in str(e):
                missing = eval(str(e).split(": ")[1])
                header  = sorted(set(header).union(missing))
                writer.file.close()
                with open(outfile, newline="", encoding="utf-8") as fh:
                    rows = list(csv.DictReader(fh))
                with open(outfile, "w", newline="", encoding="utf-8") as fh:
                    w = csv.DictWriter(fh, fieldnames=header)
                    w.writeheader()
                    for r in rows:  w.writerow({k: r.get(k, "") for k in header})
                writer = csv.DictWriter(open(outfile, "a", newline="", encoding="utf-8"),
                                        fieldnames=header)
            else:
                raise
    raise RuntimeError("Unable to write row after retries")

# ───────────────────────── main ─────────────────────────
def main():
    args = prompt_if_missing(parse_arguments())
    Entrez.email = args.email
    if args.api_key:  Entrez.api_key = args.api_key

    pairs      = load_pairs(args.pmid_study)
    outfile    = "publication_data.csv"

    core_cols = ["PMID","Study Identifier","Title","Authors","Consortium",
                 "Publication Year","Link","Sequencing","Abstract","Journal",
                 "Volume","Issue","Pages","DOI","PMCID","Publication Types",
                 "MeSH Terms","Keywords","Affiliations","Grants","Retrieved Date"]

    header     = list(core_cols)
    processed  = set()
    if os.path.exists(outfile):
        with open(outfile, newline="", encoding="utf-8") as fh:
            processed = {r["PMID"] for r in csv.DictReader(fh) if "PMID" in r}
            header    = sorted(set(header + next(csv.reader(open(outfile)), [])))

    file_exists = os.path.exists(outfile)
    fh          = open(outfile, "a" if file_exists else "w", newline="", encoding="utf-8")
    writer      = csv.DictWriter(fh, fieldnames=header)
    if not file_exists:  writer.writeheader()

    new_cnt = 0
    for pmid, study in pairs:
        if pmid in processed:  continue

        row = {k: "" for k in header}
        row.update({"PMID": pmid,
                    "Study Identifier": study,
                    "Retrieved Date": datetime.today().strftime("%Y-%m-%d")})

        if pmid.lower() == "na":
            header = write_with_retry(writer, row, outfile, header);  new_cnt += 1;  continue

        try:
            rec   = fetch_pubmed(pmid)
            art   = rec["PubmedArticle"][0]["MedlineCitation"]["Article"]
            cit   = rec["PubmedArticle"][0]["MedlineCitation"]
            pdata = rec["PubmedArticle"][0].get("PubmedData", {})

            pub_year  = art["Journal"]["JournalIssue"]["PubDate"].get("Year","")
            authors   = art.get("AuthorList", [])
            names     = ["{} {}".format(a.get("ForeName",""), a.get("LastName","")).strip()
                         for a in authors if "LastName" in a]

            # ----- core scalar fields ------------------------------------
            abstract = art.get("Abstract", {}).get("AbstractText", "")
            if isinstance(abstract, list):  abstract = " ".join(str(x) for x in abstract)

            doi = pmcid = ""
            for aid in pdata.get("ArticleIdList", []):
                t = aid.attributes.get("IdType", "")
                if t == "doi": doi = str(aid)
                if t == "pmc": pmcid = str(aid)

            mesh = [str(m["DescriptorName"]) for m in cit.get("MeshHeadingList",[])]
            kw   = [str(k) for grp in cit.get("KeywordList",[]) for k in grp]

            affiliations = []
            for a in authors:
                for aff in a.get("AffiliationInfo", []):
                    affiliations.append(aff.get("Affiliation",""))

            grants = []
            for g in cit.get("GrantList", []):
                agency = g.get("Agency", ""); gid = g.get("GrantID","")
                grants.append(f"{agency} ({gid})" if gid else agency)

            row.update({
                "Publication Year" : pub_year,
                "Title"            : art.get("ArticleTitle",""),
                "Authors"          : "; ".join(names),
                "Link"             : f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                "Journal"          : art["Journal"].get("Title",""),
                "Volume"           : art["Journal"]["JournalIssue"].get("Volume",""),
                "Issue"            : art["Journal"]["JournalIssue"].get("Issue",""),
                "Pages"            : art.get("Pagination",{}).get("MedlinePgn",""),
                "Publication Types": "; ".join(map(str, art.get("PublicationTypeList",[]))),
                "Abstract"         : abstract,
                "DOI"              : doi,
                "PMCID"            : pmcid,
                "MeSH Terms"       : "; ".join(mesh),
                "Keywords"         : "; ".join(kw),
                "Sequencing"       : detect_sequencing_type(kw, mesh),
                "Affiliations"     : "; ".join(affiliations),
                "Grants"           : "; ".join(grants)
            })

            # ----- flatten everything else --------------------------------
            flatten("Article",   art,   row)
            flatten("Citation",  cit,   row)
            flatten("PubmedData", pdata, row)

        except Exception as e:
            print(f"[ERROR] {pmid}: {e}")
            row["Title"] = "Error retrieving data"

        header = write_with_retry(writer, row, outfile, header)
        new_cnt += 1

    fh.close()
    print(f"Finished.  New PMIDs processed: {new_cnt}.  Output → {outfile}")

if __name__ == "__main__":
    main()

