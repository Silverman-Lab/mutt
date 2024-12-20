import os
import subprocess
import json
import csv
import sys

def main():
    repo_path = os.path.abspath(".")
    output_csv = "silvermanrepo_parse_results.csv"
    results = []

    headers = ["study", "scale", "metadata", "proportions", "taxa", "sequences", "phylo"]

    for entry in os.listdir(repo_path):
        full_path = os.path.join(repo_path, entry)
        
        if os.path.isdir(full_path) and entry not in [".git", ".github", "__pycache__"]:
            parse_r_path = os.path.join(full_path, "parse.R")
            if os.path.isfile(parse_r_path):
                print(f"Testing parse in folder: {entry}", file=sys.stderr)
                
                try:
                    parse_output = subprocess.check_output(
                        ["Rscript", parse_r_path],
                        cwd=full_path,
                        stderr=subprocess.STDOUT,
                        text=True
                    ).strip()
                except subprocess.CalledProcessError as e:
                    print(f"Error running parse.R in {entry}:\n{e.output}", file=sys.stderr)
                    results.append([entry, 0, 0, 0, 0])
                    continue
                
                try:
                    data = json.loads(parse_output)
                except json.JSONDecodeError:
                    print(f"Failed to decode JSON in {entry}. Output was:\n{parse_output}", file=sys.stderr)
                    results.append([entry, 0, 0, 0, 0])
                    continue

                def has_data(obj):
                    if isinstance(obj, dict) and "data" in obj:
                        return len(obj["data"]) > 0
                    if isinstance(obj, list):
                        return len(obj) > 0
                    return bool(obj)

                scale_ok = 1 if "scale" in data and has_data(data["scale"]) else 0
                metadata_ok = 1 if "metadata" in data and has_data(data["metadata"]) else 0
                proportions_ok = 1 if "proportions" in data and has_data(data["proportions"]) else 0
                taxa_ok = 1 if "tax" in data and has_data(data["tax"]) else 0
                sequences_ok = 1 if "sequences" in data and has_data(data["sequences"]) else 0
                phylo_ok = 1 if "phylo" in data and has_data(data["phylo"]) else 0

                results.append([entry, scale_ok, metadata_ok, proportions_ok, taxa_ok, sequences_ok, phylo_ok])

    with open(output_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(headers)  
        writer.writerows(results) 

    print(f"Results saved to {output_csv}")

if __name__ == "__main__":
    main()
