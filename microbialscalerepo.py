#!/usr/bin/env python3
"""
microbialscalerepo – Python analogue of the R helper

Requires:
    pip install rpy2 pandas tqdm
R side:
    install.packages(c("tidyverse", "optparse")) 
"""

import os
import pickle
from pathlib import Path
from typing import Dict, List, Mapping, MutableMapping, Optional, Sequence, Union

import pandas as pd
from tqdm import tqdm
import argparse

import rpy2.robjects as ro
from rpy2.robjects import pandas2ri

# ─────────────────────────────  R helpers  ──────────────────────────────
pandas2ri.activate()
r       = ro.r
new_env = r["new.env"]
sys_src = r["sys.source"]


def _source_parse_fun(parse_file: Union[str, Path], parser: str):
    env = new_env()
    sys_src(str(parse_file), env)
    fun_name = f"parse_{parser}"
    return env[fun_name] if fun_name in env else None


def _convert(obj):
    from rpy2.robjects.vectors import DataFrame, Matrix, ListVector

    if isinstance(obj, DataFrame):
        return pandas2ri.rpy2py(obj)
    if isinstance(obj, Matrix):
        return pd.DataFrame(
            pandas2ri.rpy2py(obj),
            index=ro.r["rownames"](obj),
            columns=ro.r["colnames"](obj),
        )
    if isinstance(obj, ListVector):
        return {k: _convert(v) for k, v in zip(obj.names, list(obj))}
    return obj  # primitives auto‑convert


# ───────────────────────────  Core function  ────────────────────────────
def microbialscalerepo(
    studies: Optional[Union[Sequence[str], Mapping[str, str]]] = None,
    base_directory: Union[str, Path] = ".",
    rawdata: bool = False,
    align_samples: bool = False,
    save_to: Optional[Union[str, Path]] = None,
) -> Dict[str, MutableMapping]:
    """See docstring in R version; behaviour replicated one‑for‑one."""
    base_directory = Path(base_directory)
    all_dirs = [
        p.name
        for p in base_directory.iterdir()
        if p.is_dir() and (p / "parse.R").exists()
    ]

    # ─── resolve selections ────────────────────────────────────────────
    if studies is None:
        selected, out_names = all_dirs, all_dirs
    elif isinstance(studies, Mapping):
        missing = set(studies.values()) - set(all_dirs)
        if missing:
            raise FileNotFoundError(f"Parser folder(s) not found: {', '.join(missing)}")
        selected, out_names = list(studies.values()), list(studies.keys())
    else:
        missing = set(studies) - set(all_dirs)
        if missing:
            raise FileNotFoundError(f"Parser folder(s) not found: {', '.join(missing)}")
        selected = out_names = list(studies)

    if not selected:
        return {}

    parsed: Dict[str, MutableMapping] = {}

    # ─── main loop ─────────────────────────────────────────────────────
    for parser in tqdm(selected, desc="Parsing studies", unit="study"):
        parse_fun = _source_parse_fun(base_directory / parser / "parse.R", parser)
        if parse_fun is None:
            print(f"[WARN] {parser}: function parse_{parser}() not found – skipped.")
            continue
        try:
            res_r = parse_fun(raw=rawdata)
        except Exception as e:
            print(f"[WARN] error inside {parser}: {e} – skipped.")
            continue

        res = _convert(res_r)
        if not isinstance(res, Mapping) or "counts" not in res:
            print(f"[WARN] {parser}: returned object lacks $counts – skipped.")
            continue

        # ─── sample alignment ──────────────────────────────────────────
        if align_samples and "scale" in res and isinstance(res["scale"], pd.DataFrame):
            id_col = res["scale"].columns[0]
            ids = res["scale"][id_col].dropna().astype(str)
            res["counts"] = res["counts"].loc[
                res["counts"].index.astype(str).isin(ids)
            ]
            res["scale"] = (
                res["scale"]
                .set_index(id_col)
                .reindex(res["counts"].index.astype(str))
                .reset_index()
            )

        parsed[out_names[selected.index(parser)]] = res

    # ─── optional pickle ───────────────────────────────────────────────
    if save_to is not None:
        save_to = Path(save_to)
        save_to.parent.mkdir(parents=True, exist_ok=True)
        with open(save_to, "wb") as fh:
            pickle.dump(parsed, fh)
        print(f"[INFO] parsed list pickled → {save_to}")

    return parsed


# ───────────────────────  Command‑line wrapper  ─────────────────────────
def _cli():
    p = argparse.ArgumentParser(
        description="Run microbialscalerepo (Python port) from the command line"
    )
    p.add_argument(
        "-s",
        "--studies",
        help="Comma‑separated list of parser folders OR of output names if --names given",
    )
    p.add_argument(
        "-n",
        "--names",
        help="Comma‑separated output names (must match --studies length)",
    )
    p.add_argument(
        "-b",
        "--base",
        default=".",
        help="Base directory containing parser sub‑folders [default %(default)s]",
    )
    p.add_argument(
        "-r", "--raw", action="store_true", help="Pass raw = TRUE to each parser"
    )
    p.add_argument(
        "-a",
        "--align",
        action="store_true",
        help="Align counts rows to sample IDs in scale",
    )
    p.add_argument(
        "-o",
        "--output",
        help="Path to save pickle (.pkl). If omitted, no file is written.",
    )

    args = p.parse_args()

    # ─── parse studies / names options ────────────────────────────────
    studs = (
        [s.strip() for s in args.studies.split(",")]
        if args.studies is not None
        else None
    )
    if args.names is not None:
        names = [n.strip() for n in args.names.split(",")]
        if len(names) != len(studs):
            p.error("--names and --studies must have the same number of entries")
        studs = dict(zip(names, studs))

    res = microbialscalerepo(
        studies=studs,
        base_directory=args.base,
        rawdata=args.raw,
        align_samples=args.align,
        save_to=args.output,
    )

    if args.output is None:
        print("Parsed studies:\n", "\n".join(res.keys()))


if __name__ == "__main__":
    _cli()
