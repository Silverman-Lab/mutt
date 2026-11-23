# convert_pkl.py
# 

import sys, types, importlib, os, zipfile, json, openpyxl

# 1) load pandas (so that its submodules get loaded into sys.modules)
import pandas as pd

# 2) capture the real numeric‐index module (if it still exists)
try:
    orig_numeric = importlib.import_module("pandas.core.indexes.numeric")
except ImportError:
    # pandas ≥2.0 may have reorganized things; fall back to pd.Index / pd.RangeIndex
    orig_numeric = types.SimpleNamespace(
        Int64Index   = pd.Index,
        UInt64Index  = pd.Index,
        Float64Index = pd.Index,
        RangeIndex   = pd.RangeIndex,
    )

# 3) build a fake module under the old name
fake = types.ModuleType("pandas.core.indexes.numeric")
fake.Int64Index   = orig_numeric.Int64Index
fake.UInt64Index  = getattr(orig_numeric, "UInt64Index", orig_numeric.Int64Index)
fake.Float64Index = getattr(orig_numeric, "Float64Index", orig_numeric.Int64Index)
fake.RangeIndex   = orig_numeric.RangeIndex

# 4) inject it so that pickle will find it
sys.modules["pandas.core.indexes.numeric"] = fake

import pickle

# 1) Set the working directory to where your ZIP lives
os.chdir("/home/maxwellkonnaris/Documents/pennstate/silvermanlab/mlscale/data_repository")

zip_path = "2021_barlow_microbiome_reimaginestudypatientsddpcr/dPCR data/dPCR_total_loads_260.pkl.zip"
csv_out  = "2021_barlow_microbiome_reimaginestudypatientsddpcr/dPCR_total_loads_260.csv"

if os.path.exists(zip_path):
    with zipfile.ZipFile(zip_path) as zf:
        name = next(f for f in zf.namelist() if f.endswith(".pkl"))
        with zf.open(name) as f:
            df = pickle.load(f)
    df.to_csv(csv_out, index=False)
    print(f"Saved CSV to {csv_out}")
else:
    print(f"⚠️  ZIP file not found: {zip_path}")
    print("I did something weird by unzipping the folder and then compressing the .pkl into a .zip (e.g. `zip dPCR_total_loads_260.pkl.zip dPCR_total_loads_260.pkl`) so do that and retry.")


# ─────── NEW: LOOP THROUGH pickle_files ───────
pickle_dir = "2021_barlow_microbiome_reimaginestudypatientsddpcr/pickle_files"
if not os.path.isdir(pickle_dir):
    raise RuntimeError(f"Directory not found: {pickle_dir} and might still be zipped so unzip the folder and retry")



for fname in os.listdir(pickle_dir):
    if not fname.lower().endswith(".pkl"):
        continue

    pkl_path = os.path.join(pickle_dir, fname)
    try:
        with open(pkl_path, "rb") as f:
            obj = pickle.load(f)
    except Exception as e:
        print(f"⚠️ Failed to load {fname}: {e}")
        continue

    used_xlsx = False

    # 1) If it's already a DataFrame (or similar), use it directly
    if hasattr(obj, "to_csv") and not isinstance(obj, (list, tuple)):
        df = obj

    # 2) If it's a list/tuple, try to frame it normally...
    elif isinstance(obj, (list, tuple)):
        try:
            df = pd.DataFrame(obj)
        except Exception:
            # ...but if that breaks (ragged), fall back to one-column JSON strings
            print(f"ℹ️  Inhomogeneous list in {fname}, writing .xlsx instead")
            rows = []
            for elem in obj:
                if hasattr(elem, "to_json"):
                    # DataFrame or Series → JSON
                    rows.append(elem.to_json())
                else:
                    try:
                        rows.append(json.dumps(elem))
                    except (TypeError, ValueError):
                        rows.append(repr(elem))
            df = pd.DataFrame({"value": rows})
            used_xlsx = True

    else:
        print(f"⚠️ Skipping unrecognized type {type(obj)} in {fname}")
        continue

    base = os.path.splitext(fname)[0]
    out_dir = pickle_dir
    if used_xlsx:
        out_path = os.path.join("2021_barlow_microbiome_reimaginestudypatientsddpcr/" + base + ".xlsx")
        df.to_excel(out_path, index=False)
        print(f"✅ Converted {fname} → {base}.xlsx")
    else:
        out_path = os.path.join("2021_barlow_microbiome_reimaginestudypatientsddpcr/" + base + ".csv")
        df.to_csv(out_path, index=False)
        print(f"✅ Converted {fname} → {base}.csv")