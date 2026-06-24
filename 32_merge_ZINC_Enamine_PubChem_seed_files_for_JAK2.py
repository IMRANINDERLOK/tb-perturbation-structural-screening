#!/usr/bin/env python3
"""
Merge ZINC / Enamine / PubChem seed-based files for JAK2 external docking.

Input:
  - Any .smi/.smiles/.txt/.csv files containing SMILES + ID/name
  - Any .sdf files from Enamine/PubChem/ZINC
  - Filenames can look like:
      CHEMBL3699463-enamine.sdf
      CHEMBL3699463-zinc.smi
      Curcumin-zinc-2.smi
      momelotinib-enamine.sdf
      Quercetin-zinc.smi
      TG-101348-enamine.sdf
      pubchem.smi

Output:
  outdir/
    01_all_parsed_molecules.csv
    02_cleaned_unique_molecules.csv
    03_filtered_ranked_external_library.csv
    04_selected_external_ligands_for_docking.csv
    05_selected_external_ligands_for_docking.sdf
    06_summary_by_source_seed.csv
    07_method_note.txt

Run in Colab:
  !pip -q install rdkit pandas numpy tqdm openpyxl matplotlib
  !python 32_merge_ZINC_Enamine_PubChem_seed_files_for_JAK2.py --input_dir . --outdir JAK2_external_merged --select_total 150
"""

import argparse
import os
import re
import glob
import math
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, Crippen, Lipinski, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker

warnings.filterwarnings("ignore")

# Default seeds: user-provided / project seed ligands
DEFAULT_SEEDS = {
    "MOMELOTINIB": "N#CCNC(=O)C1=CC=C(C2=CC=NC(NC3=CC=C(N4CCOCC4)C=C3)=N2)C=C1",
    "CHEMBL3699463": "CS(=O)(=O)NC1=CC=C(C2=NC(NC3=CC=C(N4CCC(C(=O)NCCCCCC(=O)NO)CC4)C=C3)=NC=C2)C=C1",
    "TG-101348": "CC1=CN=C(NC2=CC=C(OCCN3CCCC3)C=C2)N=C1NC1=CC(S(=O)(=O)NC(C)(C)C)=CC=C1",
    "CURCUMIN": "COC1=CC(=CC=C1O)/C=C/C(=O)CC(=O)/C=C/C2=CC(=C(C=C2)O)OC",
    "QUERCETIN": "O=C1C(O)=C(Oc2cc(O)cc(O)c2C1=O)c1ccc(O)c(O)c1",
}

# ---------- helpers ----------

def norm_text(x):
    if x is None:
        return ""
    return str(x).strip()


def infer_source_and_seed(filepath):
    base = Path(filepath).stem
    b = base.lower()
    source = "Unknown"
    if "enamine" in b:
        source = "Enamine"
    elif "zinc" in b:
        source = "ZINC"
    elif "pubchem" in b or "pub_chem" in b:
        source = "PubChem"
    elif "chembl" in b:
        source = "ChEMBL/External"

    seed = "Unknown"
    if "momelotinib" in b:
        seed = "MOMELOTINIB"
    elif "tg-101348" in b or "tg_101348" in b or "tg101348" in b:
        seed = "TG-101348"
    elif "curcumin" in b:
        seed = "CURCUMIN"
    elif "quercetin" in b:
        seed = "QUERCETIN"
    elif "chembl3699463" in b or "chembl_3699463" in b:
        seed = "CHEMBL3699463"
    elif "pubchem" in b:
        seed = "PubChem_auto"

    return source, seed


def read_smi_like(filepath):
    source, seed = infer_source_and_seed(filepath)
    rows = []
    # Try CSV first if extension csv
    ext = Path(filepath).suffix.lower()
    if ext == ".csv":
        try:
            df = pd.read_csv(filepath)
            cols_lower = {c.lower(): c for c in df.columns}
            smiles_col = None
            for k in ["smiles", "canonical_smiles", "isomeric_smiles", "smile"]:
                if k in cols_lower:
                    smiles_col = cols_lower[k]
                    break
            if smiles_col is None:
                # fall back to first column
                smiles_col = df.columns[0]
            id_col = None
            for k in ["library_id", "compound_id", "id", "zinc_id", "name", "library_name", "cid"]:
                if k in cols_lower:
                    id_col = cols_lower[k]
                    break
            for i, r in df.iterrows():
                smi = norm_text(r.get(smiles_col, ""))
                cid = norm_text(r.get(id_col, f"{Path(filepath).stem}_{i+1}")) if id_col else f"{Path(filepath).stem}_{i+1}"
                if smi:
                    rows.append({"input_id": cid, "input_name": cid, "smiles_raw": smi, "input_file": Path(filepath).name, "source": source, "seed_group": seed})
            return rows
        except Exception as e:
            print(f"[WARN] CSV parse failed for {filepath}: {e}. Trying whitespace parse.")

    with open(filepath, "r", encoding="utf-8", errors="ignore") as fh:
        for i, line in enumerate(fh):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = re.split(r"[\t, ]+", line)
            if len(parts) == 0:
                continue
            smi = parts[0].strip()
            cid = parts[1].strip() if len(parts) > 1 else f"{Path(filepath).stem}_{i+1}"
            if smi.lower() in ["smiles", "canonical_smiles"]:
                continue
            rows.append({"input_id": cid, "input_name": cid, "smiles_raw": smi, "input_file": Path(filepath).name, "source": source, "seed_group": seed})
    return rows


def read_sdf(filepath):
    source, seed = infer_source_and_seed(filepath)
    rows = []
    suppl = Chem.SDMolSupplier(filepath, sanitize=False, removeHs=False)
    for i, mol in enumerate(suppl):
        if mol is None:
            continue
        # Try to sanitize later. Extract ID/name
        props = list(mol.GetPropNames())
        name = ""
        for p in ["ID", "id", "Name", "name", "Catalog ID", "CatalogID", "ZINC_ID", "PUBCHEM_COMPOUND_CID"]:
            if mol.HasProp(p):
                name = mol.GetProp(p)
                break
        if not name:
            name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"{Path(filepath).stem}_{i+1}"
        try:
            Chem.SanitizeMol(mol)
            smi = Chem.MolToSmiles(mol, isomericSmiles=True)
        except Exception:
            try:
                smi = Chem.MolToSmiles(mol, isomericSmiles=True)
            except Exception:
                continue
        rows.append({"input_id": norm_text(name), "input_name": norm_text(name), "smiles_raw": smi, "input_file": Path(filepath).name, "source": source, "seed_group": seed})
    return rows


def standardize_smiles(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None, None
    try:
        # disconnect metals and normalize
        mol = rdMolStandardize.MetalDisconnector().Disconnect(mol)
        mol = rdMolStandardize.Normalize(mol)
        # keep largest fragment
        chooser = rdMolStandardize.LargestFragmentChooser()
        mol = chooser.choose(mol)
        mol = rdMolStandardize.Uncharger().uncharge(mol)
        Chem.SanitizeMol(mol)
        can = Chem.MolToSmiles(mol, isomericSmiles=True)
        return mol, can
    except Exception:
        try:
            Chem.SanitizeMol(mol)
            can = Chem.MolToSmiles(mol, isomericSmiles=True)
            return mol, can
        except Exception:
            return None, None


def calc_props(mol):
    return {
        "MW": Descriptors.MolWt(mol),
        "LogP": Crippen.MolLogP(mol),
        "HBA": Lipinski.NumHAcceptors(mol),
        "HBD": Lipinski.NumHDonors(mol),
        "RotB": Lipinski.NumRotatableBonds(mol),
        "TPSA": rdMolDescriptors.CalcTPSA(mol),
        "HAC": mol.GetNumHeavyAtoms(),
        "Fsp3": rdMolDescriptors.CalcFractionCSP3(mol),
        "Rings": rdMolDescriptors.CalcNumRings(mol),
        "AromaticRings": rdMolDescriptors.CalcNumAromaticRings(mol),
    }


def is_druglike(p, strict=False):
    if strict:
        return (
            250 <= p["MW"] <= 550 and
            0.5 <= p["LogP"] <= 5.0 and
            1 <= p["HBA"] <= 10 and
            0 <= p["HBD"] <= 5 and
            p["RotB"] <= 10 and
            35 <= p["TPSA"] <= 140 and
            16 <= p["HAC"] <= 45
        )
    return (
        180 <= p["MW"] <= 650 and
        -1.0 <= p["LogP"] <= 6.0 and
        p["HBA"] <= 12 and
        p["HBD"] <= 6 and
        p["RotB"] <= 12 and
        p["TPSA"] <= 170 and
        10 <= p["HAC"] <= 55
    )


def pains_filter():
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
    return FilterCatalog(params)


def fp(mol):
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)


def tanimoto(f1, f2):
    return DataStructs.TanimotoSimilarity(f1, f2)


def make_seed_fps(seed_csv=None):
    seeds = dict(DEFAULT_SEEDS)
    if seed_csv and os.path.exists(seed_csv):
        sdf = pd.read_csv(seed_csv)
        cols = {c.lower(): c for c in sdf.columns}
        smi_col = cols.get("smiles") or cols.get("seed_smiles")
        name_col = cols.get("seed_name") or cols.get("compound") or cols.get("name")
        if smi_col:
            for i, r in sdf.iterrows():
                name = norm_text(r.get(name_col, f"seed_{i+1}")) if name_col else f"seed_{i+1}"
                smi = norm_text(r.get(smi_col, ""))
                if smi:
                    seeds[name.upper()] = smi
    seed_records = []
    for name, smi in seeds.items():
        mol, can = standardize_smiles(smi)
        if mol is not None:
            seed_records.append({"seed_name": name.upper(), "smiles": can, "mol": mol, "fp": fp(mol)})
    return seed_records


def greedy_diverse_select(df, fps_by_idx, select_n, similarity_cutoff=0.75):
    selected = []
    for idx in df.index.tolist():
        if len(selected) >= select_n:
            break
        f = fps_by_idx[idx]
        too_close = False
        for sidx in selected:
            if tanimoto(f, fps_by_idx[sidx]) >= similarity_cutoff:
                too_close = True
                break
        if not too_close:
            selected.append(idx)
    # If not enough selected, fill remaining by rank
    if len(selected) < select_n:
        for idx in df.index.tolist():
            if idx not in selected:
                selected.append(idx)
            if len(selected) >= select_n:
                break
    return selected


def safe_minmax_series(s, higher_better=True):
    s = pd.to_numeric(s, errors="coerce")
    mn, mx = s.min(), s.max()
    if pd.isna(mn) or pd.isna(mx) or mx == mn:
        return pd.Series(np.ones(len(s)) * 0.5, index=s.index)
    scaled = (s - mn) / (mx - mn)
    return scaled if higher_better else 1 - scaled


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input_dir", default=".")
    ap.add_argument("--outdir", default="JAK2_external_merged")
    ap.add_argument("--seed_csv", default="JAK2_seed_ligands_READY.csv")
    ap.add_argument("--select_total", type=int, default=150)
    ap.add_argument("--select_per_seed_source", type=int, default=25)
    ap.add_argument("--min_sim", type=float, default=0.35, help="global minimum max similarity")
    ap.add_argument("--max_sim", type=float, default=0.85, help="exclude molecules too similar to known seeds")
    ap.add_argument("--strict_druglike", action="store_true")
    ap.add_argument("--diversity_cutoff", type=float, default=0.75)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    files = []
    for pat in ["*.smi", "*.smiles", "*.txt", "*.csv", "*.sdf"]:
        files.extend(glob.glob(os.path.join(args.input_dir, pat)))
    # Exclude known script/readme/seed/template/result files
    skip_words = ["seed", "template", "readme", "ranked", "selected", "method", "summary"]
    input_files = []
    for f in files:
        b = Path(f).name.lower()
        if b.endswith(".py"):
            continue
        if any(w in b for w in skip_words):
            continue
        input_files.append(f)

    print("Detected input files:")
    for f in input_files:
        print(" -", Path(f).name)
    if not input_files:
        raise SystemExit("No real input files found. Upload .smi/.sdf files first. Templates are ignored.")

    all_rows = []
    for f in input_files:
        ext = Path(f).suffix.lower()
        print(f"Reading {Path(f).name} ...")
        try:
            if ext == ".sdf":
                rows = read_sdf(f)
            else:
                rows = read_smi_like(f)
            print(f"  parsed {len(rows)} rows")
            all_rows.extend(rows)
        except Exception as e:
            print(f"[WARN] failed reading {f}: {e}")
    all_df = pd.DataFrame(all_rows)
    if all_df.empty:
        raise SystemExit("No molecules parsed from input files.")
    all_df.to_csv(os.path.join(args.outdir, "01_all_parsed_molecules.csv"), index=False)

    # Prepare filters / seeds
    seed_records = make_seed_fps(args.seed_csv if os.path.exists(args.seed_csv) else None)
    seed_df = pd.DataFrame([{k:v for k,v in s.items() if k not in ['mol','fp']} for s in seed_records])
    seed_df.to_csv(os.path.join(args.outdir, "00_seed_ligands_used.csv"), index=False)
    print("Seeds used:", ", ".join(seed_df['seed_name'].tolist()))
    cat = pains_filter()

    clean_rows = []
    fps_by_idx = {}
    mols_by_idx = {}
    print("Standardizing, filtering, and scoring molecules...")
    for i, r in tqdm(all_df.iterrows(), total=len(all_df)):
        mol, can = standardize_smiles(r["smiles_raw"])
        if mol is None:
            continue
        p = calc_props(mol)
        druglike = is_druglike(p, strict=args.strict_druglike)
        if not druglike:
            continue
        # PAINS flag; remove if pains
        pains = cat.HasMatch(mol)
        if pains:
            continue
        fpm = fp(mol)
        sims = [(s["seed_name"], tanimoto(fpm, s["fp"])) for s in seed_records]
        best_seed, maxsim = max(sims, key=lambda x: x[1]) if sims else ("NA", np.nan)
        # Filename seed similarity if available
        seed_group = str(r.get("seed_group", "Unknown")).upper()
        seed_group_sim = np.nan
        for sn, sv in sims:
            if seed_group in sn or sn in seed_group:
                seed_group_sim = sv
                break
        row = dict(r)
        row.update({"canonical_smiles": can, "best_seed": best_seed, "max_seed_similarity": maxsim, "seed_group_similarity": seed_group_sim})
        row.update(p)
        clean_rows.append(row)
        idx = len(clean_rows)-1
        fps_by_idx[idx] = fpm
        mols_by_idx[idx] = mol

    clean_df = pd.DataFrame(clean_rows)
    if clean_df.empty:
        raise SystemExit("No molecules passed standardization/drug-like/PAINS filters. Try without --strict_druglike or use larger input.")

    # Deduplicate by canonical smiles, keep highest max similarity / first source
    clean_df = clean_df.sort_values(["max_seed_similarity"], ascending=False).drop_duplicates("canonical_smiles", keep="first").reset_index(drop=True)
    # Rebuild fps and mol dict after reset
    fps_by_idx = {}
    mols_by_idx = {}
    for idx, r in clean_df.iterrows():
        mol = Chem.MolFromSmiles(r["canonical_smiles"])
        fps_by_idx[idx] = fp(mol)
        mols_by_idx[idx] = mol
    clean_df.to_csv(os.path.join(args.outdir, "02_cleaned_unique_molecules.csv"), index=False)

    # Similarity window and novelty exclusion
    filtered = clean_df[(clean_df["max_seed_similarity"] >= args.min_sim) & (clean_df["max_seed_similarity"] <= args.max_sim)].copy()
    if filtered.empty:
        print("[WARN] No molecules in requested similarity window. Relaxing min similarity to 0.25.")
        filtered = clean_df[(clean_df["max_seed_similarity"] >= 0.25) & (clean_df["max_seed_similarity"] <= args.max_sim)].copy()
    # score optimal similarity centered around 0.55/0.60; include properties
    # Moderate similarity is preferred over too-low or too-high similarity.
    target_sim = 0.58
    filtered["similarity_window_score"] = 1 - (filtered["max_seed_similarity"] - target_sim).abs() / max(target_sim-args.min_sim, args.max_sim-target_sim)
    filtered["similarity_window_score"] = filtered["similarity_window_score"].clip(0, 1)
    filtered["MW_score"] = 1 - ((filtered["MW"] - 400).abs() / 250).clip(0, 1)
    filtered["TPSA_score"] = 1 - ((filtered["TPSA"] - 85).abs() / 100).clip(0, 1)
    filtered["RotB_score"] = 1 - (filtered["RotB"] / 12).clip(0, 1)
    filtered["External_Priority_Score"] = 100 * (
        0.55*filtered["similarity_window_score"] +
        0.15*filtered["MW_score"] +
        0.15*filtered["TPSA_score"] +
        0.15*filtered["RotB_score"]
    )
    filtered = filtered.sort_values(["External_Priority_Score", "max_seed_similarity"], ascending=False).reset_index(drop=True)
    filtered.to_csv(os.path.join(args.outdir, "03_filtered_ranked_external_library.csv"), index=False)

    # Selection: balanced per source+seed_group first, then global diverse fill
    selected_indices = []
    filtered_orig_index = filtered.index.tolist()
    # Need local df with index currently 0..N; fps_by_idx was from clean_df. map to canonical smiles -> fp by clean index
    # Create filtered-local fps dict
    local_fps = {}
    local_mols = {}
    for i, r in filtered.iterrows():
        mol = Chem.MolFromSmiles(r["canonical_smiles"])
        local_fps[i] = fp(mol)
        local_mols[i] = mol

    group_cols = ["source", "seed_group"]
    for (src, seed), g in filtered.groupby(group_cols, dropna=False):
        g = g.sort_values(["External_Priority_Score"], ascending=False)
        chosen = greedy_diverse_select(g, local_fps, select_n=min(args.select_per_seed_source, len(g)), similarity_cutoff=args.diversity_cutoff)
        selected_indices.extend(chosen)
        selected_indices = list(dict.fromkeys(selected_indices))
        if len(selected_indices) >= args.select_total:
            break

    # Fill by global rank/diversity
    if len(selected_indices) < args.select_total:
        remaining = filtered.loc[[i for i in filtered.index if i not in selected_indices]].sort_values("External_Priority_Score", ascending=False)
        # Combine selected plus remaining with diverse condition relative to selected
        for idx in remaining.index.tolist():
            if len(selected_indices) >= args.select_total:
                break
            f = local_fps[idx]
            too_close = False
            for sidx in selected_indices:
                if tanimoto(f, local_fps[sidx]) >= args.diversity_cutoff:
                    too_close = True
                    break
            if not too_close:
                selected_indices.append(idx)
        # if still short, fill by rank
        if len(selected_indices) < args.select_total:
            for idx in remaining.index.tolist():
                if idx not in selected_indices:
                    selected_indices.append(idx)
                if len(selected_indices) >= args.select_total:
                    break

    selected = filtered.loc[selected_indices].copy().sort_values("External_Priority_Score", ascending=False).reset_index(drop=True)
    selected["Docking_Status"] = "Selected_for_LigPrep_Glide"
    selected.to_csv(os.path.join(args.outdir, "04_selected_external_ligands_for_docking.csv"), index=False)

    # Write SDF
    writer = Chem.SDWriter(os.path.join(args.outdir, "05_selected_external_ligands_for_docking.sdf"))
    for _, r in selected.iterrows():
        mol = Chem.MolFromSmiles(r["canonical_smiles"])
        if mol is None:
            continue
        mol = Chem.AddHs(mol)
        # Generate 3D coordinates if possible
        try:
            AllChem.EmbedMolecule(mol, randomSeed=42, maxAttempts=20)
            AllChem.UFFOptimizeMolecule(mol, maxIters=200)
        except Exception:
            mol = Chem.RemoveHs(mol)
        mol.SetProp("_Name", str(r.get("input_id", "external_ligand")))
        for prop in ["input_id", "input_name", "source", "seed_group", "best_seed", "max_seed_similarity", "External_Priority_Score", "canonical_smiles"]:
            if prop in r:
                mol.SetProp(prop, str(r[prop]))
        writer.write(mol)
    writer.close()

    summary = filtered.groupby(["source", "seed_group"], dropna=False).agg(
        molecules_ranked=("canonical_smiles", "count"),
        mean_similarity=("max_seed_similarity", "mean"),
        max_similarity=("max_seed_similarity", "max"),
        mean_priority=("External_Priority_Score", "mean")
    ).reset_index()
    sel_summary = selected.groupby(["source", "seed_group"], dropna=False).size().reset_index(name="selected_n")
    summary = summary.merge(sel_summary, on=["source", "seed_group"], how="left")
    summary["selected_n"] = summary["selected_n"].fillna(0).astype(int)
    summary.to_csv(os.path.join(args.outdir, "06_summary_by_source_seed.csv"), index=False)

    method = f"""
JAK2 external-library preparation method
======================================
Input directory: {args.input_dir}
Input files detected: {len(input_files)}
Total parsed molecules: {len(all_df)}
Cleaned unique molecules: {len(clean_df)}
Filtered ranked molecules: {len(filtered)}
Selected molecules for docking: {len(selected)}

Filtering:
- largest fragment kept, salts/metal disconnected where possible, neutralized where possible
- duplicate canonical SMILES removed
- PAINS A/B/C removed using RDKit FilterCatalog
- drug-like filter {'strict' if args.strict_druglike else 'standard'} applied
- Morgan fingerprint radius=2, 2048 bits
- max Tanimoto similarity to seed ligands calculated
- selected similarity window: {args.min_sim} to {args.max_sim}
- moderate similarity was prioritized to avoid compounds that are too distant or too close to known JAK2 seed ligands
- diverse selection used pairwise Tanimoto cutoff: {args.diversity_cutoff}

Main output for Maestro/LigPrep:
05_selected_external_ligands_for_docking.sdf

Recommended next step:
LigPrep -> Glide HTVS/SP/XP -> MM-GBSA -> select best external predicted hit for MD.
""".strip()
    with open(os.path.join(args.outdir, "07_method_note.txt"), "w", encoding="utf-8") as f:
        f.write(method + "\n")

    print("\nDONE")
    print("Parsed molecules:", len(all_df))
    print("Cleaned unique molecules:", len(clean_df))
    print("Filtered ranked molecules:", len(filtered))
    print("Selected for docking:", len(selected))
    print("\nOutput folder:", args.outdir)
    print("Important files:")
    print(" - 04_selected_external_ligands_for_docking.csv")
    print(" - 05_selected_external_ligands_for_docking.sdf")
    print(" - 06_summary_by_source_seed.csv")

if __name__ == "__main__":
    main()
