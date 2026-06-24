#!/usr/bin/env python3
"""
23_external_similarity_druglike_selection_RDKit.py

Purpose
-------
Select a focused external ligand panel for JAK2 docking by structural similarity
to known JAK2 seed ligands and simple drug-likeness filters.

Input files
-----------
1) seed_jak2_ligands.csv
   Required columns: compound_name, smiles
   Optional columns: source, notes

2) external_library.csv OR external_library.smi
   CSV required columns: compound_name, smiles
   Optional columns: source, vendor_id
   SMI format: first column SMILES, second column compound_name/id

Output
------
external_similarity_screening_results.xlsx with sheets:
- selected_for_docking
- all_valid_screened
- excluded_invalid_or_failed
- seed_ligands_used
- summary

Dependencies
------------
conda install -c conda-forge rdkit pandas openpyxl
or use Google Colab with RDKit installed.

Example
-------
python 23_external_similarity_druglike_selection_RDKit.py \
  --seeds seed_jak2_ligands.csv \
  --library external_library.csv \
  --outdir external_similarity_screening \
  --similarity_cutoff 0.35 \
  --top_n 50
"""

import argparse
import os
import sys
import pandas as pd

try:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem, Descriptors, Lipinski, Crippen, rdMolDescriptors, QED
    from rdkit.Chem.MolStandardize import rdMolStandardize
    from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
except ImportError as e:
    sys.stderr.write("ERROR: RDKit is not installed. Install with: conda install -c conda-forge rdkit pandas openpyxl\n")
    raise e


def read_ligand_file(path):
    ext = os.path.splitext(path)[1].lower()
    if ext in [".csv", ".tsv", ".txt"]:
        sep = "\t" if ext == ".tsv" else ","
        df = pd.read_csv(path, sep=sep)
        lower = {c.lower(): c for c in df.columns}
        if "smiles" not in lower:
            raise ValueError(f"{path} must contain a 'smiles' column")
        if "compound_name" not in lower:
            # try common alternatives
            for alt in ["name", "compound", "id", "zinc_id", "chembl_id", "pubchem_cid"]:
                if alt in lower:
                    df = df.rename(columns={lower[alt]: "compound_name"})
                    break
            else:
                df["compound_name"] = [f"compound_{i+1}" for i in range(len(df))]
        df = df.rename(columns={lower.get("smiles", "smiles"): "smiles"})
        return df
    elif ext in [".smi", ".smiles"]:
        rows = []
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for i, line in enumerate(f):
                line = line.strip()
                if not line:
                    continue
                parts = line.split()
                smi = parts[0]
                name = parts[1] if len(parts) > 1 else f"compound_{i+1}"
                rows.append({"compound_name": name, "smiles": smi})
        return pd.DataFrame(rows)
    else:
        raise ValueError("Input library must be .csv, .tsv, .txt, .smi, or .smiles")


normalizer = rdMolStandardize.Normalizer()
uncharger = rdMolStandardize.Uncharger()
fragment_chooser = rdMolStandardize.LargestFragmentChooser()


def standardize_mol(smiles):
    if pd.isna(smiles) or str(smiles).strip() == "":
        return None, None
    mol = Chem.MolFromSmiles(str(smiles).strip())
    if mol is None:
        return None, None
    try:
        mol = fragment_chooser.choose(mol)
        mol = normalizer.normalize(mol)
        mol = uncharger.uncharge(mol)
        Chem.SanitizeMol(mol)
        can = Chem.MolToSmiles(mol, isomericSmiles=True)
        return mol, can
    except Exception:
        return None, None


def morgan_fp(mol, radius=2, nbits=2048):
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)


# PAINS filter catalog
params = FilterCatalogParams()
params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
pains_catalog = FilterCatalog(params)


def get_descriptors(mol):
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    rotb = Lipinski.NumRotatableBonds(mol)
    rings = Lipinski.RingCount(mol)
    qed = QED.qed(mol)
    pains_match = pains_catalog.GetFirstMatch(mol)
    pains = pains_match.GetDescription() if pains_match is not None else ""
    return {
        "MW": mw,
        "cLogP": logp,
        "HBD": hbd,
        "HBA": hba,
        "TPSA": tpsa,
        "RotB": rotb,
        "RingCount": rings,
        "QED": qed,
        "PAINS_alert": pains,
    }


def druglike_pass(row):
    # Broad kinase-inhibitor-compatible filters, not too strict.
    return (
        250 <= row["MW"] <= 650 and
        -1 <= row["cLogP"] <= 6 and
        row["HBD"] <= 5 and
        row["HBA"] <= 12 and
        row["TPSA"] <= 160 and
        row["RotB"] <= 12 and
        row["PAINS_alert"] == ""
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seeds", default="seed_jak2_ligands.csv")
    parser.add_argument("--library", default="external_library.csv")
    parser.add_argument("--outdir", default="external_similarity_screening")
    parser.add_argument("--similarity_cutoff", type=float, default=0.35)
    parser.add_argument("--top_n", type=int, default=50)
    parser.add_argument("--radius", type=int, default=2)
    parser.add_argument("--nbits", type=int, default=2048)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    seeds_raw = read_ligand_file(args.seeds)
    lib_raw = read_ligand_file(args.library)

    seed_rows = []
    seed_fps = []
    for _, r in seeds_raw.iterrows():
        mol, can = standardize_mol(r.get("smiles"))
        if mol is None:
            continue
        fp = morgan_fp(mol, args.radius, args.nbits)
        name = str(r.get("compound_name", "seed"))
        seed_rows.append({"seed_name": name, "seed_smiles": can, "source": r.get("source", "seed")})
        seed_fps.append((name, fp))

    if not seed_fps:
        raise ValueError("No valid seed ligands found. Check seed_jak2_ligands.csv")

    valid_rows = []
    bad_rows = []
    seen = set()
    for _, r in lib_raw.iterrows():
        name = str(r.get("compound_name", "compound"))
        smi = r.get("smiles")
        mol, can = standardize_mol(smi)
        if mol is None:
            bad_rows.append({"compound_name": name, "input_smiles": smi, "reason": "invalid_smiles_or_standardization_failed"})
            continue
        if can in seen:
            bad_rows.append({"compound_name": name, "input_smiles": smi, "standardized_smiles": can, "reason": "duplicate_standardized_smiles"})
            continue
        seen.add(can)
        fp = morgan_fp(mol, args.radius, args.nbits)
        sims = [(seed_name, DataStructs.TanimotoSimilarity(fp, seed_fp)) for seed_name, seed_fp in seed_fps]
        best_seed, best_sim = max(sims, key=lambda x: x[1])
        desc = get_descriptors(mol)
        row = {
            "compound_name": name,
            "source": r.get("source", "external"),
            "vendor_id": r.get("vendor_id", r.get("id", "")),
            "standardized_smiles": can,
            "best_seed_match": best_seed,
            "max_tanimoto_to_seed": best_sim,
        }
        row.update(desc)
        row["druglike_filter_pass"] = druglike_pass(row)
        valid_rows.append(row)

    valid = pd.DataFrame(valid_rows)
    bad = pd.DataFrame(bad_rows)
    seeds = pd.DataFrame(seed_rows)

    if len(valid) == 0:
        raise ValueError("No valid external molecules found after standardization")

    selected = valid[(valid["druglike_filter_pass"] == True) & (valid["max_tanimoto_to_seed"] >= args.similarity_cutoff)].copy()
    selected = selected.sort_values(["max_tanimoto_to_seed", "QED"], ascending=[False, False]).head(args.top_n)

    # Rank labels
    selected.insert(0, "selection_rank", range(1, len(selected) + 1))
    selected["selection_reason"] = (
        "External molecule: drug-like, no PAINS alert, and structurally similar to a known JAK2 seed ligand"
    )

    summary = pd.DataFrame({
        "Metric": [
            "Seed ligands valid",
            "External compounds input",
            "External compounds valid unique",
            "Drug-like valid compounds",
            f"Selected compounds similarity >= {args.similarity_cutoff}",
            "Top N requested",
        ],
        "Value": [
            len(seeds), len(lib_raw), len(valid), int(valid["druglike_filter_pass"].sum()), len(selected), args.top_n
        ]
    })

    out_xlsx = os.path.join(args.outdir, "external_similarity_screening_results.xlsx")
    out_csv = os.path.join(args.outdir, "selected_external_ligands_for_docking.csv")
    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
        selected.to_excel(writer, sheet_name="selected_for_docking", index=False)
        valid.sort_values("max_tanimoto_to_seed", ascending=False).to_excel(writer, sheet_name="all_valid_screened", index=False)
        bad.to_excel(writer, sheet_name="excluded_invalid_or_failed", index=False)
        seeds.to_excel(writer, sheet_name="seed_ligands_used", index=False)
        summary.to_excel(writer, sheet_name="summary", index=False)
    selected.to_csv(out_csv, index=False)

    print("Done.")
    print(f"Valid seed ligands: {len(seeds)}")
    print(f"Valid unique external compounds: {len(valid)}")
    print(f"Selected compounds: {len(selected)}")
    print(f"Output: {out_xlsx}")
    print(f"CSV: {out_csv}")


if __name__ == "__main__":
    main()
