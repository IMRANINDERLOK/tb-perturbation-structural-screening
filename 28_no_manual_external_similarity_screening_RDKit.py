#!/usr/bin/env python3
"""
28_no_manual_external_similarity_screening_RDKit.py

Purpose
-------
Automatic external-library similarity screening for JAK2 ATP-pocket docking.
No manual seed curation is required: default seeds are embedded in the script.

Main idea
---------
1. Use known/reference JAK2-pocket seeds to retain ATP-pocket compatibility.
2. Use non-classical CMap seeds to retain TB-signature/CMap novelty.
3. Exclude exact seed molecules and molecules that are too similar to known JAK2 drugs.
4. Select external compounds in a moderate-similarity novelty window for docking.

Required input
--------------
external_library.csv with at least one SMILES column. Acceptable SMILES column names:
    smiles, SMILES, canonical_smiles, Canonical_SMILES, structure, mol_smiles
Optional columns:
    library_id, id, compound_id, zinc_id, pubchem_cid
    library_name, name, compound_name
    source

Google Colab example
--------------------
!pip -q install rdkit-pypi pandas numpy matplotlib tqdm openpyxl
!python 28_no_manual_external_similarity_screening_RDKit.py \
    --library external_library.csv \
    --outdir external_similarity_no_manual \
    --select_n 100 \
    --top_n 500

Outputs
-------
01_default_seed_compounds_used.csv
02_external_library_standardized_filtered.csv
03_similarity_ranked_all.csv
04_selected_external_ligands_for_docking.csv
05_selected_external_ligands_for_docking.sdf
06_similarity_distribution.png
07_similarity_source_scatter.png
08_method_note_no_manual_external_similarity.txt
"""

import argparse
import os
import sys
import math
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, Crippen, Lipinski, rdMolDescriptors, QED
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

# -------------------------------------------------------------------------
# Default seeds: no manual curation needed.
# User-provided JAK2-pocket seeds plus non-classical CMap seeds.
# -------------------------------------------------------------------------
DEFAULT_SEEDS = [
    {
        "seed_name": "Momelotinib",
        "seed_role": "Native/reference JAK2 ligand",
        "seed_group": "known_JAK2_reference",
        "smiles": "N#CCNC(=O)C1=CC=C(C2=CC=NC(NC3=CC=C(N4CCOCC4)C=C3)=N2)C=C1",
    },
    {
        "seed_name": "TG-101348",
        "seed_role": "CMap/JAK2-linked comparator; fedratinib-like",
        "seed_group": "known_JAK2_reference",
        "smiles": "CC1=CN=C(NC2=CC=C(OCCN3CCCC3)C=C2)N=C1NC1=CC(S(=O)(=O)NC(C)(C)C)=CC=C1",
    },
    {
        "seed_name": "CHEMBL3699463",
        "seed_role": "ChEMBL JAK2 docking/MM-GBSA hit",
        "seed_group": "known_JAK2_reference",
        "smiles": "CS(=O)(=O)NC1=CC=C(C2=NC(NC3=CC=C(N4CCC(C(=O)NCCCCCC(=O)NO)CC4)C=C3)=NC=C2)C=C1",
    },
    {
        "seed_name": "Curcumin",
        "seed_role": "Non-classical CMap/TB-signature candidate",
        "seed_group": "nonclassical_CMap",
        "smiles": "COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O",
    },
    {
        "seed_name": "Quercetin",
        "seed_role": "Non-classical CMap/TB-signature candidate",
        "seed_group": "nonclassical_CMap",
        "smiles": "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12",
    },
]

SMILES_COL_CANDIDATES = [
    "smiles", "SMILES", "canonical_smiles", "Canonical_SMILES", "canonical SMILES",
    "structure", "mol_smiles", "Molecule_SMILES", "isomeric_smiles"
]
ID_COL_CANDIDATES = ["library_id", "id", "compound_id", "zinc_id", "ZINC_ID", "pubchem_cid", "CID", "chembl_id", "ChEMBL_ID"]
NAME_COL_CANDIDATES = ["library_name", "name", "compound_name", "Compound_Name", "drug_name", "Drug name"]
SOURCE_COL_CANDIDATES = ["source", "Source", "library_source"]

# -------------------------------------------------------------------------
# Standardization helpers
# -------------------------------------------------------------------------
def pick_col(columns, candidates):
    lower_map = {str(c).lower(): c for c in columns}
    for cand in candidates:
        if cand in columns:
            return cand
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    return None


def largest_fragment(mol):
    try:
        chooser = rdMolStandardize.LargestFragmentChooser()
        return chooser.choose(mol)
    except Exception:
        frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
        if not frags:
            return None
        return max(frags, key=lambda m: m.GetNumHeavyAtoms())


def standardize_smiles(smiles):
    if pd.isna(smiles):
        return None, None
    smiles = str(smiles).strip()
    if smiles == "" or smiles.lower() in {"nan", "none", "na"}:
        return None, None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None, None
        mol = largest_fragment(mol)
        if mol is None:
            return None, None
        # basic cleanup/uncharging
        try:
            mol = rdMolStandardize.Cleanup(mol)
            uncharger = rdMolStandardize.Uncharger()
            mol = uncharger.uncharge(mol)
        except Exception:
            pass
        Chem.SanitizeMol(mol)
        can = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        return mol, can
    except Exception:
        return None, None


def morgan_fp(mol, radius=2, nbits=2048):
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)


def calc_descriptors(mol):
    try:
        return {
            "MW": float(Descriptors.MolWt(mol)),
            "LogP": float(Crippen.MolLogP(mol)),
            "HBD": int(Lipinski.NumHDonors(mol)),
            "HBA": int(Lipinski.NumHAcceptors(mol)),
            "TPSA": float(rdMolDescriptors.CalcTPSA(mol)),
            "RotB": int(Lipinski.NumRotatableBonds(mol)),
            "Rings": int(rdMolDescriptors.CalcNumRings(mol)),
            "HeavyAtoms": int(mol.GetNumHeavyAtoms()),
            "QED": float(QED.qed(mol)),
        }
    except Exception:
        return {"MW": np.nan, "LogP": np.nan, "HBD": np.nan, "HBA": np.nan, "TPSA": np.nan,
                "RotB": np.nan, "Rings": np.nan, "HeavyAtoms": np.nan, "QED": np.nan}


def lipinski_veber_score(row):
    score = 0
    # Soft scoring, not hard medicinal claim
    if row["MW"] <= 500: score += 1
    if row["LogP"] <= 5: score += 1
    if row["HBD"] <= 5: score += 1
    if row["HBA"] <= 10: score += 1
    if row["TPSA"] <= 140: score += 1
    if row["RotB"] <= 10: score += 1
    return score / 6.0


def make_pains_catalog():
    params = FilterCatalogParams()
    try:
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
        return FilterCatalog(params)
    except Exception:
        return None


def pains_flag(mol, catalog):
    if catalog is None:
        return False
    try:
        return bool(catalog.HasMatch(mol))
    except Exception:
        return False

# -------------------------------------------------------------------------
# Similarity scoring
# -------------------------------------------------------------------------
def window_score(sim, low, high, optimal=None):
    """Return 0-1 score favoring moderate similarity window."""
    if sim is None or np.isnan(sim):
        return 0.0
    sim = float(sim)
    if optimal is None:
        optimal = (low + high) / 2.0
    if sim < low:
        return max(0.0, sim / low * 0.5)
    if low <= sim <= high:
        # peak at optimal; still acceptable near low/high
        width = max(optimal - low, high - optimal, 1e-6)
        return max(0.5, 1.0 - 0.5 * abs(sim - optimal) / width)
    # above high: too similar, drop quickly
    return max(0.0, 0.5 * (1.0 - (sim - high) / max(1.0 - high, 1e-6)))


def compute_selection_score(row):
    # Moderate JAK2 similarity for pocket compatibility.
    jak2_component = window_score(row["Max_Known_JAK2_Similarity"], 0.30, 0.75, optimal=0.55)
    # Non-classical CMap similarity for novelty/phenotype bridge.
    cmap_component = window_score(row["Max_Nonclassical_CMap_Similarity"], 0.25, 0.80, optimal=0.50)
    druglike_component = row.get("Lipinski_Veber_Score", 0.0)
    qed_component = row.get("QED", 0.0)
    # Penalty if too similar to known JAK2 seed.
    too_similar_penalty = 0.0
    if row["Max_Known_JAK2_Similarity"] >= 0.85:
        too_similar_penalty = 0.35
    if row.get("PAINS_Flag", False):
        too_similar_penalty += 0.10
    score = (
        0.42 * jak2_component +
        0.22 * cmap_component +
        0.18 * druglike_component +
        0.18 * qed_component -
        too_similar_penalty
    )
    return max(0.0, min(1.0, float(score)))

# -------------------------------------------------------------------------
# Main
# -------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="No-manual-curation external similarity screening for JAK2 docking.")
    parser.add_argument("--library", required=True, help="External library CSV with SMILES column")
    parser.add_argument("--outdir", default="external_similarity_no_manual", help="Output directory")
    parser.add_argument("--top_n", type=int, default=500, help="Number of ranked molecules to keep in top table")
    parser.add_argument("--select_n", type=int, default=100, help="Number of selected molecules for docking")
    parser.add_argument("--min_jak2_sim", type=float, default=0.30, help="Minimum JAK2-reference similarity for automatic selection")
    parser.add_argument("--max_jak2_sim", type=float, default=0.85, help="Maximum allowed similarity to known JAK2 seeds")
    parser.add_argument("--allow_pains", action="store_true", help="Do not exclude PAINS-flagged compounds from selected docking set")
    parser.add_argument("--strict_druglike", action="store_true", help="Apply stricter drug-like filters")
    parser.add_argument("--make_3d", action="store_true", help="Generate 3D conformers in SDF; otherwise SDF is 2D only")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Prepare seeds
    seed_rows = []
    seed_mols = []
    seed_fps = []
    seed_canonical_set = set()
    for s in DEFAULT_SEEDS:
        mol, can = standardize_smiles(s["smiles"])
        if mol is None:
            print(f"WARNING: Seed failed to parse: {s['seed_name']}", file=sys.stderr)
            continue
        fp = morgan_fp(mol)
        row = dict(s)
        row["canonical_smiles"] = can
        row.update(calc_descriptors(mol))
        seed_rows.append(row)
        seed_mols.append((s, mol))
        seed_fps.append((s, fp))
        seed_canonical_set.add(can)

    seeds_df = pd.DataFrame(seed_rows)
    seeds_df.to_csv(outdir / "01_default_seed_compounds_used.csv", index=False)

    if seeds_df.empty:
        raise RuntimeError("No valid seeds. Check embedded seed SMILES.")

    # Load library
    lib = pd.read_csv(args.library)
    smiles_col = pick_col(lib.columns, SMILES_COL_CANDIDATES)
    id_col = pick_col(lib.columns, ID_COL_CANDIDATES)
    name_col = pick_col(lib.columns, NAME_COL_CANDIDATES)
    source_col = pick_col(lib.columns, SOURCE_COL_CANDIDATES)

    if smiles_col is None:
        raise ValueError(f"No SMILES column found. Columns available: {list(lib.columns)}")

    print(f"Using SMILES column: {smiles_col}")
    print(f"Using ID column: {id_col if id_col else 'generated'}")
    print(f"Using name column: {name_col if name_col else 'generated'}")

    rows = []
    catalog = make_pains_catalog()

    for idx, r in tqdm(lib.iterrows(), total=len(lib), desc="Standardizing library"):
        mol, can = standardize_smiles(r[smiles_col])
        if mol is None:
            continue
        desc = calc_descriptors(mol)
        # broad filters before similarity
        if desc["HeavyAtoms"] < 12:
            continue
        if desc["MW"] < 180 or desc["MW"] > 750:
            continue
        if desc["HBA"] > 15 or desc["HBD"] > 8:
            continue
        if desc["RotB"] > 18:
            continue
        if desc["TPSA"] > 220:
            continue
        if args.strict_druglike:
            if not (250 <= desc["MW"] <= 650):
                continue
            if desc["LogP"] > 6 or desc["HBA"] > 12 or desc["HBD"] > 6 or desc["RotB"] > 12 or desc["TPSA"] > 160:
                continue
        fp = morgan_fp(mol)
        sims = []
        for seed, sfp in seed_fps:
            sims.append({
                "seed_name": seed["seed_name"],
                "seed_group": seed["seed_group"],
                "similarity": float(DataStructs.TanimotoSimilarity(fp, sfp))
            })
        sims_df = pd.DataFrame(sims)
        best = sims_df.sort_values("similarity", ascending=False).iloc[0]
        max_jak2 = sims_df.loc[sims_df["seed_group"] == "known_JAK2_reference", "similarity"].max()
        max_cmap = sims_df.loc[sims_df["seed_group"] == "nonclassical_CMap", "similarity"].max()
        pains = pains_flag(mol, catalog)
        out = {
            "Library_ID": str(r[id_col]) if id_col else f"LIB_{idx+1}",
            "Library_Name": str(r[name_col]) if name_col else f"compound_{idx+1}",
            "Library_Source": str(r[source_col]) if source_col else "external_library",
            "Input_SMILES": str(r[smiles_col]),
            "Canonical_SMILES": can,
            "Best_Seed": best["seed_name"],
            "Best_Seed_Group": best["seed_group"],
            "Max_Seed_Similarity": float(best["similarity"]),
            "Max_Known_JAK2_Similarity": float(max_jak2),
            "Max_Nonclassical_CMap_Similarity": float(max_cmap),
            "Exact_Seed_Match": can in seed_canonical_set,
            "PAINS_Flag": pains,
        }
        out.update(desc)
        rows.append(out)

    df = pd.DataFrame(rows)
    if df.empty:
        raise RuntimeError("No compounds remained after standardization and broad filters.")

    # Deduplicate by canonical smiles
    df = df.sort_values(["Max_Seed_Similarity"], ascending=False).drop_duplicates("Canonical_SMILES", keep="first").reset_index(drop=True)
    df["Lipinski_Veber_Score"] = df.apply(lipinski_veber_score, axis=1)
    df["Selection_Score"] = df.apply(compute_selection_score, axis=1)

    # Selection flags
    df["Similarity_Window"] = np.where(
        (df["Max_Known_JAK2_Similarity"] >= args.min_jak2_sim) &
        (df["Max_Known_JAK2_Similarity"] < args.max_jak2_sim),
        "Moderate_JAK2_like", "Outside_window"
    )
    df["Novelty_Class"] = np.where(
        df["Max_Known_JAK2_Similarity"] >= 0.85, "Too similar to known JAK2 seed",
        np.where(df["Max_Known_JAK2_Similarity"] >= 0.65, "High similarity; check novelty",
                 np.where(df["Max_Known_JAK2_Similarity"] >= 0.35, "Moderate similarity; preferred novelty window",
                          "Low JAK2 similarity; exploratory"))
    )

    # Base selected pool
    select_mask = (
        (~df["Exact_Seed_Match"]) &
        (df["Max_Known_JAK2_Similarity"] >= args.min_jak2_sim) &
        (df["Max_Known_JAK2_Similarity"] < args.max_jak2_sim)
    )
    if not args.allow_pains:
        select_mask = select_mask & (~df["PAINS_Flag"])

    selected = df.loc[select_mask].sort_values("Selection_Score", ascending=False).head(args.select_n).copy()

    # If selection is too small, relax with CMap similarity support
    if len(selected) < min(20, args.select_n):
        relaxed_mask = (
            (~df["Exact_Seed_Match"]) &
            (df["Max_Known_JAK2_Similarity"] < args.max_jak2_sim) &
            ((df["Max_Known_JAK2_Similarity"] >= 0.20) | (df["Max_Nonclassical_CMap_Similarity"] >= 0.30))
        )
        if not args.allow_pains:
            relaxed_mask = relaxed_mask & (~df["PAINS_Flag"])
        selected = df.loc[relaxed_mask].sort_values("Selection_Score", ascending=False).head(args.select_n).copy()

    ranked = df.sort_values("Selection_Score", ascending=False).reset_index(drop=True)
    ranked["Rank"] = np.arange(1, len(ranked) + 1)
    selected = selected.sort_values("Selection_Score", ascending=False).reset_index(drop=True)
    selected["Docking_Set_Rank"] = np.arange(1, len(selected) + 1)

    # Save tables
    filtered_cols_first = [
        "Library_ID", "Library_Name", "Library_Source", "Canonical_SMILES",
        "Best_Seed", "Best_Seed_Group", "Max_Seed_Similarity",
        "Max_Known_JAK2_Similarity", "Max_Nonclassical_CMap_Similarity",
        "Selection_Score", "Similarity_Window", "Novelty_Class", "Exact_Seed_Match", "PAINS_Flag",
        "MW", "LogP", "HBD", "HBA", "TPSA", "RotB", "Rings", "HeavyAtoms", "QED", "Lipinski_Veber_Score"
    ]
    cols = [c for c in filtered_cols_first if c in ranked.columns] + [c for c in ranked.columns if c not in filtered_cols_first]
    ranked_top = ranked[cols].head(args.top_n)
    df[cols].to_csv(outdir / "02_external_library_standardized_filtered.csv", index=False)
    ranked_top.to_csv(outdir / "03_similarity_ranked_all.csv", index=False)
    selected[cols + ["Docking_Set_Rank"] if "Docking_Set_Rank" not in cols else cols].to_csv(outdir / "04_selected_external_ligands_for_docking.csv", index=False)

    # Write SDF for selected compounds
    sdf_path = outdir / "05_selected_external_ligands_for_docking.sdf"
    writer = Chem.SDWriter(str(sdf_path))
    for _, row in selected.iterrows():
        mol, can = standardize_smiles(row["Canonical_SMILES"])
        if mol is None:
            continue
        mol = Chem.AddHs(mol)
        if args.make_3d:
            try:
                AllChem.EmbedMolecule(mol, randomSeed=42, maxAttempts=50)
                AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
            except Exception:
                pass
        mol.SetProp("_Name", str(row["Library_Name"]))
        for prop in ["Library_ID", "Library_Source", "Selection_Score", "Max_Known_JAK2_Similarity",
                     "Max_Nonclassical_CMap_Similarity", "Best_Seed", "Novelty_Class", "Canonical_SMILES"]:
            if prop in row:
                mol.SetProp(prop, str(row[prop]))
        writer.write(mol)
    writer.close()

    # Plots
    plt.figure(figsize=(7.5, 5.0), dpi=200)
    plt.hist(ranked["Max_Known_JAK2_Similarity"].dropna(), bins=40, alpha=0.8, label="Known JAK2 seed similarity")
    plt.hist(ranked["Max_Nonclassical_CMap_Similarity"].dropna(), bins=40, alpha=0.6, label="Non-classical CMap seed similarity")
    plt.axvline(args.min_jak2_sim, linestyle="--", linewidth=1)
    plt.axvline(args.max_jak2_sim, linestyle="--", linewidth=1)
    plt.xlabel("Tanimoto similarity")
    plt.ylabel("Compound count")
    plt.title("External library similarity distribution")
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(outdir / "06_similarity_distribution.png", dpi=600)
    plt.close()

    plt.figure(figsize=(7.0, 5.5), dpi=200)
    plt.scatter(ranked["Max_Known_JAK2_Similarity"], ranked["Max_Nonclassical_CMap_Similarity"], s=16, alpha=0.35)
    if len(selected) > 0:
        plt.scatter(selected["Max_Known_JAK2_Similarity"], selected["Max_Nonclassical_CMap_Similarity"], s=28, alpha=0.9, marker="o", label="Selected for docking")
    plt.axvline(args.min_jak2_sim, linestyle="--", linewidth=1)
    plt.axvline(args.max_jak2_sim, linestyle="--", linewidth=1)
    plt.xlabel("Max similarity to known JAK2 seeds")
    plt.ylabel("Max similarity to non-classical CMap seeds")
    plt.title("Similarity space for external-library selection")
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(outdir / "07_similarity_source_scatter.png", dpi=600)
    plt.close()

    note = f"""
No-manual-curation external similarity screening method note
===========================================================

Input library: {args.library}
Valid standardized external molecules after broad filtering: {len(df)}
Selected molecules for docking: {len(selected)}

Default seed design:
- Known/reference JAK2-pocket seeds: Momelotinib, TG-101348, CHEMBL3699463.
- Non-classical CMap/TB-signature seeds: Curcumin, Quercetin.

Selection logic:
- Exact seed matches were excluded.
- Molecules too similar to known JAK2 seeds were penalized/excluded using max known-JAK2 Tanimoto similarity cutoff < {args.max_jak2_sim}.
- Preferred external hits were moderately similar to known JAK2 seeds, using max known-JAK2 Tanimoto similarity >= {args.min_jak2_sim}.
- Additional ranking used similarity to non-classical CMap seeds, QED, and soft Lipinski/Veber drug-like score.
- PAINS-flagged compounds were {'allowed' if args.allow_pains else 'excluded from selected docking set'}.

Recommended use:
- Use 05_selected_external_ligands_for_docking.sdf for LigPrep/Maestro docking.
- After docking and MM-GBSA, select one best external hit for MD together with Momelotinib/native control and one CMap comparator.
- Do not claim external hits are validated JAK2 inhibitors; describe them as external-library predicted candidates selected by moderate similarity, drug-like filtering, docking, and MM-GBSA.
"""
    with open(outdir / "08_method_note_no_manual_external_similarity.txt", "w", encoding="utf-8") as f:
        f.write(textwrap.dedent(note).strip() + "\n")

    print("\nDone.")
    print(f"Output directory: {outdir.resolve()}")
    print(f"Selected for docking: {len(selected)}")
    print(f"Main docking SDF: {sdf_path}")


if __name__ == "__main__":
    main()
