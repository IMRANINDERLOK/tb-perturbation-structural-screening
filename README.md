# Perturbation-guided host-response analysis identifies JAK2-focused candidate compounds for active tuberculosis

🧬 **Host-response signature | Perturbation screening | CMap/CLUE | JAK2 docking | Molecular dynamics**

This repository contains processed datasets and analysis outputs associated with the manuscript:

**“Perturbation-guided host-response analysis identifies JAK2-focused candidate compounds for active tuberculosis”**

**Graphical Abstract**: <img width="1491" height="1055" alt="Graphical Abstract" src="https://github.com/user-attachments/assets/3ead53a9-41cf-4e5b-8d56-f18490486b92" />


## 🔎 Study focus

Active tuberculosis blood transcriptional signatures are commonly used to describe disease state, treatment response and diagnostic classification. In this study, an **ATB_pre-to-ATB_12m directional host-response signature** was used as a functional input to identify perturbagens predicted to oppose the active-disease transcriptional pattern.

The analysis connected active-state signature reversal with candidate curation, CMap/CLUE connectivity, target-axis mapping and JAK2-focused structural assessment. The final structure-based evaluation supported **EN84** and **ZINC000149253998** as computationally prioritised compounds for follow-up testing in tuberculosis-relevant host-cell models.

## 🧪 Analysis summary

The study used the ATB_pre-to-ATB_12m host-response signature to identify compound perturbation profiles predicted to reverse the active tuberculosis expression state. Directional perturbation analysis showed that suppression of the active-state gene set was the dominant matching signal.

Candidate compounds were then refined using CMap/CLUE connectivity and target-axis mapping. The prioritised candidate space was mainly linked to kinase, inflammatory and interferon-response-related axes. Because the active-state signature was centred on interferon-associated genes, STAT1 and GBP1 were used as biological anchors. JAK2 was selected as a structure-testable kinase node connected to the STAT1-related signalling route.

The JAK2-focused analysis used the 8BXH structure, native momelotinib redocking, molecular docking, MM-GBSA rescoring, external-library screening and 300 ns molecular dynamics simulation summaries.

## 🧬 Input signature

The host-response signature used in this study was generated from the publicly available Gene Expression Omnibus dataset:

**GSE40553**

Processed datasets from the previous signature-generation study are available at:

**Previous study repository:** `insert_previous_repository_link_here`

## 💊 JAK2 bioactivity data

JAK2 bioactivity data were obtained from the ChEMBL database using the target identifier:

**CHEMBL2971**

The ChEMBL-derived dataset was processed to retain JAK2 small-molecule bioactivity records relevant for downstream compound comparison and modelling.

## 🧩 Key outputs

This repository provides processed files supporting:

- perturbation terms predicted to oppose the active-TB host-response signature;
- curated candidate-level compound tables;
- CMap/CLUE connectivity classes;
- kinase, inflammatory and interferon-response target-axis mapping;
- JAK2 docking and MM-GBSA comparisons;
- external-library compound screening outputs;
- molecular dynamics summary data for JAK2–momelotinib, JAK2–EN84 and JAK2–ZINC000149253998 complexes.

## 🎯 Prioritised compounds

The current analysis prioritised:

- **EN84**
- **ZINC000149253998**

These compounds showed JAK2 ATP-pocket engagement and retained binding-site occupancy during 300 ns molecular dynamics simulations. These findings are computational and require experimental testing in relevant mycobacterial host-cell models.

## 📊 Data availability

Processed datasets from the previous signature-generation study are available in the GitHub repository:

**Previous study repository:** `https://github.com/IMRANINDERLOK/PTLD-ML-omics.git`

## ⚠️ Note

This repository provides computational analysis outputs only. The prioritised compounds have not been experimentally validated for anti-tuberculosis or host-directed therapeutic activity in this study.

## 📬 Contact

For questions about the data or analysis files, please contact:

**Mohd Imran**  
**Email:** imran.pchem@gmail.com
