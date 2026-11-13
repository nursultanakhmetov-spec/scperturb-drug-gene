# scPerturb: dataset inspection and a first drug-to-gene comparison

Scripts that inspect four scPerturb datasets and test whether a drug's
expression fingerprint matches the knockdown of its target gene
(sci-Plex 3 K562 drugs vs. Replogle K562 CRISPRi knockdowns).

## Input data (not included, too large for GitHub)

Download into this folder:

From scPerturb (scperturb.org):
- NadigOConner2024_hepg2.h5ad
- NormanWeissman2019_filtered.h5ad
- SrivatsanTrapnell2020_sciplex2.h5ad
- SrivatsanTrapnell2020_sciplex3.h5ad

From Replogle et al. 2022 (Figshare / gwps.wi.mit.edu):
- K562_gwps_normalized_bulk_01.h5ad

## Setup

    pip install anndata scanpy pandas

## Scripts (run in order)

| Script | What it does | Main output |
|---|---|---|
| 01_inspect.py | Reads metadata of the four scPerturb files | inspect_output.txt, *_obs.csv.gz, *_var.csv |
| 02_details.py | Cells per perturbation, doses, cell lines | details_output.txt, sciplex3_drug_targets.csv |
| 03_sciplex2_fingerprints.py | sci-Plex 2 fingerprints, check of known genes | sciplex2_check.txt |
| 04_sciplex3_k562.py | K562 drug fingerprints from sci-Plex 3 | sciplex3_k562_check.txt, sciplex3_k562_meta.csv |
| 05_replogle_inspect.py | Layout of the Replogle file | replogle_inspect.txt |
| 06_compare.py | First comparison on hand-picked drugs | compare_check.txt |
| 07_compare_all.py | Full comparison, all drugs and doses | compare_all_summary.txt, compare_all_table.csv |

Files ending in .h5ad and .csv.gz are excluded by .gitignore.
04_sciplex3_k562.py needs about 8 GB of free memory.

## Main result

For 7 of 126 testable drugs the known target gene ranks at or near the top
of 5,293 knockdowns (e.g. nilotinib and dasatinib -> BCR at rank 1,
temsirolimus -> MTOR, flavopiridol -> CDK9, GSK-LSD1 -> KDM1A).
See compare_all_summary.txt for details and limitations.