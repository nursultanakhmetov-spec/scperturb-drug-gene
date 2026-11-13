import pandas as pd

with open("details_output.txt", "w", encoding="utf-8") as out:
    def p(*args):
        print(*args); print(*args, file=out)

    # ---------- Nadig ----------
    o = pd.read_csv("nadig_obs.csv.gz", index_col=0, low_memory=False)
    p("=== NADIG ===")
    p("метка контроля:", o.loc[o["nperts"] == 0, "perturbation"].unique()[:5])
    vc = o.loc[o["nperts"] == 1, "perturbation"].value_counts()
    p("генов-мишеней:", len(vc))
    p("клеток на ген:", vc.describe().round(1).to_dict())
    p("генов с >=50 клеток:", (vc >= 50).sum(), " с >=100:", (vc >= 100).sum(), " с >=200:", (vc >= 200).sum())
    vc.to_csv("nadig_perturbations.csv")

    # ---------- Norman ----------
    o = pd.read_csv("norman_obs.csv.gz", index_col=0, low_memory=False)
    p("\n=== NORMAN ===")
    per = o.drop_duplicates("perturbation")
    p("одиночных:", (per["nperts"] == 1).sum(), " парных:", (per["nperts"] == 2).sum())
    good = o[(o["good_coverage"] == True) & (o["number_of_cells"] == 1)]
    p("клеток после фильтра качества:", len(good), "из", len(o))
    vc = good["perturbation"].value_counts()
    p("клеток на возмущение:", vc.describe().round(1).to_dict())
    p("возмущений с >=50:", (vc >= 50).sum(), " с >=200:", (vc >= 200).sum())
    singles = sorted(per.loc[per["nperts"] == 1, "perturbation"])
    p("одиночные гены:", singles)
    vc.to_csv("norman_perturbations.csv")

    # ---------- sciplex2 ----------
    o = pd.read_csv("sciplex2_obs.csv.gz", index_col=0, low_memory=False)
    p("\n=== SCIPLEX2: препарат x доза ===")
    p(pd.crosstab(o["perturbation"], o["dose_value"]).to_string())
    p("клеток без метки лунки (hash_umis пусто):", o["hash_umis"].isna().sum())
    p("их возмущения:", o.loc[o["hash_umis"].isna(), "perturbation"].value_counts().to_dict())

    # ---------- sciplex3 ----------
    o = pd.read_csv("sciplex3_obs.csv.gz", index_col=0, low_memory=False)
    p("\n=== SCIPLEX3 ===")
    p("линия x время:"); p(pd.crosstab(o["cell_line"].fillna("нет"), o["time"].fillna(-1)).to_string())
    p("линия x доза:");  p(pd.crosstab(o["cell_line"].fillna("нет"), o["dose_value"]).to_string())
    k = o[o["cell_line"] == "K562"]
    p("K562: препаратов", k["perturbation"].nunique() - 1, " контрольных клеток", (k["perturbation"] == "control").sum())
    vc = k.groupby(["perturbation", "dose_value"]).size()
    p("K562, клеток на пару препарат-доза:", vc.describe().round(1).to_dict())
    tab = o.drop_duplicates("perturbation")[["perturbation", "target", "pathway", "pathway_level_1", "pathway_level_2"]]
    tab.sort_values("target").to_csv("sciplex3_drug_targets.csv", index=False)
    p("все классы мишеней:", sorted(o["target"].dropna().unique()))