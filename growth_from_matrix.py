import csv
import math
from pathlib import Path

import pandas as pd

# Ako želiš bez seaborn-a, ostavi samo matplotlib:
import matplotlib.pyplot as plt
USE_SEABORN = True
try:
    import seaborn as sns
except ImportError:
    USE_SEABORN = False

# ==== USER CONFIG ===========================================================
MATRIX_CSV = r"D:\1 Rad\BioINGLab\AI OCTE Chapter\Literature\octe_ai_lit_matrix.csv"

# Izaberi koje kolone (upite) iz matrice želiš da prikažeš
CATEGORIES = [
    "octe", "scaffolds", "comp_model", "ai",
    "ai_scaffolds", "ai_octe", "scaffolds_octe",
    "comp_scaffolds", "ai_comp", "ai_scaffolds_octe",
    "ai_comp_scaffolds", "comp_scaffolds_octe",
    "ai_comp_octe", "ai_comp_scaffolds_octe"
]

# Godišnji opseg – možeš ostaviti None pa će se uzeti min/max iz CSV-a
YEAR_START = 2010
YEAR_END   = 2025

# Snimanje figura
OUTDIR = r"D:\1 Rad\BioINGLab\AI OCTE Chapter\Literature"
SAVE_FIGS = True

# (Opcionalno) boje za kategorije (kao u tvom starom primeru)
CUSTOM_PALETTE = {
    # primer: "ai": "#4878CF",
}
# ===========================================================================


def load_matrix(path):
    df = pd.read_csv(path, encoding="utf-8")
    # Osigurajmo da "Year" bude int i sortirano
    df["Year"] = df["Year"].astype(int)
    df = df.sort_values("Year")
    return df

def subset_years(df, y0=None, y1=None):
    if y0 is None: y0 = int(df["Year"].min())
    if y1 is None: y1 = int(df["Year"].max())
    return df[(df["Year"] >= y0) & (df["Year"] <= y1)].copy()

def first_nonzero_pair(series_years, series_values):
    """
    Vrati (y_start, v_start, y_end, v_end) gde je v_start prvi nenulti,
    v_end je poslednji (po vremenu) nenulti. Ako nema dovoljno podataka, vrati None.
    """
    # nađi prvi nenulti
    start_idx = None
    for i, v in enumerate(series_values):
        if v and v > 0:
            start_idx = i
            break
    if start_idx is None:
        return None

    # nađi poslednji nenulti
    end_idx = None
    for i in range(len(series_values)-1, -1, -1):
        if series_values[i] and series_values[i] > 0:
            end_idx = i
            break
    if end_idx is None or end_idx == start_idx:
        return None

    return (series_years[start_idx], series_values[start_idx],
            series_years[end_idx],   series_values[end_idx])

def cagr(p0, pt, t_years):
    if t_years <= 0 or p0 <= 0 or pt <= 0:
        return 0.0
    return (pt / p0) ** (1.0 / t_years) - 1.0

def compute_cagr_for_categories(df, categories):
    years = df["Year"].tolist()
    result = []
    for cat in categories:
        if cat not in df.columns:
            result.append((cat, 0.0))
            continue
        vals = df[cat].fillna(0).astype(int).tolist()
        pair = first_nonzero_pair(years, vals)
        if pair is None:
            result.append((cat, 0.0))
            continue
        y0, v0, y1, v1 = pair
        rate = cagr(v0, v1, y1 - y0) * 100.0
        result.append((cat, rate))
    # sortiraj opadajuće
    result.sort(key=lambda x: x[1], reverse=True)
    return pd.DataFrame(result, columns=["Category", "AverageAnnualGrowthRate(%)"])

def plot_cagr_bar(df_cagr, title="Average Annual Growth Rate (%) by Category",
                  savepath=None):
    fig, ax = plt.subplots(figsize=(10, 5))
    if USE_SEABORN:
        if CUSTOM_PALETTE:
            pal = [CUSTOM_PALETTE.get(cat, None) for cat in df_cagr["Category"]]
            # ako neka boja nije zadana, seaborn će dodeliti podrazumevanu
            sns.barplot(data=df_cagr, y="Category", x="AverageAnnualGrowthRate(%)",
                        palette=pal)
        else:
            sns.barplot(data=df_cagr, y="Category", x="AverageAnnualGrowthRate(%)")
    else:
        # čist matplotlib
        y = range(len(df_cagr))
        x = df_cagr["AverageAnnualGrowthRate(%)"].values
        ax.barh(y, x)
        ax.set_yticks(y)
        ax.set_yticklabels(df_cagr["Category"])

    ax.set_xlabel("Average Annual Growth Rate (%)", fontsize=14)
    ax.set_ylabel("Category", fontsize=14)
    ax.set_title(title, fontsize=16)
    plt.tight_layout()
    if savepath:
        plt.savefig(savepath, dpi=300)
    plt.show()

def plot_lines(df, categories, title="Publications by Year", savepath=None):
    fig, ax = plt.subplots(figsize=(12, 6))
    for cat in categories:
        if cat in df.columns:
            ax.plot(df["Year"], df[cat], label=cat)
    ax.set_xlabel("Year")
    ax.set_ylabel("Publications (count)")
    ax.set_title(title)
    ax.legend(loc="upper left", ncol=2, fontsize=8)
    plt.tight_layout()
    if savepath:
        plt.savefig(savepath, dpi=300)
    plt.show()


if __name__ == "__main__":
    df = load_matrix(MATRIX_CSV)
    df = subset_years(df, YEAR_START, YEAR_END)

    # 1) CAGR tabela
    df_cagr = compute_cagr_for_categories(df, CATEGORIES)
    print(df_cagr.to_string(index=False))

    # 2) Bar plot (CAGR)
    outdir = Path(OUTDIR) if OUTDIR else Path(".")
    outdir.mkdir(parents=True, exist_ok=True)

    bar_png = str(outdir / "cagr_by_category.png") if SAVE_FIGS else None
    plot_cagr_bar(df_cagr, savepath=bar_png)

    # 3) Line plot (opciono)
    line_png = str(outdir / "lines_by_year.png") if SAVE_FIGS else None
    plot_lines(df, CATEGORIES, savepath=line_png)
