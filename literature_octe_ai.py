from Bio import Entrez
from collections import defaultdict, OrderedDict
from tqdm import tqdm
import time
import csv
import re

# =============== USER SETTINGS ===============
# Obavezno postavi svoj mejl zbog NCBI pravila
Entrez.email = "your.email@example.com"

# Godine koje analiziramo (OCTE paper scope 2010–2025)
YEARS = list(range(2010, 2026))

# Pauza između poziva da budemo "gentle" prema NCBI
SLEEP_SEC = 0.5

# PubMed vrsta publikacija (pregledi i originalni radovi)
PUBTYPES = '(review[Publication Type] OR journal article[Publication Type])'

# Folder za rezultate (Windows raw string zbog backslash)
OUTDIR = r"D:\1 Rad\BioINGLab\AI OCTE Chapter\Literature"
# ============================================


# ----------------- HELPERI ------------------
def ta(x):
    """Wrap term list for [Title/Abstract] OR spoj."""
    if isinstance(x, (list, tuple, set)):
        parts = [f'"{t}"[Title/Abstract]' for t in x]
        return "(" + " OR ".join(parts) + ")"
    return f'("{x}"[Title/Abstract])'

def mk_query(blocks, pubtypes=PUBTYPES):
    """
    blocks: lista stringova (već grupisanih sa zagradama) koje se AND-uju
    """
    blocks = [b.strip() for b in blocks if b and b.strip()]
    core = " AND ".join(blocks)
    if pubtypes:
        return f"({core}) AND {pubtypes}"
    return f"({core})"

def year_filter(y):
    return f'("{y}"[Date - Publication])'

def run_count_by_year(term):
    """Vrati dict {year: count} za dati PubMed upit."""
    counts = {}
    for y in tqdm(YEARS, desc="Years"):
        q = f'({term}) AND {year_filter(y)}'
        handle = Entrez.esearch(db="pubmed", term=q, retmax=0)
        record = Entrez.read(handle)
        handle.close()
        counts[y] = int(record.get("Count", "0"))
        time.sleep(SLEEP_SEC)
    return counts

def save_year_series_csv(name, counts):
    fn = rf"{OUTDIR}\{name}.csv"
    with open(fn, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["Year", "Count"])
        for y in YEARS:
            w.writerow([y, counts.get(y, 0)])
    print(f"✅ Saved: {fn}")

# --------------- DOMENE & SINONIMI ----------
# 1) Osteochondral tissue engineering (OCTE)
octe_terms = [
    "osteochondral tissue engineering",
    "osteochondral repair",
    "osteochondral regeneration",
    "osteochondral defect",
    "osteochondral unit",
    "cartilage-bone interface"
]

# 2) Scaffolds
scaffold_terms = [
    "scaffold",
    "biomaterial scaffold",
    "porous scaffold",
    "composite scaffold",
    "3D scaffold",
    "3D printed scaffold",
    "scaffold design"
]

# 3) Computational modelling (klasično)
compmod_terms = [
    "finite element",
    "finite element analysis",
    "computational modelling",
    "computational modeling",  # US varijanta
    "numerical simulation",
    "multiscale modelling",
    "multiscale modeling",
    "in silico modelling",
    "in silico modeling"
]

# 4) Artificial intelligence (AI)
ai_terms = [
    "artificial intelligence",
    "machine learning",
    "deep learning",
    "neural network",
    "predictive modelling",
    "predictive modeling",
    "data-driven"
]

# Pretvaramo termine u ( ... OR ... ) blokove za [Title/Abstract]
OCTE = ta(octe_terms)
SCAFF = ta(scaffold_terms)
COMPMOD = ta(compmod_terms)
AI = ta(ai_terms)

# ------------- DEFINICIJA UPITA -------------
# Jednostruki domeni
queries = OrderedDict()
queries["octe"]         = mk_query([OCTE])
queries["scaffolds"]    = mk_query([SCAFF])
queries["comp_model"]   = mk_query([COMPMOD])
queries["ai"]           = mk_query([AI])

# Dvojne kombinacije
queries["ai_scaffolds"]   = mk_query([AI, SCAFF])
queries["ai_octe"]        = mk_query([AI, OCTE])
queries["scaffolds_octe"] = mk_query([SCAFF, OCTE])
queries["comp_scaffolds"] = mk_query([COMPMOD, SCAFF])
queries["ai_comp"]        = mk_query([AI, COMPMOD])

# Trojne kombinacije
queries["ai_scaffolds_octe"]   = mk_query([AI, SCAFF, OCTE])
queries["ai_comp_scaffolds"]   = mk_query([AI, COMPMOD, SCAFF])
queries["comp_scaffolds_octe"] = mk_query([COMPMOD, SCAFF, OCTE])
queries["ai_comp_octe"]        = mk_query([AI, COMPMOD, OCTE])

# Četvorostruki preseci
queries["ai_comp_scaffolds_octe"] = mk_query([AI, COMPMOD, SCAFF, OCTE])

# ============== IZVRSAVANJE =================
if __name__ == "__main__":
    # 1) Snimi pojedinačne CSV serije po upitu
    all_counts = {}
    for name, term in queries.items():
        print(f"\n🔎 Running query: {name}")
        counts = run_count_by_year(term)
        all_counts[name] = counts
        save_year_series_csv(name, counts)

    # 2) Napravi "matrix" CSV (godine × svi upiti) — zgodno za heatmap
    matrix_fn = rf"{OUTDIR}\octe_ai_lit_matrix.csv"
    with open(matrix_fn, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        header = ["Year"] + list(queries.keys())
        w.writerow(header)
        for y in YEARS:
            row = [y] + [all_counts[q].get(y, 0) for q in queries.keys()]
            w.writerow(row)
    print(f"🧩 Saved matrix: {matrix_fn}")

    print("\n🎉 All queries processed successfully!")
