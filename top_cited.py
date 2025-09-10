# top_cited_bulk.py
# PubMed (PMIDs -> DOI) + OpenAlex (citations) => Top N most cited per query
# Python 3.10+

import time
import csv
import warnings
from urllib.parse import quote
from pathlib import Path
from typing import List, Dict, Tuple, Optional

import pandas as pd
import requests
from requests.adapters import HTTPAdapter, Retry
from Bio import Entrez

# ===================== USER CONFIG =====================
# OBAVEZNO postavi tvoj mejl (NCBI pravilo) i mailto za OpenAlex
Entrez.email = "marko.zivanovic@uni.kg.ac.rs"
OPENALEX_MAILTO = "marko.zivanovic@uni.kg.ac.rs"

# Godišnji opseg
YEAR_START = 2010
YEAR_END   = 2025

# Koliko najcitiranijih radova po upitu vadimo
TOP_N = 20

# Gde da snimimo rezultate
OUTDIR = r"D:\1 Rad\BioINGLab\AI OCTE Chapter\Literature\top_cited"
COMBINED_XLSX = "top_cited_all_queries.xlsx"  # Excel sa više sheet-ova

# Pauza između poziva (budimo gentle prema API-jima)
SLEEP_SEC = 0.34

# ====== Sinonim-blokovi (Title/Abstract) ======
def TA_OR(terms):
    return "(" + " OR ".join([f'"{t}"[Title/Abstract]' for t in terms]) + ")"

OCTE_TA   = TA_OR(["osteochondral tissue engineering","osteochondral repair","osteochondral regeneration","osteochondral defect","osteochondral unit","cartilage-bone interface"])
SCAF_TA   = TA_OR(["scaffold","biomaterial scaffold","porous scaffold","composite scaffold","3D scaffold","3D printed scaffold","scaffold design"])
COMP_TA   = TA_OR(["finite element","finite element analysis","computational modelling","computational modeling","numerical simulation","multiscale modelling","multiscale modeling","in silico modelling","in silico modeling"])
AI_TA     = TA_OR(["artificial intelligence","machine learning","deep learning","neural network","predictive modelling","predictive modeling","data-driven"])

PUBTYPES  = '(review[Publication Type] OR journal article[Publication Type])'

def Q(*blocks):
    return "(" + " AND ".join(blocks) + f") AND {PUBTYPES}"

# Lista upita: naziv → PubMed Boolean (Title/Abstract)
QUERIES = {
    # 1) Osnovni domeni
    "octe":               Q(OCTE_TA),
    "scaffolds":          Q(SCAF_TA),
    "comp_model":         Q(COMP_TA),
    "ai":                 Q(AI_TA),

    # 2) Dvojne kombinacije
    "ai_scaffolds":       Q(AI_TA, SCAF_TA),
    "ai_octe":            Q(AI_TA, OCTE_TA),
    "scaffolds_octe":     Q(SCAF_TA, OCTE_TA),
    "comp_scaffolds":     Q(COMP_TA, SCAF_TA),
    "ai_comp":            Q(AI_TA, COMP_TA),

    # 3) Trostruke kombinacije
    "ai_scaffolds_octe":  Q(AI_TA, SCAF_TA, OCTE_TA),
    "ai_comp_scaffolds":  Q(AI_TA, COMP_TA, SCAF_TA),
    "comp_scaffolds_octe":Q(COMP_TA, SCAF_TA, OCTE_TA),
    "ai_comp_octe":       Q(AI_TA, COMP_TA, OCTE_TA),

    # 4) Četvorostruka kombinacija
    "ai_comp_scaffolds_octe": Q(AI_TA, COMP_TA, SCAF_TA, OCTE_TA),
}
# ======================================================


# ================ TEHNIČKE POSTAVKE ===================
warnings.filterwarnings("ignore", category=UserWarning, module="ssl")
OPENALEX_BASE = "https://api.openalex.org/works"

def session_with_retry() -> requests.Session:
    s = requests.Session()
    retries = Retry(
        total=6,
        backoff_factor=0.6,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset(["GET"])
    )
    s.mount("https://", HTTPAdapter(max_retries=retries))
    s.headers.update({"User-Agent": f"TopCitedScript/1.0 (mailto:{OPENALEX_MAILTO})"})
    return s
# ======================================================


# ===================== PUBMED ==========================
def pubmed_pmids(query: str, y0: int, y1: int) -> List[str]:
    term = f'({query}) AND ("{y0}"[Date - Publication] : "{y1}"[Date - Publication])'
    h = Entrez.esearch(db="pubmed", term=term, retmax=100000)
    rec = Entrez.read(h); h.close()
    return rec.get("IdList", [])

def pubmed_meta_with_doi(pmids: List[str]) -> List[Dict]:
    results = []
    BATCH=200
    for i in range(0, len(pmids), BATCH):
        chunk = pmids[i:i+BATCH]
        if not chunk:
            continue
        h = Entrez.efetch(db="pubmed", id=",".join(chunk), rettype="medline", retmode="xml")
        rec = Entrez.read(h); h.close()
        for art in rec.get("PubmedArticle", []):
            artd = {"pmid": None, "title": None, "year": None, "journal": None, "doi": None}
            artd["pmid"] = str(art["MedlineCitation"]["PMID"])

            # Title
            try:
                artd["title"] = str(art["MedlineCitation"]["Article"].get("ArticleTitle",""))
            except: pass

            # Year
            year = None
            try:
                dp = art["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"]
                year = dp.get("Year") or dp.get("MedlineDate","")[:4]
            except: pass
            artd["year"] = year

            # Journal
            try:
                artd["journal"] = art["MedlineCitation"]["Article"]["Journal"]["Title"]
            except: pass

            # DOI
            try:
                ids = art["PubmedData"]["ArticleIdList"]
                for aid in ids:
                    if aid.attributes.get("IdType") == "doi":
                        artd["doi"] = str(aid)
                        break
            except: pass

            results.append(artd)
        time.sleep(SLEEP_SEC)
    return results
# ======================================================


# ===================== OPENALEX ========================
def openalex_bulk_citations(dois: List[str], sess: Optional[requests.Session] = None,
                            sleep_sec: float = SLEEP_SEC) -> Dict[str, Tuple[int, Optional[str]]]:
    if sess is None:
        sess = session_with_retry()

    out: Dict[str, Tuple[int, Optional[str]]] = {}
    BATCH = 180  # bezbedno ispod limita 200
    for i in range(0, len(dois), BATCH):
        chunk = [d.lower() for d in dois[i:i+BATCH] if d]
        if not chunk:
            continue

        filt = "filter=" + ",".join([f"doi:{quote(d)}" for d in chunk])
        url = f"{OPENALEX_BASE}?{filt}&per-page=200&mailto={quote(OPENALEX_MAILTO)}"

        r = sess.get(url, timeout=30)
        if r.status_code == 200:
            j = r.json()
            for item in j.get("results", []):
                doi = (item.get("doi") or "").lower()
                out[doi] = (item.get("cited_by_count", 0), item.get("id"))
        else:
            # fallback pojedinačno
            for d in chunk:
                u = f"{OPENALEX_BASE}/doi:{quote(d)}?mailto={quote(OPENALEX_MAILTO)}"
                rr = sess.get(u, timeout=30)
                if rr.status_code == 200:
                    jj = rr.json()
                    out[d] = (jj.get("cited_by_count", 0), jj.get("id"))
                time.sleep(0.2)

        time.sleep(sleep_sec)
    return out
# ======================================================


# ===================== PIPELINE ========================
def process_query(name: str, query: str, y0: int, y1: int, top_n: int,
                  outdir: Path) -> Path:
    print(f"\n=== QUERY: {name} ===")
    print("🔎 PubMed search…")
    pmids = pubmed_pmids(query, y0, y1)
    print(f"PMIDs: {len(pmids)}")

    print("📥 Fetch metadata + DOIs…")
    arts = pubmed_meta_with_doi(pmids)

    # Ako nema ničega – snimi prazan CSV sa standardnim kolonama i izađi
    if not arts:
        out_csv = outdir / f"top_cited_{name}.csv"
        pd.DataFrame(columns=[
            "rank","citations","title","year","journal","doi","pmid","openalex_id"
        ]).to_csv(out_csv, index=False, encoding="utf-8")
        print(f"⚠️  No records for query '{name}'. Saved empty CSV: {out_csv}")
        return out_csv

    # DOI-jevi
    dois = sorted({(a.get("doi") or "").strip().lower() for a in arts if a.get("doi")})
    print(f"DOIs with value: {len(dois)}")

    print("🔗 Query OpenAlex citations (bulk)…")
    sess = session_with_retry()
    doi2cit = openalex_bulk_citations(dois, sess=sess)

    # Sastavi redove (uvek kreiramo kolone, čak i ako budu prazne)
    out_rows = []
    for a in arts:
        d = (a.get("doi") or "").strip().lower()
        ccount, oa_id = doi2cit.get(d, (0, None))
        out_rows.append({
            "title": a.get("title",""),
            "year": a.get("year",""),
            "journal": a.get("journal",""),
            "doi": a.get("doi",""),
            "pmid": a.get("pmid",""),
            "openalex_id": oa_id or "",
            "citations": ccount
        })

    df = pd.DataFrame(out_rows, columns=[
        "title","year","journal","doi","pmid","openalex_id","citations"
    ])

    if df.empty:
        out_csv = outdir / f"top_cited_{name}.csv"
        pd.DataFrame(columns=[
            "rank","citations","title","year","journal","doi","pmid","openalex_id"
        ]).to_csv(out_csv, index=False, encoding="utf-8")
        print(f"⚠️  No DOI-mapped records for '{name}'. Saved empty CSV: {out_csv}")
        return out_csv

    df["citations"] = pd.to_numeric(df["citations"], errors="coerce").fillna(0).astype(int)
    # sekundarni sort po godini (stariji gore kad je isti broj citata)
    df = df.sort_values(["citations","year"], ascending=[False, True])

    top = df.head(top_n).copy()
    top.insert(0, "rank", range(1, len(top)+1))

    out_csv = outdir / f"top_cited_{name}.csv"
    top.to_csv(out_csv, index=False, encoding="utf-8")
    print(f"✅ Saved: {out_csv}")
    return out_csv
# ======================================================


if __name__ == "__main__":
    outdir = Path(OUTDIR)
    outdir.mkdir(parents=True, exist_ok=True)

    csv_paths = []
    for name, q in QUERIES.items():
        p = process_query(name, q, YEAR_START, YEAR_END, TOP_N, outdir)
        csv_paths.append(p)

    # Excel sa više sheet-ova
    xlsx_path = outdir / COMBINED_XLSX
    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
        for p in csv_paths:
            df = pd.read_csv(p, encoding="utf-8")
            sheet = Path(p).stem.replace("top_cited_", "")[:31]
            df.to_excel(writer, index=False, sheet_name=sheet)
    print(f"📗 Combined Excel saved: {xlsx_path}")
