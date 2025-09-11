[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17093170.svg)](https://doi.org/10.5281/zenodo.17093170)
# octe-ai-mining-framework

This repository provides a standardized and reproducible framework for literature mining, citation analysis, and visualization in the context of AI, scaffolds, computational modeling, and OCTE.

## Contents
- `literature_octe_ai.py` — Automated PubMed mining (NCBI Entrez API).  
- `top_cited.py` — Extraction of most cited papers (via OpenAlex API).  
- `growth_from_matrix.py` — Visualization of trends (growth curves, heatmaps, pie charts).  

## Installation
Requires **Python 3.9+** and the following packages:
```bash
pip install pandas matplotlib requests
## Citation
If you use this framework in your research, please cite it as:

Živanović, M. (2025). *octe-ai-mining-framework (v1.0.1)* [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.17093171

