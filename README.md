# Genomic Analysis - E. coli

This project performs an exploratory analysis of the **Escherichia coli** genome using FASTA and GFF files, with interactive visualization via Streamlit.

## Features

- Reads FASTA (genome) and GFF (annotations) files
- Basic genome statistics (ID, description, length, GC%, AT%)
- Nitrogenous base counts (A, T, G, C)
- Detection of GC-rich and AT-rich regions
- GC content plot across the genome (sliding window)
- Interactive dashboard with Streamlit

## How to Use

### 1. Requirements

- Python 3.8+
- [Biopython](https://biopython.org/)
- [Streamlit](https://streamlit.io/)
- [matplotlib](https://matplotlib.org/)
- [numpy](https://numpy.org/)

Install dependencies with:

```sh
pip install biopython streamlit matplotlib numpy
```

### 2. File Structure

```
Genomic-Analysis/
│
├── Analysis/
│   ├── Data/
│   │   ├── ecoli.fasta
│   │   └── ecoli.gff
│   └── Notebooks/
│       ├── Main.ipynb
│       └── Dashboard.py
└── README.md
```

### 3. Running the Dashboard

In the terminal, run:

```sh
streamlit run Analysis/Notebooks/Dashboard.py
```

The dashboard will open automatically in your browser (http://localhost:8501).

### 4. Customization

- To analyze another genome, replace the `ecoli.fasta` file in `Analysis/Data/`.
- Adjust the window size in the dashboard for different GC content resolutions.

## Example Usage

- View basic genome statistics.
- Visualize the GC content plot across the genome.
- Explore base counts and GC/AT-rich regions.

## Credits

- Project developed by Felipe Tavares.
- Example data: NCBI.

---

Feel free to contribute or adapt for other genomes!