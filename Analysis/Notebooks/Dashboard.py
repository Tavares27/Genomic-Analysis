# dashboard.py
import streamlit as st
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np

# Carrega o genoma
fasta_path = "Analysis/Data/ecoli.fasta"
record = SeqIO.read(fasta_path, "fasta")
sequence = str(record.seq)

st.title("🔬 E. coli Genomic Dashboard")
st.write(f"**ID:** {record.id}")
st.write(f"**Description:** {record.description}")
st.write(f"**Lenght:** {len(record.seq):,} bases")

# GC e AT content globais
gc_count = sequence.count("G") + sequence.count("C")
gc_percent = gc_count / len(sequence) * 100
at_count = sequence.count("A") + sequence.count("T")
at_percent = at_count / len(sequence) * 100

st.metric("GC (%)", f"{gc_percent:.2f}")
st.metric("AT (%)", f"{at_percent:.2f}")

# Gráfico de GC content ao longo do genoma
window_size = st.slider("Window lenght (bases)", 100, 5000, 1000, 100)
gc_content = []
for i in range(0, len(sequence) - window_size + 1, window_size):
    window = sequence[i:i+window_size]
    gc = (window.count("G") + window.count("C")) / window_size * 100
    gc_content.append(gc)
positions = np.arange(0, len(gc_content)) * window_size

fig, ax = plt.subplots(figsize=(10,4))
ax.plot(positions, gc_content, color='green')
ax.set_xlabel("Inicial Position (bases)")
ax.set_ylabel("GC (%)")
ax.set_title("GC content through the genome")
st.pyplot(fig)

# Contagem de bases
base_counts = {base: sequence.count(base) for base in "ATGC"}
st.write("### Nitrogenous Bases Count")
st.write(base_counts)