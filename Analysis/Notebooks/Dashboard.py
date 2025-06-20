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
st.metric("Adenine", base_counts["A"])
st.metric("Thymine", base_counts["T"])
st.metric("Guanine", base_counts["G"])
st.metric("Cytosine", base_counts["C"])

# Bar chart for nitrogenous base counts
st.write("### Nitrogenous Bases Count (Bar Chart)")

fig2, ax2 = plt.subplots()
ax2.bar(base_counts.keys(), base_counts.values(), color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728'])
ax2.set_xlabel("Base")
ax2.set_ylabel("Count")
ax2.set_title("Nitrogenous Base Counts")
st.pyplot(fig2)

# Contagem de pares de bases
st.write("### Base Pair Counts")
pares = {
    "A-T": record.seq.count("AT") + record.seq.count("TA"),
    "G-C": record.seq.count("GC") + record.seq.count("CG"),
}
st.metric("A-T (Adenine-Thymine)", pares["A-T"])
st.metric("G-C (Guanine-Cytosine)", pares["G-C"])   

# Gráfico de pares de bases
st.write("### Linked Pairs Count (Bar Chart)")
fig3, ax3 = plt.subplots()
ax3.bar(pares.keys(), pares.values(), color=['#9467bd', '#8c564b'])
ax3.set_xlabel("Pair")
ax3.set_ylabel("Count")
ax3.set_title("Linked Pairs Counts")
st.pyplot(fig3)

# Search for regions rich in GC and AT content (fixed windows of 1000 bases)
gc_regions = []
at_regions = []
window = 1000
threshold = 70

for i in range(0, len(sequence) - window + 1, window):
    subseq = sequence[i:i+window]
    gc_rcontent = (subseq.count("G") + subseq.count("C")) / window * 100
    at_rcontent = (subseq.count("A") + subseq.count("T")) / window * 100
    if gc_rcontent > threshold:
        gc_regions.append(i)
    elif at_rcontent > threshold:
        at_regions.append(i)

# Bar chart: number of GC-rich and AT-rich regions

st.metric("GC-rich regions", len(gc_regions))
st.metric("AT-rich regions", len(at_regions))

for i in range(len(gc_regions)):
    st.write(f"GC-rich region {i+1}: Start at base {gc_regions[i]}")
for i in range(len(at_regions)):
    st.write(f"AT-rich region {i+1}: Start at base {at_regions[i]}")
st.write("### GC/AT-rich Regions (Bar Chart)")
region_counts = {"GC-rich": len(gc_regions), "AT-rich": len(at_regions)}
fig4, ax4 = plt.subplots()
ax4.bar(region_counts.keys(), region_counts.values(), color=['#2ca02c', '#ff7f0e'])
ax4.set_ylabel("Number of Regions")
ax4.set_title(f"Regions with >{threshold}% GC or AT (window={window})")
st.pyplot(fig4)

# Pie chart for base composition in a selected region

st.write("### Base Proportion in a Selected Region")

# Número total de janelas
num_windows = (len(sequence) - window) // window + 1

# Seleção da janela (região)
region_idx = st.number_input(
    "Select region number (0 = first window)", 
    min_value=0, 
    max_value=num_windows-1, 
    value=0, 
    step=1
)
region_start = region_idx * window
region_seq = sequence[region_start:region_start+window]

region_base_counts = {base: region_seq.count(base) for base in "ATGC"}

fig5, ax5 = plt.subplots()
ax5.pie(region_base_counts.values(), labels=region_base_counts.keys(), autopct='%1.1f%%', colors=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728'])
ax5.set_title(f"Base Proportion (Region {region_start} - {region_start+window})")
st.pyplot(fig5)

st.write("### Finding specific sequences")

seq = st.text_input("Enter a sequence to search for (e.g., 'ATG')", "ATG").upper()
if seq:
    if not all(base in "ATGC" for base in seq):
        st.error("Invalid sequence! Please enter a sequence containing only A, T, G, and C.")
    else:
        st.success(f"Searching for sequence: {seq}")
        # Encontrar todas as posições da sequência
        positions = []
        idx = sequence.find(seq)
        while idx != -1:
            positions.append(idx)
            idx = sequence.find(seq, idx + 1)
        st.metric(f"Occurrences of '{seq}'", len(positions))
        if positions:
            st.write(f"### Distribution of '{seq}' occurrences per window along the genome")
            # Use o mesmo window_size do gráfico de GC content
            num_windows = (len(sequence) - window_size) // window_size + 1
            occurrences_per_window = []
            for i in range(num_windows):
                start = i * window_size
                end = start + window_size
                window_seq = sequence[start:end]
                count = window_seq.count(seq)
                occurrences_per_window.append(count)
            window_positions = [i * window_size for i in range(num_windows)]

            fig7, ax7 = plt.subplots(figsize=(10, 3))
            ax7.plot(window_positions, occurrences_per_window, marker='o', color='purple')
            ax7.set_xlabel("Window start position (bases)")
            ax7.set_ylabel(f"Occurrences of '{seq}'")
            ax7.set_title(f"Occurrences of '{seq}' per window ({window_size} bases)")
            st.pyplot(fig7)
        else:
            st.info("Sequence not found in the genome.")