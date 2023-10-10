import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram

# Load protein sequences from a CSV file with two columns: Name and Sequence
csv_file = "Protein sequences.csv"
# Read the CSV file into a pandas DataFrame
df = pd.read_csv(csv_file)
# Extract sequence names and sequences from the DataFrame
sequence_names = df["Name"].tolist()
sequences = df["protein_sequence"].tolist()
# Calculate the distance matrix (similarity matrix) based on sequence alignment
# You can use different alignment methods and scoring matrices here
# For example, you can use pairwise2 module from Bio.Align to perform alignments
from Bio import pairwise2
alignment_matrix = np.zeros((len(sequences), len(sequences)))
for i in range(len(sequences)):
    for j in range(i + 1, len(sequences)):
        alignment_score = pairwise2.align.globalxx(sequences[i], sequences[j], one_alignment_only=True)[0].score
        alignment_matrix[i, j] = alignment_score
        alignment_matrix[j, i] = alignment_score
# Perform hierarchical clustering
linkage_matrix = linkage(alignment_matrix, method='average', metric='euclidean')
# Create a dendrogram plot
plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)  # Create a subplot for the dendrogram
dendrogram(linkage_matrix, labels=sequence_names, orientation="top", leaf_font_size=8)
plt.xticks(rotation=90)
plt.xlabel("Genes")
plt.ylabel("Distance")
plt.title("Hierarchical Clustering Dendrogram")
# Create a heatmap plot for the distance matrix
plt.subplot(1, 2, 2)  # Create a subplot for the heatmap
plt.imshow(alignment_matrix, cmap='hot', interpolation='nearest', aspect='auto')
plt.colorbar(label='Alignment Score')
plt.xticks(range(len(sequence_names)), sequence_names, rotation=90)
plt.yticks(range(len(sequence_names)), sequence_names)
plt.xlabel("Genes")
plt.ylabel("Genes")
plt.title("Distance Matrix Heatmap")
plt.tight_layout()
plt.show()





