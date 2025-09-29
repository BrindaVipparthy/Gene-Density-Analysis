# STEP 1: CLEAN FILE
import pandas as pd # use Pandas to clean the file and calculate gene density


# Set the filename
filename = "human_genome.txt"

# Define the expected column names that are in the file
columns = [
    "#bin", "name", "chrom", "strand", "txStart", "txEnd",
    "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds",
    "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames"
]
# Chromosome sizes (GRCh38/hg38 reference genome) These are fixed and according to the reference genome
chrom_sizes = {
    'chr1': 248956422,
    'chr2': 242193529,
    'chr3': 198295559,
    'chr4': 190214555,
    'chr5': 181538259,
    'chr6': 170805979,
    'chr7': 159345973,
    'chr8': 145138636,
    'chr9': 138394717,
    'chr10': 133797422,
    'chr11': 135086622,
    'chr12': 133275309,
    'chr13': 114364328,
    'chr14': 107043718,
    'chr15': 101991189,
    'chr16': 90338345,
    'chr17': 83257441,
    'chr18': 80373285,
    'chr19': 58617616,
    'chr20': 64444167,
    'chr21': 46709983,
    'chr22': 50818468,
    'chrX': 156040895,
    'chrY': 57227415,
    'chrM': 16569
}


# Load the file AS DataFrame and read
df = pd.read_csv(filename, sep="\t", names=columns, header=0, low_memory=False)

# Drop unnecessary columns, we only want the columns necessary for analysis purposes
df.drop(columns=["#bin", "score", "cdsStartStat", "cdsEndStat", "exonFrames"], inplace=True)

# Preview the cleaned data to make sure we have necessary columns
print("Cleaned Data:")
print(df.head())
print("Updated Columns:")
print(df.columns)

# STEP 2: GENE DENSITY ANALYSIS
# first calculate gene lengths = difference bet. transcript end and start
df["gene_length"] = df["txEnd"] - df["txStart"]

# Keep only standard chromosomes from the main chromosomes
main_chroms = set(chrom_sizes.keys())
df = df[df["chrom"].isin(main_chroms)] # refer to the chromosome dict

# Group by chromosome and calculate unique genes, avg gene length, and sum of gene lengths
grouped = df.groupby("chrom").agg(
    gene_count=("name2", "nunique"),      # number of unique gene names present in chr
    avg_gene_length=("gene_length", "mean"),
    total_gene_length=("gene_length", "sum")
).reset_index()

# Preview the result from above so we can see the stats for each chr
print(grouped.head())

# Add chromosome size
grouped["chrom_size"] = grouped["chrom"].map(chrom_sizes)

# Drop rows where chrom_size is missing (non-standard chromosomes)
grouped = grouped.dropna(subset=["chrom_size"])

# Calculate gene density per megabase
grouped["gene_density"] = grouped["gene_count"] / (grouped["chrom_size"] / 1e6)

# Final preview with info about gene density
print(grouped[["chrom", "gene_count", "chrom_size", "gene_density"]].sort_values("gene_density", ascending=False).head())

# convert file to csv so we can analyze in Tableau
df.to_csv("gene_density_summary.csv", index=False)

# STEP 3: PCA AND CLUSTERING
# import all libraries to process the data
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import AgglomerativeClustering
import matplotlib.pyplot as plt
import seaborn as sns

# Select the features for PCA
features = grouped[["gene_count", "chrom_size", "gene_density"]]

# Standardize the features
scaler = StandardScaler()
scaled_features = scaler.fit_transform(features)

# Perform PCA reduce to 2 components for visualization
pca = PCA(n_components=2)
pca_components = pca.fit_transform(scaled_features)

# Add PCA results to the grouped dataframe
grouped["PC1"] = pca_components[:, 0]
grouped["PC2"] = pca_components[:, 1]

# Perform clustering and group chromosomes into 3 clusters
clustering = AgglomerativeClustering(n_clusters=3)
grouped["cluster"] = clustering.fit_predict(scaled_features)

# Plot the PCA + clustering result
plt.figure(figsize=(10, 6))
sns.scatterplot(data=grouped, x="PC1", y="PC2", hue="cluster", palette="Set2", s=100)

# Label each point with chromosome name
for i in range(len(grouped)):
    plt.text(grouped["PC1"][i] + 0.05, grouped["PC2"][i], grouped["chrom"][i], fontsize=8)

plt.title("Chromosome Clustering via PCA")
plt.xlabel("Principal Component 1")
plt.ylabel("Principal Component 2")
plt.grid(True)
plt.tight_layout()
plt.show()

# convert to csv so we can visualize in tableu
# Save PCA and clustering results to CSV
grouped.to_csv("gene_density_pca_clusters.csv", index=False)
