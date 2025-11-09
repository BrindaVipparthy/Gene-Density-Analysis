# Gene Density and Chromosome Clustering Analysis

## Overview:

This study investigates gene distribution across the human genome and analyzes gene
density and gene count within each chromosome. The raw data was retrieved from the UCSC
Genome Browser and processed to conduct gene density, clustering, and PCA analysis using
Python programming. 

## Objectives:
1. Clean and prepare genomic data
2. Compute gene density per chromosome
3. Perform PCA and Clustering
4. Generate summaries and visuals for downstream analysis

## Libraries Used:
* pandas
* numpy
* scikit-learn
* matplotlib 
* seaborn

## How to run program:

1. Place your annotation file (e.g., human_genome.txt) in the same directory as the script.
The file should be tab-delimited with transcript start and end coordinates.

2. Run the script:
   - python gene_density_analysis.py
3. The script will generate:

   - gene_density_summary.csv 
   - gene_density_pca_clusters.csv 
   - A PCA scatterplot visualization
