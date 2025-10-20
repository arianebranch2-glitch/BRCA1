Author: A'Riane Branch
Date: 15 Oct 2025

Overview:
I used the Ensembl REST API to get human BRCA1 (ENSG00000012048) variants. I computed GC-content for a representative allele and the genomic length for each variant, then plotted both.

Step 1 — Data retrieval
lookup/symbol/homo_sapiens/BRCA1 → Ensembl gene ID
overlap/id/ENSG00000012048?feature=variation → variant records (JSON)
Saved as brca1_variants.json
Total variants: 57,556

Step 2 — Processing and analysis
For each record:
Picked the first non-deletion allele from alleles (skipped “-”).
GC% = (G + C) / length × 100 for that allele.
Length = end − start + 1 (inclusive).
Wrote sequence_summary.csv with columns: variant_id,length,gc_percent.
Count with length = 100 bp: 1

Step 3 — Visualizations
Figure 1: Histogram of allele GC% (integer part).
Figure 2: Histogram of variant genomic lengths.
Saved as figure1.png and figure2.png.
Insert the two images here.

Results (short)
GC% shows discrete values (many very short alleles). See Figure 1.
Lengths are mostly small (SNPs/short indels) with a long tail. See Figure 2.

Files produced
final_analysis.py
brca1_variants.json
sequence_summary.csv
figure1.png, figure2.png
analysis_report.pdf (this document)
