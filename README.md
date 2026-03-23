# **Gut_microbiome_diet_metagenomics**

Author: Vian Lelo
Date created: March 20th, 2026, Last Updated: March 23rd, 2026

# Shotgun Metagenomics: Vegan vs Omnivore Gut Microbiome Analysis

---

# **1.0 Introduction**

The human gut microbiome is a complex community of microorganisms that plays a critical role in host metabolism, immune function, and overall health. Diet is one of the most influential modulators of gut microbiome composition, with plant-based and animal-based diets shown to differentially shape microbial diversity and species abundance (Sonnenburg & Bäckhed, 2016; David et al., 2014). This project investigates whether dietary pattern (vegan vs. omnivore) significantly influences gut microbiome composition at the species level, using shotgun metagenomics data from Fragiadakis et al. (2020). Unlike 16S metabarcoding, shotgun metagenomics sequences all DNA in a sample, providing higher taxonomic resolution and avoiding amplification bias associated with primer selection (Durazzi et al., 2021).

Taxonomic classification was performed using **Kraken2** (Wood et al., 2019), a k-mer-based approach that is substantially faster than alignment-based tools such as BLAST or Kaiju while maintaining high sensitivity at an appropriate confidence threshold. A confidence threshold of 0.15 was applied above the default of 0 to reduce false positive classifications — a known limitation of k-mer methods at default settings. **Bracken** (Lu et al., 2017) was used downstream to redistribute genus-level reads to species level and normalize for genome size, providing more accurate species-level abundance estimates than Kraken2 alone. Quality control was performed with **fastp**, which performs adapter trimming and quality filtering in a single efficient pass. Diversity analysis used **phyloseq** and **vegan**, the standard R packages for microbiome community analysis (McMurdie & Holmes, 2013; Oksanen et al., 2020). **ANCOM-BC2** (Lin & Peddada, 2020) was selected for differential abundance over DESeq2 because it is specifically designed for compositional microbiome data and corrects for unequal sampling fractions, making it more appropriate than RNA-seq-derived methods applied to metagenomics. Interactive taxonomic visualizations were generated using **KronaTools** to provide an intuitive overview of community composition across samples.

---

# **2.0 Methods**

## **2.1 Environment Setup**
All analyses were performed on the **Narval HPC cluster** (Digital Alliance of Canada) using SLURM for job scheduling. Downloads were performed directly on the login node inside `tmux` sessions, as compute nodes do not have external internet access.

```bash
module load sra-toolkit/3.0.9
module load fastp/0.24.0
module load kraken2/2.1.6
module load bracken/3.0
module load kronatools/2.8.1
```

### **2.1.1 Directory Structure**

## **2.2 Data & Metadata**

### **2.2.1 Data Source**
Raw shotgun metagenomics sequences were obtained from:

> Fragiadakis, G.K. et al. (2020). Links between environment, diet, and the hunter-gatherer microbiome. *Cell Host & Microbe*, 27(3), 380–391.
> NCBI SRA Accession: [SRP126540](https://www.ncbi.nlm.nih.gov/sra/?term=SRP126540)

### **2.2.2 Sample Information**

| SRR Accession | Diet Group | Read Length |
|---|---|---|
| SRR8146951 | Vegan | 150 bp |
| SRR8146952 | Vegan | 150 bp |
| SRR8146944 | Vegan | 150 bp |
| SRR8146936 | Omnivore | 150 bp |
| SRR8146935 | Omnivore | 150 bp |
| SRR8146938 | Omnivore | 150 bp |

### **2.2.3 Data Download**
Samples were downloaded on the Narval login node using `fasterq-dump` inside a `tmux` session:

```bash
tmux new -s downloads
module load sra-toolkit/3.0.9

for SRR in SRR8146951 SRR8146952 SRR8146944 SRR8146936 SRR8146935 SRR8146938; do
    fasterq-dump $SRR --split-files --threads 6 --skip-technical
    gzip ${SRR}*.fastq
done
```
### **2.2.4 Kraken2 Database**
The Kraken2 **Standard-8 database** (February 2026) was downloaded and extracted on the login node. This database contains RefSeq bacteria, archaea, viral, plasmid, human, and UniVec_Core sequences capped at 8 GB.

```bash
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08_GB_20260226.tar.gz
tar -xzf k2_standard_08_GB_20260226.tar.gz -C kraken_db/
```

## **2.3 Quality Control**


### **2.3.1 fastp**

Quality control and adapter trimming were performed using `fastp/0.24.0`. Adapters were auto-detected for paired-end data, polyG tails were trimmed, and reads shorter than 50 bp were discarded. Each sample produced an HTML and JSON report.

```bash
fastp --in1 data/${SRR}_1.fastq.gz --in2 data/${SRR}_2.fastq.gz \
      --out1 qc_output/${SRR}_1_clean.fastq.gz \
      --out2 qc_output/${SRR}_2_clean.fastq.gz \
      --detect_adapter_for_pe --trim_poly_g --length_required 50 \
      --thread 8
```

## **2.4 Taxonomic Classification**
### **2.4.1 Kraken2**

Taxonomic classification was performed using `kraken2/2.1.6` with the Standard-8 database. A confidence threshold of 0.15 was applied to reduce false positives. The database was loaded into RAM (200 GB allocated) for faster processing. Raw `.kraken` output files were deleted after each sample to conserve storage.

```bash
kraken2 --db $DB --confidence 0.15 --threads 32 \
        --report kraken_output/${SRR}.report \
        --paired qc_output/${SRR}_1_clean.fastq.gz qc_output/${SRR}_2_clean.fastq.gz
```

### **2.4.2 Bracken**
Bracken (`bracken/3.0`) re-estimated species-level abundances from Kraken2 reports by redistributing genus-level reads and normalizing for genome size. Read length was set to 150 bp to match confirmed sequencing length. All six `.bracken` output files were combined into a single abundance matrix using `combine_bracken_outputs.py`.

```bash
bracken -d $DB -i kraken_output/${SRR}.report \
        -o bracken_output/${SRR}.bracken -r 150 -l S -t 10
```
### **2.4.3 Krona Visualization**

Interactive taxonomic sunburst plots were generated using `kronatools/2.8.1` and `kreport2krona.py` (KrakenTools). All six samples were combined into a single interactive HTML file (`krona_output/all_samples_krona.html`).

## **2.5 Diversity Analysis**

All diversity analyses were performed in R using `phyloseq`, `vegan`, `ggplot2`, and `ANCOMBC`. The combined Bracken abundance matrix was imported and only read count columns were used (fraction columns excluded).

### **2.5.1 Alpha Diversity**
Observed species richness, Shannon index, and Simpson index were calculated using `plot_richness()` in phyloseq and compared between diet groups.

### **2.5.2 Beta Diversity**
Bray-Curtis dissimilarity was calculated on relative abundance data and visualized using PCoA ordination (`ordinate()`, method = "PCoA"). A PERMANOVA test (`adonis2`) was performed to assess statistical significance of group separation.

### **2.5.3 Differential Abundance**
Differential abundance was assessed using **ANCOM-BC2** (`ANCOMBC` R package), which corrects for compositional bias and unequal sampling fractions. Due to the small sample size (n=3 per group), no taxa reached FDR < 0.05. The top 20 taxa ranked by p-value are presented as exploratory results. Rarefaction curves were generated using `vegan::rarecurve()` to confirm sufficient sequencing depth.

---

# **3.0 Results**

## **3.1 Quality Control Summary**

## **3.2 Taxonomic Abundance**

## **3.3 Alpha Diversity**

## **3.4 Beta Diversity**

## **3.5 Differential Abundance**

---

# **4.0 Discussion**

## **4.1 Biological Interpretation**

## **4.2 Key Taxa**

## **4.3 Limitations**
# 5.0 References

---

# **6.0 References**
