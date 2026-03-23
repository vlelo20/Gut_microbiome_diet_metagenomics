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

### **2.2.2 Sample Information**

### **2.2.3 Data Download**

### **2.2.4 Kraken2 Database**

## **2.3 Quality Control**

### **2.3.1 fastp**

### **2.3.2 MultiQC**

## **2.4 Taxonomic Classification**

### **2.4.1 Kraken2**

### **2.4.2 Bracken**

## **2.5 Diversity Analysis**

### **2.5.1 Alpha Diversity**

### **2.5.2 Beta Diversity**

### **2.5.3 Differential Abundance**

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
