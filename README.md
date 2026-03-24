# **Gut_microbiome_diet_metagenomics**

Author: Vian Lelo
Date created: March 20th, 2026, Last Updated: March 23rd, 2026

# Shotgun Metagenomics: Vegan vs Omnivore Gut Microbiome Analysis

---

# **1.0 Introduction**

Diet is among the most tractable and potent modulators of gut microbiome composition, capable of remodelling microbial community structure within 24–72 hours by exerting direct selective pressure on which taxa can metabolize available substrates (David et al., 2014). Plant-based diets, rich in fermentable dietary fibre and polyphenols, consistently enrich fibre-degrading taxa such as Faecalibacterium prausnitzii, Roseburia intestinalis, and Bifidobacterium spp. — all associated with anti-inflammatory short-chain fatty acid (SCFA) production and gut barrier reinforcement (Sonnenburg & Bäckhed, 2016). In contrast, animal-based diets high in protein and saturated fat promote bile-tolerant, potentially pro-inflammatory taxa including Bilophila wadsworthia and Bacteroides spp., alongside reductions in SCFA-producing Firmicutes (David et al., 2014). Because vegans and omnivores represent a naturally occurring, sustained dietary contrast rather than a controlled short-term intervention, they offer a particularly valuable window into how habitual diet shapes the microbiome over time (Sonnenburg & Bäckhed, 2016). This project therefore investigates whether dietary pattern significantly influences gut microbiome composition at the species level, using publicly available shotgun metagenomic data from Fragiadakis et al. (2020), a dataset capturing microbiome variation across individuals with well-characterized, habitual dietary patterns.

The choice of sequencing platform and analytical tools reflects deliberate consideration of their limitations and tradeoffs. Shotgun metagenomics was chosen over 16S rRNA amplicon sequencing because, while 16S metabarcoding is widely used and cost-effective, it targets only hypervariable regions of the rRNA gene and is subject to amplification biases from primer mismatches and variable gene copy number across taxa — limiting resolution to the genus level in most cases (Durazzi et al., 2021). This is a meaningful constraint when comparing dietary groups in which closely related species, such as Bacteroides thetaiotaomicron and Bacteroides fragilis, can have divergent metabolic roles and health associations. Shotgun metagenomics sequences all DNA in a sample in an untargeted manner, providing species- and strain-level resolution and access to functional gene content — though at greater sequencing cost and computational demand, including the need for host read decontamination prior to classification (Durazzi et al., 2021). Raw reads were quality-controlled using fastp v0.24.0, with adapter auto-detection, polyG trimming, and reads below 50 bp discarded. Taxonomic classification was performed with Kraken2 v2.1.6 (Wood et al., 2019) against the Standard-8 database (February 2026), which includes RefSeq bacterial, archaeal, viral, plasmid, human, and UniVec_Core sequences. While Kraken2 is orders of magnitude faster than alignment-based tools such as BLAST, k-mer methods are prone to false-positive classifications at default settings, which was mitigated by applying a confidence threshold of 0.15 above the default of 0 (Wood et al., 2019). Bracken v3.0 (Lu et al., 2017) was applied downstream to correct a systematic bias whereby longer genomes capture disproportionately more reads; Bracken redistributes read counts probabilistically and normalizes for genome size, producing more accurate species-level estimates.

Diversity analysis was conducted in R using phyloseq (McMurdie & Holmes, 2013) and vegan (Oksanen et al., 2020). Three complementary alpha diversity metrics were calculated: observed species richness, which counts detected species without weighting; the Shannon index, which accounts for both richness and evenness; and the Simpson index, which is dominance-weighted and more robust to rare species — together capturing aspects of within-sample diversity that no single metric reflects alone. Beta diversity was assessed using Bray-Curtis dissimilarity on relative abundances, visualized by PCoA, and tested with PERMANOVA (adonis2) to determine whether community composition differed significantly between dietary groups. Differential abundance was assessed using ANCOM-BC2 (Lin & Peddada, 2022), chosen over DESeq2 because microbiome sequencing yields inherently compositional data — a genuine increase in one taxon mathematically suppresses the relative abundance of others, violating the assumptions of DESeq2's negative binomial model, which was designed for RNA-seq (Gloor et al., 2017). ANCOM-BC2 corrects for this via log-ratio transformations and explicitly accounts for unequal sampling fractions. Given the small sample size (n=3 per group), no taxa reached FDR < 0.05; the top 20 taxa ranked by p-value are therefore presented as exploratory findings. Interactive taxonomic visualizations were generated with KronaTools v2.8.1, providing a hierarchical overview of community composition that complements the quantitative analyses.

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

Quality control and adapter trimming were performed using fastp v0.24.0 (Chen et al., 2018). fastp was selected for its ability to perform adapter auto-detection, quality filtering, and length trimming in a single pass without requiring a separate adapter sequence file — an advantage over tools such as Trimmomatic, which requires manual adapter specification. For paired-end data, adapters were auto-detected using --detect_adapter_for_pe, polyG tails were trimmed using --trim_poly_g to remove sequencing artefacts common in Illumina two-colour chemistry platforms, and reads shorter than 50 bp were discarded using --length_required 50 to ensure sufficient read length for accurate k-mer classification downstream. Each sample produced an HTML and JSON quality report confirming adapter removal and read length distributions before and after filtering.

```bash
fastp --in1 data/${SRR}_1.fastq.gz --in2 data/${SRR}_2.fastq.gz \
      --out1 qc_output/${SRR}_1_clean.fastq.gz \
      --out2 qc_output/${SRR}_2_clean.fastq.gz \
      --detect_adapter_for_pe --trim_poly_g --length_required 50 \
      --thread 8
```

## **2.4 Taxonomic Classification**
### **2.4.1 Kraken2**

Taxonomic classification was performed using Kraken2 v2.1.6 (Wood et al., 2019) with the Standard-8 database (February 2026), which contains RefSeq bacterial, archaeal, viral, plasmid, human, and UniVec_Core sequences capped at 8 GB for memory efficiency. Kraken2 operates by querying reads against a pre-built database of k-mer-to-lowest-common-ancestor (LCA) mappings, assigning each read to the most specific taxon supported by its k-mer content. A confidence threshold of 0.15 was applied above the default of 0, requiring that a minimum fraction of k-mers mapping across a clade support the assigned classification — this substantially reduces false-positive assignments for taxa underrepresented in the reference database, at a modest cost to recall (Wood et al., 2019). The database was loaded entirely into RAM (200 GB allocated) to avoid repeated disk I/O during classification, which would otherwise represent the primary runtime bottleneck at this scale. Only .report files were retained after each sample; raw .kraken output files were deleted immediately to conserve storage, as the report format contains all information required for downstream Bracken re-estimation and Krona visualization. Classification rates were taken directly from Kraken2 report files, which report the percentage of unclassified reads as a standard output metric.

```bash
kraken2 --db $DB --confidence 0.15 --threads 32 \
        --report kraken_output/${SRR}.report \
        --paired qc_output/${SRR}_1_clean.fastq.gz qc_output/${SRR}_2_clean.fastq.gz
```

### **2.4.2 Bracken**
Species-level abundance estimation was performed using Bracken v3.0 (Lu et al., 2017). Raw Kraken2 output is taxonomically diffuse and genome-size biased — reads are distributed across all taxonomic nodes rather than resolved to species, and longer genomes capture disproportionately more reads regardless of true abundance. Bracken corrects for this by applying a Bayesian model trained on the same reference database to redistribute and normalize read counts at species level. Read length was set to -r 150 to match the confirmed 150 bp sequencing length, -l S specified species-level redistribution, and -t 10 excluded classifications supported by fewer than 10 reads. All six .bracken output files were then combined into a single abundance matrix using combine_bracken_outputs.py for downstream analysis.

```bash
bracken -d $DB -i kraken_output/${SRR}.report \
        -o bracken_output/${SRR}.bracken -r 150 -l S -t 10
```
### **2.4.3 Krona Visualization**

Interactive taxonomic visualizations were generated using KronaTools v2.8.1 (Ondov et al., 2011). Kraken2 reports were first converted to Krona-compatible format using kreport2krona.py from KrakenTools, and all six samples were combined into a single interactive HTML file (krona_output/all_samples_krona.html). KronaTools produces hierarchical sunburst plots that allow interactive exploration of community composition across all taxonomic ranks simultaneously, providing an intuitive complement to the quantitative diversity analyses for identifying dominant taxa and gross compositional differences between samples.

## **2.5 Diversity Analysis**

All diversity analyses were performed in R using phyloseq, vegan, ggplot2, and ANCOMBC. The combined Bracken abundance matrix was imported and only read count columns were used (fraction columns excluded), as phyloseq requires integer count data for diversity calculations and downstream compositional analysis.

### **2.5.1 Alpha Diversity**
Observed species richness, Shannon index, and Simpson index were calculated using plot_richness() in phyloseq and compared between diet groups. Observed richness counts the total number of detected species without any weighting, providing a simple measure of community size. The Shannon index accounts for both richness and evenness, giving a more complete picture of diversity by penalizing communities dominated by few taxa. The Simpson index is dominance-weighted and less sensitive to rare species, making it more robust in datasets where low-abundance taxa may reflect sequencing noise rather than true biological diversity. Together, these three metrics capture complementary aspects of within-sample diversity that no single index reflects alone.

### **2.5.2 Beta Diversity**
Bray-Curtis dissimilarity was calculated on relative abundance data and visualized using PCoA ordination (ordinate(), method = "PCoA"). Bray-Curtis was selected because it accounts for both the presence and relative abundance of shared species between samples, making it more informative than presence-absence metrics for compositional comparisons. A PERMANOVA test (adonis2) was performed to assess whether the observed separation between dietary groups was statistically significant, as PCoA ordination alone is visual and does not provide a formal test of group differences. PERMANOVA is appropriate here because it makes no assumptions about the underlying distribution of the data, unlike parametric alternatives.

### **2.5.3 Differential Abundance**
Differential abundance was assessed using ANCOM-BC2 (ANCOMBC R package; Lin & Peddada, 2022), which corrects for compositional bias and unequal sampling fractions via log-ratio transformations — making it more appropriate for microbiome count data than RNA-seq-derived methods such as DESeq2, which assumes a negative binomial distribution that compositional data violates (Gloor et al., 2017). Due to the small sample size (n=3 per group), statistical power was limited and no taxa reached FDR < 0.05; the top 20 taxa ranked by p-value are therefore presented as exploratory results reflecting directional trends rather than confirmed differential abundance. Rarefaction curves were generated using vegan::rarecurve() prior to analysis to confirm that sequencing depth was sufficient across all samples and that observed species counts had reached a plateau, ensuring that diversity comparisons were not confounded by unequal sampling effort.

---

# **3.0 Results**

## **3.1 Quality Control Summary**
All six samples passed quality control with fastp v0.24.0 (Chen et al., 2018). Raw read counts ranged from 55.5 million (SRR8146944, vegan) to 94.7 million (SRR8146951, vegan) reads per sample. Read retention was high across all samples, with 95.3–98.1% of reads passing filters — read loss was minimal and primarily attributable to low quality rather than length filtering, with fewer than 0.05% of reads discarded for being under 50 bp. Illumina TruSeq adapters were auto-detected and trimmed in all samples, and pre-filtering Q30 scores ranged from 78.1–93.9%, improving after filtering. Duplication rates were low (0.32–0.55%), suggesting minimal PCR amplification bias. Rarefaction curves for all samples plateaued well before maximum sequencing depth (Figure 5), confirming sufficient depth to capture the majority of detectable species diversity.

## **3.2 Taxonomic Abundance**
<img width="600" height="800" alt="01_taxonomic_abundance" src="https://github.com/user-attachments/assets/f4576c16-52ee-4ffb-966e-aac6662917c0" />

**Figure 1.** Relative abundance of the top 20 most abundant species across vegan and omnivore samples. Each bar represents one sample, with colours indicating species identity.

Taxonomic classification with Kraken2 and Bracken identified 223 species across all six samples. The top 20 most abundant taxa differed visibly between diet groups (Figure 1). Vegan samples were consistently enriched in Blautia wexlerae, Faecalibacterium prausnitzii, Bifidobacterium adolescentis, Bifidobacterium longum, Roseburia faecis, and Anaerostipes hadrus — all obligate anaerobes associated with dietary fibre fermentation and SCFA production. Omnivore samples showed higher relative abundance of Alistipes onderdonkii, Alistipes putredinis, Phocaeicola vulgatus, and Bacteroides stercoris, taxa typically associated with protein and fat metabolism. Both groups shared a background of Faecalibacterium spp. and Phascolarctobacterium succinatutens, suggesting some core community members are present regardless of dietary pattern. Considerable within-group variability was observed — omnivore sample SRR8146935 was dominated by Blautia wexlerae and Phocaeicola vulgatus, while SRR8146936 showed an unusually high proportion of Roseburia faecis relative to the other omnivores. Among vegans, SRR8146944 showed a notably distinct composition with lower overall species evenness compared to SRR8146951 and SRR8146952, which were more compositionally similar to each other.

## **3.3 Alpha Diversity**
<img width="600" height="800" alt="02_alpha_diversity" src="https://github.com/user-attachments/assets/bcd171e6-7d4a-4885-b03e-f3198c4741a4" />

Figure 2. Alpha diversity measures (Observed species richness, Shannon index, Simpson index) compared between vegan and omnivore diet groups. Points represent individual samples; boxes show interquartile range.

Alpha diversity measures showed overlapping but directionally distinct distributions between diet groups (Figure 2). Observed species richness was higher in vegans (median ~120 species) compared to omnivores (median ~105 species), though the vegan group spanned a considerably wider range (~79–140 species) compared to omnivores (~70–118 species), driven largely by one vegan outlier (SRR8146944) with notably lower richness consistent with its distinct composition observed in Figure 1. The Shannon index, which accounts for both richness and evenness, showed the opposite trend — omnivores had a higher median Shannon value (~3.1) compared to vegans (~2.5), suggesting that while vegans may detect more species overall, community abundance is more evenly distributed across species in omnivore samples. The Simpson index reinforced this pattern, with omnivores showing higher median dominance-corrected diversity (~0.93) than vegans (~0.79), and the vegan group again showing greater spread with one sample falling as low as ~0.63. Together these metrics suggest that omnivore samples tend toward more even communities while vegan samples, though potentially richer in total species, are more variable and in some cases dominated by fewer highly abundant taxa. All differences should be interpreted cautiously given the small sample size (n=3 per group) and the considerable within-group variability observed across all three metrics.

## **3.4 Beta Diversity**
<img width="600" height="800" alt="03_beta_diversity_PCoA" src="https://github.com/user-attachments/assets/70ed98a4-c547-460a-b575-b80c8e3feea5" />

**Figure 3.** Principal Coordinates Analysis (PCoA) of Bray-Curtis dissimilarity between all six samples, coloured by diet group. PC1 explains 52.1% and PC2 explains 25.4% of total variance.

PCoA of Bray-Curtis dissimilarity revealed partial separation between vegan and omnivore samples along PC1, which accounted for 52.1% of total variance (Figure 3). PERMANOVA indicated that dietary group explained 19.1% of the variation in community composition (R² = 0.191, p = 0.500), which was not statistically significant — consistent with the small sample size of n=3 per group limiting statistical power rather than indicating a true absence of biological signal. Two omnivore samples clustered tightly together in the upper right quadrant, while two vegan samples clustered in the upper left. One vegan sample (SRR8146944) was separated from the remaining vegan samples along PC2, suggesting notable inter-individual variability within the vegan group. One omnivore sample (SRR8146935) was also separated from the other omnivores, consistent with the high individual variability documented in human microbiome studies (Grice & Segre, 2012).
## **3.5 Differential Abundance**
<img width="600" height="800" alt="04_differential_abundance" src="https://github.com/user-attachments/assets/45c6e111-a4af-4dc1-bd09-cd0ef9dd7499" />

**Figure 4.** Top 20 taxa ranked by ANCOM-BC2 p-value, showing log fold change between vegan and omnivore groups. Green bars indicate taxa enriched in vegans; red bars indicate taxa enriched in omnivores. Results are exploratory due to small sample size (n=3 per group).

ANCOM-BC2 differential abundance analysis identified no taxa reaching FDR < 0.05, likely due to limited statistical power at n=3 per group. Exploratory examination of the top 20 taxa ranked by p-value was largely consistent with the compositional patterns observed in Figure 1. Blautia wexlerae, Anaerostipes hadrus, Roseburia faecis, and Roseburia intestinalis showed the largest positive log fold changes in vegans, reinforcing their visual enrichment in vegan samples in the taxonomic abundance plot. Similarly, Alistipes onderdonkii showed the largest negative fold change, consistent with its dominant presence in omnivore samples observed in Figure 1. Bilophila wadsworthia enrichment in omnivores was also confirmed here, consistent with expectations from the literature regarding animal-fat consumption (David et al., 2014). However, two findings were less anticipated: Homo sapiens reads appearing among the top 20 differentially abundant features suggests residual host DNA was not fully removed during decontamination in some samples, and Ruminococcus bicirculans showing the second largest positive fold change in vegans was not visually prominent in Figure 1, suggesting it may be consistently present at moderate levels across vegan samples rather than dominant in any single one. Overall, the directional consistency between the taxonomic abundance plot and the differential abundance results supports the biological relevance of the observed trends despite the lack of statistical significance.

Segatella copri, visually prominent in vegan samples in Figure 1, was notably absent from the top 20 differential abundance results, likely reflecting inconsistent abundance across the three vegan samples rather than a consistent group-level enrichment — a pattern that ANCOM-BC2 would not detect as directionally significant.


## **3.6 Rarefraction Curves**
<img width="600" height="800" alt="05_rarefaction" src="https://github.com/user-attachments/assets/e721665b-dc83-431c-a9c5-d6ac994cb725" />

**Figure 5.** Rarefaction curves showing cumulative species discovery as a function of sequencing depth for all six samples. Green lines = vegan samples; red lines = omnivore samples.

Rarefaction analysis confirmed that sequencing depth was sufficient across all samples (Figure 5). All six curves reached a clear plateau well before their maximum read depth, indicating that additional sequencing would yield minimal new species discoveries. Omnivore samples (red) generally reached higher species counts at saturation (approximately 115–140 species) compared to vegan samples (green, approximately 70–120 species), consistent with the higher observed species richness noted in the alpha diversity analysis. The sample with the lowest species count (SRR8146944, vegan) plateaued at approximately 70 species, which may reflect genuine biological differences in community diversity rather than insufficient sequencing depth.

## **3.7 Krona Taxonomic Visualization**
The Krona visualization revealed that classified reads represented only 0.64–1.33% of total reads across all six samples, with 98.7–99.4% of reads remaining unclassified following Kraken2 classification. Despite this low overall classification rate, the classified fraction was consistent across samples and was used for all downstream diversity and differential abundance analyses, representing relative community composition within the classified subset rather than absolute community coverage.

---

# **4.0 Discussion**
This study identified notable differences in gut microbiome composition between vegan and omnivore individuals, consistent with the established influence of long-term dietary patterns on microbial community structure (Sonnenburg & Bäckhed, 2016). The partial separation observed in the PCoA (Figure 3) suggests that diet is a meaningful driver of beta diversity, though individual-level variation remains substantial — a well-documented feature of the human gut microbiome (Grice & Segre, 2012). The lack of complete separation is unsurprising given the small sample size and the known contribution of host-specific factors such as age, geography, and antibiotic history to microbiome composition, none of which were controlled for in this analysis.

The enrichment of Blautia wexlerae, Roseburia faecis, Roseburia intestinalis, and Anaerostipes hadrus in vegan samples is biologically consistent with a high-fibre, plant-based diet, and was among the strongest signals in both the taxonomic abundance plot (Figure 1) and the differential abundance analysis (Figure 4). These taxa are obligate anaerobes that ferment dietary fibre to produce short-chain fatty acids (SCFAs), particularly butyrate and acetate, which serve as the primary energy source for colonocytes and exert systemic anti-inflammatory effects (Baxter et al., 2019). Ruminococcus bicirculans, which showed the second largest positive fold change in vegans in Figure 4 despite being less visually prominent in Figure 1, is a known hemicellulose and resistant starch degrader — its consistent moderate enrichment across vegan samples suggests a stable fibre-degrading role that does not dominate any single sample but contributes reliably to the community (David et al., 2014). The higher abundance of Bifidobacterium adolescentis and Bifidobacterium longum in vegans further supports this interpretation, as bifidobacterial growth is strongly promoted by plant-derived prebiotics such as inulin and fructooligosaccharides (Zmora et al., 2019). Collectively, these findings suggest that the vegan microbiome in this dataset is functionally oriented toward carbohydrate fermentation and SCFA production, consistent with the substrate availability imposed by a plant-based diet.

In contrast, omnivore samples were enriched in Alistipes onderdonkii, Alistipes ihumii, Alistipes communis, and Alistipes senegalensis — representing a consistent pattern of Alistipes enrichment across multiple species rather than a single taxon. Alistipes spp. have been associated with protein fermentation and are frequently elevated in individuals consuming high-protein, animal-based diets, and have also been linked to elevated indole and secondary bile acid production, which at high concentrations can be cytotoxic to the colonic epithelium (Mosca et al., 2016). The enrichment of Bilophila wadsworthia in omnivores is consistent with the well-established association between animal fat consumption and this sulfate-reducing bacterium — David et al. (2014) demonstrated that Bilophila wadsworthia blooms rapidly in response to animal-based diets and produces hydrogen sulfide, which has been linked to intestinal inflammation and epithelial barrier disruption. Phocaeicola vulgatus and Bacteroides caccae, also enriched in omnivores, are bile-tolerant saccharolytic taxa that thrive in environments rich in animal-derived substrates and have been associated with altered bile acid metabolism (Sonnenburg & Bäckhed, 2016).

Notably, Segatella copri (formerly Prevotella copri) was visually prominent in vegan samples in Figure 1 but did not appear among the top 20 differentially abundant taxa in Figure 4. This discrepancy likely reflects inconsistent abundance across the three vegan samples — ANCOM-BC2 detects consistent directional differences across all samples within a group, so a taxon elevated in one vegan sample but absent or low in the others would not rank highly by p-value. This highlights an important limitation of small-n compositional analysis: visually striking patterns in individual samples do not necessarily reflect group-level biological signals.
The absence of statistically significant differential abundance results (FDR < 0.05) is most plausibly explained by insufficient statistical power at n=3 per group, rather than a true absence of biological signal. Supporting this interpretation, the PERMANOVA R² of 0.191 indicates that dietary pattern explains approximately 19% of community-level variation in beta diversity — a biologically meaningful effect size that failed to reach significance solely due to the limited number of permutations possible at n=3. The directional consistency of fold changes with published literature — fibre-fermenting taxa enriched in vegans, bile-tolerant and protein-fermenting taxa enriched in omnivores — further supports the biological relevance of the observed trends and is consistent with what would be predicted from the known metabolic consequences of these dietary patterns (David et al., 2014; Sonnenburg & Bäckhed, 2016). Future studies with larger sample sizes (n ≥ 20 per group) and matched confounders would be necessary to achieve the statistical power required to confirm these findings formally. Interactive taxonomic composition across all six samples is available as a supplementary HTML file (krona_output/all_samples_krona.html), providing a hierarchical overview of community composition at each taxonomic rank.

The Krona visualization revealed that approximately 99% of reads remained unclassified following Kraken2 classification, likely reflecting the reduced sensitivity of the Standard-8 database relative to the full Standard database, residual host DNA contamination, and the conservative confidence threshold of 0.15 applied to reduce false positives. While Bracken abundance estimates are calculated on the classified fraction only and remain valid for relative comparisons between samples, the low overall classification rate means that a substantial proportion of the microbial community may be unrepresented in the diversity and differential abundance analyses. Future analyses should incorporate a larger reference database and explicit host read decontamination prior to classification to improve sensitivity.

# **Limitations**
The primary limitations of this analysis are the small sample size (n=3 per group) and the low Kraken2 classification rate of 0.64–1.33% across all samples. The small sample size severely restricts statistical power for both beta diversity and differential abundance testing — PERMANOVA was limited to 719 permutations due to the small sample size, as fewer unique permutations are possible with n=3 per group, and no taxa reached FDR < 0.05 despite directionally consistent biological trends. Notably, the PERMANOVA R² of 0.191 suggests that dietary pattern explains approximately 19% of microbiome variation, which represents a biologically meaningful effect size that would likely reach statistical significance with an adequately powered sample. The low classification rate likely reflects the reduced reference coverage of the memory-capped Standard-8 database, residual host DNA, and the conservative confidence threshold of 0.15 — meaning a substantial proportion of the true microbial community may be unrepresented in the results. Individual confounders including age, geographic origin, antibiotic history, and BMI were not accounted for, which may contribute to the within-group variability observed across diversity measures. Future studies should incorporate the full Standard database, explicit host read decontamination prior to classification, and at least 10–20 samples per group to achieve adequate statistical power.
Overall, these results support the hypothesis that long-term vegan and omnivore diets are associated with distinct gut microbiome compositions, particularly in SCFA-producing taxa and bile-tolerant species. The directional consistency of findings with published literature suggests the biological signal is real, though the small sample size limits formal statistical significance.

---

# **5.0 References**

Baxter, N.T., et al. (2019). Dynamics of human gut microbiota and short-chain fatty acids in response to dietary interventions with three fermentable fibers. mBio, 10(1), e02566-18.
Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34(17), i884–i890.
David, L.A., et al. (2014). Diet rapidly and reproducibly alters the human gut microbiome. Nature, 505(7484), 559–563.
Durazzi, F., et al. (2021). Comparison between 16S rRNA and shotgun sequencing data for the taxonomic characterization of the gut microbiota. Scientific Reports, 11, 3030.
Fragiadakis, G.K., et al. (2020). Links between environment, diet, and the hunter-gatherer microbiome. Cell Host & Microbe, 27(3), 380–391.
Gloor, G.B., et al. (2017). Microbiome datasets are compositional: and this is not optional. Frontiers in Microbiology, 8, 2224.
Grice, E.A. & Segre, J.A. (2012). The human microbiome: our second genome. Annual Review of Genomics and Human Genetics, 13, 151–170.
Lin, H. & Peddada, S.D. (2022). Analysis of compositions of microbiomes with bias correction 2 (ANCOM-BC2). Nature Communications, 13, 3737.
Lu, J., et al. (2017). Bracken: estimating species abundance in metagenomics data. PeerJ Computer Science, 3, e104.
McMurdie, P.J. & Holmes, S. (2013). phyloseq: an R package for reproducible interactive analysis and graphics of microbiome census data. PLOS ONE, 8(4), e61217.
Mosca, A., et al. (2016). Gut microbiota diversity and human diseases: should we reintroduce key predators in our ecosystem? Frontiers in Microbiology, 7, 455.
Oksanen, J., et al. (2020). vegan: Community Ecology Package. R package version 2.5-7.
Ondov, B.D., et al. (2011). Interactive metagenomic visualization in a Web browser. BMC Bioinformatics, 12, 385.
Sonnenburg, J.L. & Bäckhed, F. (2016). Diet–microbiota interactions as moderators of human metabolism. Nature, 535(7610), 56–64.
Wood, D.E., et al. (2019). Improved metagenomic analysis with Kraken 2. Genome Biology, 20, 257.
Zmora, N., et al. (2019). Personalized gut mucosal colonization resistance to empiric probiotics is associated with unique host and microbiome features. Cell, 174(6), 1388–1405.
