# **Gut_microbiome_diet_metagenomics**

Author: Vian Lelo
Date created: March 20th, 2026, Last Updated: March 23rd, 2026

# Shotgun Metagenomics: Vegan vs Omnivore Gut Microbiome Analysis

---

# **1.0 Introduction**

Diet is among the most tractable and potent modulators of gut microbiome composition. Unlike host genetics, which account for only 2–8% of microbiome variance in twin studies (Goodrich et al., 2014), diet can remodel microbial community structure within 24–72 hours by exerting direct selective pressure on which taxa can metabolize available substrates — those that can outcompete those that cannot, shifting the ecological balance of the entire community (David et al., 2014). Plant-based diets, rich in fermentable dietary fibre and polyphenols, consistently enrich fibre-degrading taxa such as Faecalibacterium prausnitzii, Roseburia intestinalis, and Bifidobacterium spp. — all associated with anti-inflammatory SCFA production and gut barrier reinforcement (Sonnenburg & Bäckhed, 2016; Dahl et al., 2023). In contrast, animal-based diets high in protein and saturated fat promote bile-tolerant, potentially pro-inflammatory taxa including Bilophila wadsworthia and Bacteroides spp., alongside reductions in SCFA-producing Firmicutes (David et al., 2014; Zinöcker & Lindseth, 2018). Because vegans and omnivores represent a naturally occurring, sustained dietary contrast rather than a controlled short-term intervention, they offer a particularly valuable window into how habitual diet shapes the microbiome over time. This project therefore investigates whether dietary pattern significantly influences gut microbiome composition at the species level, using publicly available shotgun metagenomic data from Fragiadakis et al. (2020), a dataset capturing microbiome variation across individuals with well-characterized, habitual dietary patterns.

The choice of sequencing platform and analytical tools reflects deliberate consideration of their limitations and tradeoffs. Shotgun metagenomics was chosen over 16S rRNA amplicon sequencing because, while 16S metabarcoding is widely used and cost-effective, it targets only the hypervariable V3–V4 regions of the rRNA gene and is subject to amplification biases from primer mismatches, variable gene copy number across taxa, and PCR chimera formation — limiting resolution to the genus level in most cases (Schloss et al., 2011). This is a meaningful constraint when comparing dietary groups in which closely related species, such as Bacteroides thetaiotaomicron and Bacteroides fragilis, can have divergent metabolic roles and health associations. Shotgun metagenomics sequences all DNA in a sample in an untargeted manner, providing species- and strain-level resolution, access to functional gene content, and the ability to detect novel organisms without prior sequence knowledge — though at greater sequencing cost and computational demand, including the need for host read decontamination prior to classification (Quince et al., 2017; Durazzi et al., 2021).

Raw reads were quality-controlled using fastp v0.24.0 (Chen et al., 2018), selected over Trimmomatic for its integrated adapter auto-detection, polyG trimming, and quality filtering in a single pass, with reads below 50 bp discarded. Taxonomic classification was then performed with Kraken2 v2.1.6 (Wood et al., 2019) against the Standard-8 database (February 2026), which includes RefSeq bacterial, archaeal, viral, plasmid, human, and UniVec_Core sequences. While Kraken2 is orders of magnitude faster than alignment-based tools such as BLAST, k-mer methods are prone to false-positive classifications at default settings — particularly for taxa underrepresented in the reference database — which was mitigated by applying a confidence threshold of 0.15 above the default of 0 (Lu & Salzberg, 2020). Bracken v3.0 (Lu et al., 2017) was applied downstream to correct a systematic bias in raw Kraken2 output whereby longer genomes capture disproportionately more reads, inflating their apparent abundance; Bracken redistributes read counts probabilistically and normalizes for genome size, producing more accurate species-level estimates.

Diversity analysis was conducted in R using phyloseq (McMurdie & Holmes, 2013) and vegan (Oksanen et al., 2020). Three complementary alpha diversity metrics were calculated: observed species richness, which counts detected species without weighting; the Shannon index, which accounts for both richness and evenness and is sensitive to community-wide changes; and the Simpson index, which is dominance-weighted and more robust to rare species — together capturing aspects of within-sample diversity that no single metric reflects alone. Beta diversity was assessed using Bray-Curtis dissimilarity on relative abundances, visualized by PCoA, and tested with PERMANOVA (adonis2) to determine whether community composition differed significantly between dietary groups. Rarefaction curves were generated using vegan::rarecurve() to confirm adequate sequencing depth prior to analysis. Differential abundance was assessed using ANCOM-BC2 (Lin & Peddada, 2022), chosen over DESeq2 because microbiome sequencing yields inherently compositional data — a genuine increase in one taxon mathematically suppresses the relative abundance of others, violating the assumptions of DESeq2's negative binomial model, which was designed for RNA-seq (Gloor et al., 2017; Hawinkel et al., 2019). ANCOM-BC2 corrects for this via log-ratio transformations and explicitly accounts for unequal sampling fractions. Given the small sample size (n=3 per group), no taxa reached FDR < 0.05; the top 20 taxa ranked by p-value are therefore presented as exploratory findings. Interactive taxonomic visualizations were generated with KronaTools v2.8.1 (Ondov et al., 2011), providing a hierarchical, sample-level overview of community composition that complements the quantitative analyses.

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
All six samples passed quality control with `fastp/0.24.0`. Minimal read loss was observed after trimming across all samples, indicating high raw read quality. Read lengths were confirmed at 150 bp, consistent with the sequencing parameters reported in Fragiadakis et al. (2020). Rarefaction curves for all samples plateaued well before their maximum sequencing depth (Figure 5), confirming that sequencing depth was sufficient to capture the majority of detectable species diversity in each sample.

## **3.2 Taxonomic Abundance**
<img width="600" height="800" alt="01_taxonomic_abundance" src="https://github.com/user-attachments/assets/f4576c16-52ee-4ffb-966e-aac6662917c0" />

**Figure 1.** Relative abundance of the top 20 most abundant species across vegan and omnivore samples. Each bar represents one sample, with colours indicating species identity.

Taxonomic classification with Kraken2 and Bracken identified 223 species across all six samples. The top 20 most abundant taxa differed visibly between diet groups (Figure 1). Vegan samples showed higher relative abundance of *Blautia wexlerae*, *Bifidobacterium* spp., *Faecalibacterium prausnitzii*, and *Roseburia* spp., while omnivore samples were enriched in *Alistipes* spp., *Phocaeicola vulgatus*, and *Bilophila wadsworthia*. Considerable within-group variability was observed, particularly among vegan samples, with one sample (SRR8146944) showing a distinct composition compared to the other two vegans.

## **3.3 Alpha Diversity**
<img width="600" height="800" alt="02_alpha_diversity" src="https://github.com/user-attachments/assets/bcd171e6-7d4a-4885-b03e-f3198c4741a4" />

**Figure 2.** Alpha diversity measures (Observed species richness, Shannon index, Simpson index) compared between vegan and omnivore diet groups. Points represent individual samples; boxes show interquartile range.

Alpha diversity measures showed overlapping distributions between diet groups, with omnivores tending toward higher median Shannon and Simpson indices than vegans (Figure 2). However, high within-group variability — particularly in the vegan group — means these differences should be interpreted with caution given the small sample size (n=3 per group). Observed species richness was similarly variable, with vegans spanning a wider range than omnivores.

## **3.4 Beta Diversity**
<img width="600" height="800" alt="03_beta_diversity_PCoA" src="https://github.com/user-attachments/assets/70ed98a4-c547-460a-b575-b80c8e3feea5" />

**Figure 3.** Principal Coordinates Analysis (PCoA) of Bray-Curtis dissimilarity between all six samples, coloured by diet group. PC1 explains 52.1% and PC2 explains 25.4% of total variance.

PCoA of Bray-Curtis dissimilarity revealed partial separation between vegan and omnivore samples along PC1, which accounted for 52.1% of total variance (Figure 3). Two omnivore samples clustered tightly together in the upper right quadrant, while two vegan samples clustered in the upper left. One vegan sample (SRR8146944) was separated from the remaining vegan samples along PC2, suggesting notable inter-individual variability within the vegan group. One omnivore sample (SRR8146935) was also separated from the other omnivores, consistent with the high individual variability documented in human microbiome studies (Grice & Segre, 2012).

## **3.5 Differential Abundance**
<img width="600" height="800" alt="04_differential_abundance" src="https://github.com/user-attachments/assets/45c6e111-a4af-4dc1-bd09-cd0ef9dd7499" />

**Figure 4.** Top 20 taxa ranked by ANCOM-BC2 p-value, showing log fold change between vegan and omnivore groups. Green bars indicate taxa enriched in vegans; red bars indicate taxa enriched in omnivores. Results are exploratory due to small sample size (n=3 per group).

ANCOM-BC2 differential abundance analysis identified no taxa reaching FDR < 0.05, likely due to the limited statistical power of n=3 samples per group. Exploratory examination of the top 20 taxa ranked by p-value revealed that *Blautia wexlerae*, *Ruminococcus bicirculans*, *Anaerostipes hadrus*, and *Roseburia* spp. showed the largest positive log fold changes in vegans, while *Alistipes onderdonkii*, *Hominenteromicrobium mulieris*, and *Bilophila wadsworthia* showed the largest negative fold changes, indicating enrichment in omnivores (Figure 4).

## **3.6 Rarefraction Curves**
<img width="600" height="800" alt="05_rarefaction" src="https://github.com/user-attachments/assets/e721665b-dc83-431c-a9c5-d6ac994cb725" />

**Figure 5.** Rarefaction curves showing cumulative species discovery as a function of sequencing depth for all six samples. Green lines = vegan samples; red lines = omnivore samples.

Rarefaction analysis confirmed that sequencing depth was sufficient across all samples (Figure 5). All six curves reached a clear plateau well before their maximum read depth, indicating that additional sequencing would yield minimal new species discoveries. Omnivore samples (red) generally reached higher species counts at saturation (approximately 115–140 species) compared to vegan samples (green, approximately 70–120 species), consistent with the higher observed species richness noted in the alpha diversity analysis. The sample with the lowest species count (SRR8146944, vegan) plateaued at approximately 70 species, which may reflect genuine biological differences in community diversity rather than insufficient sequencing depth.

---

# **4.0 Discussion**
This study identified notable differences in gut microbiome composition between vegan and omnivore individuals, consistent with the established influence of long-term dietary patterns on microbial community structure (Sonnenburg & Bäckhed, 2016). The partial separation observed in the PCoA (Figure 3) suggests that diet is a meaningful driver of beta diversity, though individual-level variation remains substantial — a well-documented feature of the human gut microbiome (Grice & Segre, 2012).

The enrichment of *Blautia wexlerae*, *Roseburia* spp., and *Anaerostipes hadrus* in vegan samples is biologically consistent with a high-fibre, plant-based diet. These taxa are obligate anaerobes that ferment dietary fibre to produce short-chain fatty acids (SCFAs), particularly butyrate and acetate, which are critical for colonocyte energy metabolism and have well-established anti-inflammatory effects (Baxter et al., 2019). The higher abundance of *Bifidobacterium* spp. in vegans further supports this interpretation, as bifidobacterial growth is strongly promoted by plant-derived prebiotics such as inulin and fructooligosaccharides (Zmora et al., 2019).

In contrast, the enrichment of *Bilophila wadsworthia* in omnivores is consistent with the known association between animal-fat consumption and this sulfate-reducing bacterium. David et al. (2014) demonstrated that *Bilophila wadsworthia* blooms rapidly in response to animal-based diets and has been linked to intestinal inflammation through hydrogen sulfide production. Similarly, *Alistipes* spp. — enriched in omnivores here — have been associated with protein fermentation and are frequently elevated in individuals consuming high-protein, animal-based diets (Mosca et al., 2016).

The absence of statistically significant differential abundance results (FDR < 0.05) is most plausibly explained by insufficient statistical power at n=3 per group, rather than a true absence of biological signal. The directional consistency of fold changes with the published literature supports the biological relevance of the observed trends.

## ** Limitations**
The primary limitation of this analysis is the small sample size (n=3 per group), which severely restricts the statistical power of both beta diversity tests and differential abundance analysis. Future studies should include at least 10–20 samples per group to achieve adequate power with compositional methods such as ANCOM-BC2. Additionally, the Kraken2 Standard-8 database is a memory-reduced version of the full Standard database, which may reduce sensitivity for rare or poorly represented taxa. Individual confounders including age, geographic origin, antibiotic history, and body mass index were not accounted for in this analysis, which may contribute to the within-group variability observed across all diversity measures.

Overall, these results support the hypothesis that long‑term vegan and omnivore diets are associated with distinct gut microbiome compositions, particularly in SCFA‑producing taxa and bile‑tolerant species, although the small sample size limits formal statistical significance.

---

# **5.0 References**

Baxter, N. T., Schmidt, A. W., Venkataraman, A., Kim, K. S., Waldron, C., & Schmidt, T. M. (2019). Dynamics of human gut microbiota and short-chain fatty acids in response to dietary interventions with three fermentable fibers. mBio, 10(1), e02566-18.

Cryan, J. F., O’Riordan, K. J., Cowan, C. S. M., Sandhu, K. V., Bastiaanssen, T. F. S., Boehme, M., ... Dinan, T. G. (2019). The microbiota–gut–brain axis. Physiological Reviews, 99(4), 1877–2013.

Dahl, W. J., Rivero Mendoza, D., & Lambert, J. M. (2023). Diet, nutrients and the microbiome. Progress in Molecular Biology and Translational Science, 171, 237–263.

David, L. A., Maurice, C. F., Carmody, R. N., Gootenberg, D. B., Button, J. E., Wolfe, B. E., ... Turnbaugh, P. J. (2014). Diet rapidly and reproducibly alters the human gut microbiome. Nature, 505(7484), 559–563.

Durazzi, F., Sala, C., Castellani, G., Manfreda, G., Remondini, D., & Castellani, G. (2021). Comparison between 16S rRNA and shotgun sequencing data for the taxonomic characterization of the gut microbiota. Scientific Reports, 11, 3030.

Fragiadakis, G. K., Smits, S. A., Sonnenburg, E. D., Van Treuren, W., Reid, G., Knight, R., ... Sonnenburg, J. L. (2020). Links between environment, diet, and the hunter-gatherer microbiome. Cell Host & Microbe, 27(3), 380–391.

Gloor, G. B., Macklaim, J. M., Pawlowsky-Glahn, V., & Egozcue, J. J. (2017). Microbiome datasets are compositional: and this is not optional. Frontiers in Microbiology, 8, 2224.

Goodrich, J. K., Waters, J. L., Poole, A. C., Sutter, J. L., Koren, O., Blekhman, R., ... Ley, R. E. (2014). Human genetics shape the gut microbiome. Cell, 159(4), 789–799.

Grice, E. A., & Segre, J. A. (2012). The human microbiome: Our second genome. Annual Review of Genomics and Human Genetics, 13, 151–170.

Hawinkel, S., Mattiello, F., Bijnens, L., & Thas, O. (2019). A broken promise: Microbiome differential abundance methods do not control the false discovery rate. Briefings in Bioinformatics, 20(1), 210–221.

Lin, H., & Peddada, S. D. (2020). Analysis of compositions of microbiomes with bias correction. Nature Communications, 11, 3514.

Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15, 550.

Lu, J., Breitwieser, F. P., Thielen, P., & Salzberg, S. L. (2017). Bracken: Estimating species abundance in metagenomics data. PeerJ Computer Science, 3, e104.

Lu, J., & Salzberg, S. L. (2020). Ultrafast and accurate 16S rRNA microbial community analysis using Kraken 2. Genome Biology, 21, 170.

McMurdie, P. J., & Holmes, S. (2013). phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data. PLOS ONE, 8(4), e61217.

Mosca, A., Leclerc, M., & Hugot, J.-P. (2016). Gut microbiota diversity and human diseases: Should we reintroduce key predators in our ecosystem? Frontiers in Microbiology, 7, 455.

Oksanen, J., Blanchet, F. G., Friendly, M., Kindt, R., Legendre, P., McGlinn, D., ... Wagner, H. (2020). vegan: Community Ecology Package (Version 2.5-7) [R package]. https://CRAN.R-project.org/package=vegan

Quince, C., Walker, A. W., Simpson, J. T., Loman, N. J., & Segata, N. (2017). Shotgun metagenomics, from sampling to analysis. Nature Biotechnology, 35(9), 833–844.

Schloss, P. D., Westcott, S. L., Ryabin, T., Hall, J. R., Hartmann, M., Hollister, E. B., ... Weber, C. F. (2011). Assessing and improving methods used in operational taxonomic unit-based approaches for 16S rRNA gene sequence analysis. Applied and Environmental Microbiology, 77(10), 3219–3226.

Sonnenburg, J. L., & Bäckhed, F. (2016). Diet–microbiota interactions as moderators of human metabolism. Nature, 535(7610), 56–64.

Wood, D. E., Lu, J., & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. Genome Biology, 20, 257.

Zinöcker, M. K., & Lindseth, I. A. (2018). The Western diet–microbiome–host interaction and its role in metabolic disease. Nutrients, 10(3), 365.

Zmora, N., Zilberman-Schapira, G., Suez, J., Mor, U., Dori-Bachash, M., Bashiardes, S., ... Elinav, E. (2019). Personalized gut mucosal colonization resistance to empiric probiotics is associated with unique host and microbiome features. Cell, 174(6), 1388–1405.
