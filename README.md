## Tumor Agnostic UMI

A tumor-agnostic approach for detection of colorectal neoplasm in ultra-deep targeted sequencing of a narrow subset of the genome.

## Quick start

1. **Requirements:** R with the packages listed in
   [DEPENDENCIES.md](DEPENDENCIES.md) (Bioconductor's `deepSNV`/
   `shearwater` for variant calling, plus ~29 CRAN/Bioconductor packages
   used across the analysis scripts), and Python with `pysam`/`numpy`/
   `gwf` for the fragment-length extraction step.
2. **This is analysis code developed on Aarhus University's GenomeDK,**
   a high-performance computing (HPC) cluster.
3. **`workflow.py`** shows the intended data-processing order: it runs
   `UMI_seq_fragment_length.py` (fragment-length extraction from
   consensus BAM files against `data/NEW_METHOD_hg38_08feb2016_capture_targets.bed`)
   as a `gwf` cluster workflow. The various `.R` scripts (`fileInput.R`,
   `readWes*.R`, `blacklisting.R`, `mahalanobis*.R`, `errorRates.R`, etc.)
   perform the downstream statistical analysis and scoring described
   below.

## Repository contents

- **`data/`** -- `IMPROVEptList`, the WES mutation CSV, the capture-panel
  BED file.
- **`figures/`** -- supplementary image(s).
- **R analysis scripts** -- each is a standalone script, meant to be run
  interactively / edited per-analysis rather than executed end-to-end as
  a pipeline.
- **`workflow.py` / `UMI_seq_fragment_length.py`** -- the Python side,
  for extracting fragment-length data from BAM files via a `gwf` cluster
  workflow.
- **License:** see [LICENSE](LICENSE) -- research/educational use.

## About

I developed a tumor-agnostic approach for the detection of colorectal neoplasm based on the DNA fragment length distribution on- and off-panel. Tumors exhibit clonal architecture with multiple mutations present in a subset of the tumor cell population. Sub-clonal mutations provide key insights into tumor evolution and sub-clonal variants can be detected using the deep coverage (over 10,000x) of next-generation sequencing (NGS) data, but their distinction from sequencing errors, library preparation and alignment artifacts depends on the noise level.

We perform variant calling via a package available at Bioconductor, which provides quantitative detection of sub-clonal mutations in ultra-deep sequencing data. The deepSNV algorithm is used for a comparative setup with a control experiment of the same loci and uses a beta-binomial model as well as a likelihood ratio test to discriminate sequencing errors and sub-clonal SNVs. The shearwater (SW) algorithm computes a Bayes classifier based on a beta-binomial model for variant calling with multiple samples for precisely estimating the detection error rates and the variance dispersion using prior knowledge, such as variation data from the COSMIC database. It computes a model based on the coverage, error rate, and a dispersion factor. The determination of statistical significance is made by comparing the expected number of nucleotide counts per read with the stringency level of the selected cut-off value based on computing a posterior probability as a function of a Bayesian factor.

**Analysis strategy:** We perform both tumor-informed and non-tumor informed analysis. We use a prior in the SW call to improve specificity by making it harder to call mutations not seen in COSMIC. In the tumor-informed analysis of the Dilution Series (DS) of a mixed tumor DNA set, there is, in principle, no need for using COSMIC as a prior, as we know that the mutations are in the tumor, and hence in the plasma too. We assess how different posterior cut off values impact specificity by analyzing the Cancer Research UK (CRUK) and Qiagen healthy samples using the DS mutation profile of 90 mutations with the objective to identify a tradeoff between sensitivity and specificity, and to make sure none of the control samples is called positive. We also assess the specificity by switching the mutation profiles of patients and assess the rate with which we call plasma samples positive when using the mutation profiles of other patient samples from pre-operative (preOP) clinical trial samples.

In non-tumor informed analysis, our requirements for specificity are much higher. Here, we search all 15,465 positions (the number will be lowered to around 15,000 based on the exact regions sequenced) in the panel and therefore it becomes critical not to make false positive calls. Using COSMIC as a prior is one strategy we use to improve specificity. Another option we have is to require a mutation to be observed on both strands (the AND algorithm).

**Limit of detection:** In every sample, there are a certain number of cfDNA molecules, which, after library preparation, are sequenced with a certain depth. In plasma, the tumors shed both mutated and non-mutated fragments, which is further mixed with normal cell fragments. We analyze the distributions of the error counts per transition to identify those resulting from clonal hematopoiesis vs the true noise on the panel. We further utilize background noise models based on using patient samples, patient-control samples, samples from healthy individuals, and dilution series to perform error suppression, which improves the accuracy of our analysis. For samples for which we have whole exome sequencing (WES) tumor tissue sequencing to guide the analysis, we require at least one patient mutation to be detected in the plasma to call it positive. For every patient tumor, we ask the question: What is the lowest ratio of tumor fragments over the total molecules for which we detect >2 tumor fragments?

I investigated how the patterns of length characteristics of the tumor fragments differ from those originating from normal cells for particular fragments lengths and genomic regions. I am interested in computing a comprehensive score and in classifying samples according to the likelihood of detecting cancer.

**Tumor informed strategy:** To quantitatively evaluate the detection limits of our assay, we have generated samples consisting of DNA from 36 tumors with known mutations, which have been pooled together. For these 36 tumors, we have an a-priori knowledge that 12 samples carry a KRAS mutation, 12 samples carry a BRAF mutation, and 12 samples carry a PIK3CA mutation, among other 90 mutations representative of colorectal cancer (CRC). In addition, we have gradually lowered the tumor concentration by mixing the tumor with normal DNA and keeping a constant amount of DNA corresponding to 25,000 GE. This way, we have created new titrated samples and have gradually reduced the frequency of the mutations with the objective to benchmark the performance of our algorithm. We have evaluated the mutation detection error rates for (i) different read coverage levels, (ii) mutation detection stringency cut-off levels, and (iii) ranges of detectability in terms of variant allele frequency (VAF). We combined the 8 samples and considered six different ranges ([0.1-0.5] with 2 values, [0.01-0.1) with 49 values, [0.001-0.01) with 163 values, [0.0001-0.001) with 255 values, (0-0.0001) with 105 values, and 146 entries have VAF=0) in terms of variant allele frequency (VAFs) and evaluated the success rate of detecting the 720 (90x8) mutations from our panel; see the table below. To evaluate the predictive power of our analysis, we compute specificity as True Negative / (True Negative + False Positive) and sensitivity as True Positive / (True Positive + False Negative). In the long run, it will be our objective to perform liquid biopsies only. Specifically, the current tumor-informed procedure includes a WES step of resected tissue to identify the mutations for each patient, which could be omitted and instead the entire unique molecular identifiers sequencing (UMI-seq) panel could be analyzed for significant mutations.

**Aggregating evidence from multiple mutations (non-tumor informed mutation score):** For classification purposes, I proposed to compute a score for each patient sample based on aggregating evidence from multiple mutations present in the plasma. I assign different weights for each mutation based on a metric (Matov et al. 2011), such as the Mahalanobis distance, which accounts for the heterogeneity of each mutation in samples from healthy individuals. In addition, I assign another weighting factor based on the likelihood to have a mutation at the concrete position with the gene and the probability the gene to be mutated based on CRC statistics and added a weight for short/long vs middle-length reads to further distinguish reads originating from tumor cells. I formulate the score as a function of multiple mutations, between five and eight, which collectively determine the presence of colorectal neoplasm. Each mutation is presented as its MAF weighted by the dispersion of the read counts for the same position in a cohort of normal samples. For positions with zero variance across the panel of normal (PON) samples, the zero is replaced with a very small value calculated as the minimal value of the variance matrix divided by 1e+07. Additionally, I utilize a second weight term based on the probability of each of the mutated genes to carry a mutation in that position. For the analysis of patient samples, I compare the Poisson distributions of the score values for a cohort of healthy individuals (normal) and CRC patients to determine a cut-off threshold score value. Further, I add weights based on the DNA fragment length distribution for particular on-panel and off-panel genomic positions.

For detailed information, see: https://www.researchgate.net/publication/382900905_Colorectal_Cancer_DNA_Fragmentation_Patterns_as_Image_Datasets


