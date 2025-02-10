# RNA-seq paper intermediate processing

<!-- ## Dependencies

* Singularity (v. 3)
* NextFlow (v. >= 21.04.0)
* python (v. 3) -->

<!-- 
## Setup
 
1. Update the `ROOT` variable at the top of the Makefile to the current working directory.
2. Update the `--bind` paths in the `SIF` variable on the second line of the Makefile; bind any paths that contain input data.
3. Update `nextflow.config` as necessary for your compute environment.
4. Get singularity container: `make singularity` -->

## Data

The `data/` directory has the following structure:

```bash
data
├── ancestry
│   └── ancestries.txt
├── crossmap-gencode-v30
│   ├── snp_mappability_100mer_2mismatch.bed.gz -> /net/topmed10/working/porchard/rnaseq/work/crossmap-gencode-v30/results/snp-mappability/snp_mappability_100mer_2mismatch.bed.gz
│   └── snp_mappability_100mer_2mismatch.bed.gz.tbi -> /net/topmed10/working/porchard/rnaseq/work/crossmap-gencode-v30/results/snp-mappability/snp_mappability_100mer_2mismatch.bed.gz.tbi
├── ensembl-regulatory-build
│   └── homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20221007.gff.gz # wget https://ftp.ensembl.org/pub/current_regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20221007.gff.gz
├── fasta
│   └── hg38
│       └── Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta
├── gene-counts
│   ├── all-cohorts.gene_counts.hdf5
│   └── all-cohorts.gene_tpm.hdf5
├── genotype-pca
│   ├── PC-scores.txt
│   └── PC-variance-explained.txt
├── genotypes
│   ├── plink-pass-filter -> /net/topmed10/working/porchard/rnaseq/work/subset-topmed-bcf/freeze-alpha/results/plink-pass-filter
│   └── vcfs-updated-ids-pass-filter -> /net/topmed10/working/porchard/rnaseq/work/subset-topmed-bcf/freeze-alpha/results/vcfs-updated-ids-pass-filter
├── gtf
│   └── gencode.v30.GRCh38.ERCC.genes.collapsed_only.gtf.gz
├── metadata
│   └── metadata.tm.txt
├── roadmap
│   ├── E029_15_coreMarks_hg38lift_mnemonics.bed.gz # wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E029_15_coreMarks_hg38lift_mnemonics.bed.gz
│   ├── E034_15_coreMarks_hg38lift_mnemonics.bed.gz
│   ├── E062_15_coreMarks_hg38lift_mnemonics.bed.gz
│   ├── E096_15_coreMarks_hg38lift_mnemonics.bed.gz
│   └── E114_15_coreMarks_hg38lift_mnemonics.bed.gz
├── scan-results -> /net/topmed11/working/porchard/rnaseq-2024-10-16-sample-update/qtl-scan-output
├── scan-samples # lists of samples to use for QTL scans (one library ID per line, no header)
│   ├── samples-to-use-for-scan.Lung.tsv
│   ├── samples-to-use-for-scan.Monocyte.tsv
│   ├── samples-to-use-for-scan.Nasal_epithelial.tsv
│   ├── samples-to-use-for-scan.PBMC.tsv
│   ├── samples-to-use-for-scan.T_cell.tsv
│   └── samples-to-use-for-scan.Whole_blood.tsv
├── sites # topmed freeze 9b BCF files
│   ├── freeze.9b.chr10.pass_and_fail.sites.bcf -> /net/topmed2/working/gt-release/exchange-area/freeze.9b/sites/freeze.9b.chr10.pass_and_fail.sites.bcf
│   ├── freeze.9b.chr10.pass_and_fail.sites.bcf.csi -> /net/topmed2/working/gt-release/exchange-area/freeze.9b/sites/freeze.9b.chr10.pass_and_fail.sites.bcf.csi
│   ├── freeze.9b.chr11.pass_and_fail.sites.bcf -> /net/topmed2/working/gt-release/exchange-area/freeze.9b/sites/freeze.9b.chr11.pass_and_fail.sites.bcf
│   ├── freeze.9b.chr11.pass_and_fail.sites.bcf.csi -> /net/topmed2/working/gt-release/exchange-area/freeze.9b/sites/freeze.9b.chr11.pass_and_fail.sites.bcf.csi
│   ├── ...
│   ├── freeze.9b.chr9.pass_and_fail.sites.bcf -> /net/topmed2/working/gt-release/exchange-area/freeze.9b/sites/freeze.9b.chr9.pass_and_fail.sites.bcf
│   ├── freeze.9b.chr9.pass_and_fail.sites.bcf.csi -> /net/topmed2/working/gt-release/exchange-area/freeze.9b/sites/freeze.9b.chr9.pass_and_fail.sites.bcf.csi
│   ├── freeze.9b.chrX.pass_and_fail.sites.bcf -> /net/topmed2/working/gt-release/exchange-area/freeze.9b/sites/freeze.9b.chrX.pass_and_fail.sites.bcf
│   └── freeze.9b.chrX.pass_and_fail.sites.bcf.csi -> /net/topmed2/working/gt-release/exchange-area/freeze.9b/sites/freeze.9b.chrX.pass_and_fail.sites.bcf.csi
├── tensorqtl-in # Gene expression/splicing phenotype files and covariates formatted for tensorQTL
│   ├── ancestry-specific
│   │   ├── cis-eqtl
│   │   │   ├── Lung___EUR.tensorqtl-in.20.covariates.tsv
│   │   │   ├── Lung___EUR.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── Monocyte___EUR.tensorqtl-in.15.covariates.tsv
│   │   │   ├── Monocyte___EUR.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── Nasal_epithelial___EUR.tensorqtl-in.20.covariates.tsv
│   │   │   ├── Nasal_epithelial___EUR.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── PBMC___AFR.tensorqtl-in.20.covariates.tsv
│   │   │   ├── PBMC___AFR.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── PBMC___EAS.tensorqtl-in.10.covariates.tsv
│   │   │   ├── PBMC___EAS.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── PBMC___EUR.tensorqtl-in.20.covariates.tsv
│   │   │   ├── PBMC___EUR.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── T_cell___EUR.tensorqtl-in.15.covariates.tsv
│   │   │   ├── T_cell___EUR.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── Whole_blood___AFR.tensorqtl-in.40.covariates.tsv
│   │   │   ├── Whole_blood___AFR.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── Whole_blood___EUR.tensorqtl-in.40.covariates.tsv
│   │   │   └── Whole_blood___EUR.tensorqtl-in.phenotypes.bed.gz
│   │   ├── cis-sqtl
│   │   │   ├── Lung___EUR.tensorqtl-in.10.covariates.tsv
│   │   │   ├── Lung___EUR.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── Monocyte___EUR.tensorqtl-in.5.covariates.tsv
│   │   │   ├── Monocyte___EUR.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── Nasal_epithelial___EUR.tensorqtl-in.10.covariates.tsv
│   │   │   ├── Nasal_epithelial___EUR.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── PBMC___AFR.tensorqtl-in.5.covariates.tsv
│   │   │   ├── PBMC___AFR.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── PBMC___EAS.tensorqtl-in.5.covariates.tsv
│   │   │   ├── PBMC___EAS.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── PBMC___EUR.tensorqtl-in.10.covariates.tsv
│   │   │   ├── PBMC___EUR.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── T_cell___EUR.tensorqtl-in.5.covariates.tsv
│   │   │   ├── T_cell___EUR.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── Whole_blood___AFR.tensorqtl-in.10.covariates.tsv
│   │   │   ├── Whole_blood___AFR.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── Whole_blood___EUR.tensorqtl-in.10.covariates.tsv
│   │   │   └── Whole_blood___EUR.tensorqtl-in.phenotypes.bed.gz
│   │   ├── trans-eqtl
│   │   │   ├── Lung___EUR.tensorqtl-in.20.covariates.tsv
│   │   │   ├── Lung___EUR.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── Monocyte___EUR.tensorqtl-in.15.covariates.tsv
│   │   │   ├── Monocyte___EUR.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── Nasal_epithelial___EUR.tensorqtl-in.20.covariates.tsv
│   │   │   ├── Nasal_epithelial___EUR.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── PBMC___AFR.tensorqtl-in.20.covariates.tsv
│   │   │   ├── PBMC___AFR.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── PBMC___EAS.tensorqtl-in.10.covariates.tsv
│   │   │   ├── PBMC___EAS.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── PBMC___EUR.tensorqtl-in.20.covariates.tsv
│   │   │   ├── PBMC___EUR.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── T_cell___EUR.tensorqtl-in.15.covariates.tsv
│   │   │   ├── T_cell___EUR.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── Whole_blood___AFR.tensorqtl-in.40.covariates.tsv
│   │   │   ├── Whole_blood___AFR.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── Whole_blood___EUR.tensorqtl-in.40.covariates.tsv
│   │   │   └── Whole_blood___EUR.tensorqtl-in.phenotypes.bed.gz
│   │   └── trans-sqtl
│   │       ├── Lung___EUR.tensorqtl-in.10.covariates.tsv
│   │       ├── Lung___EUR.tensorqtl-in.phenotypes.bed.gz
│   │       ├── Monocyte___EUR.tensorqtl-in.5.covariates.tsv
│   │       ├── Monocyte___EUR.tensorqtl-in.phenotypes.bed.gz
│   │       ├── Nasal_epithelial___EUR.tensorqtl-in.10.covariates.tsv
│   │       ├── Nasal_epithelial___EUR.tensorqtl-in.phenotypes.bed.gz
│   │       ├── PBMC___AFR.tensorqtl-in.5.covariates.tsv
│   │       ├── PBMC___AFR.tensorqtl-in.phenotypes.bed.gz
│   │       ├── PBMC___EAS.tensorqtl-in.5.covariates.tsv
│   │       ├── PBMC___EAS.tensorqtl-in.phenotypes.bed.gz
│   │       ├── PBMC___EUR.tensorqtl-in.10.covariates.tsv
│   │       ├── PBMC___EUR.tensorqtl-in.phenotypes.bed.gz
│   │       ├── T_cell___EUR.tensorqtl-in.5.covariates.tsv
│   │       ├── T_cell___EUR.tensorqtl-in.phenotypes.bed.gz
│   │       ├── Whole_blood___AFR.tensorqtl-in.10.covariates.tsv
│   │       ├── Whole_blood___AFR.tensorqtl-in.phenotypes.bed.gz
│   │       ├── Whole_blood___EUR.tensorqtl-in.10.covariates.tsv
│   │       └── Whole_blood___EUR.tensorqtl-in.phenotypes.bed.gz
│   ├── joint
│   │   ├── cis-eqtl
│   │   │   ├── Lung.tensorqtl-in.75.covariates.tsv
│   │   │   ├── Lung.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── Monocyte.tensorqtl-in.30.covariates.tsv
│   │   │   ├── Monocyte.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── Nasal_epithelial.tensorqtl-in.30.covariates.tsv
│   │   │   ├── Nasal_epithelial.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── PBMC.tensorqtl-in.30.covariates.tsv
│   │   │   ├── PBMC.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── T_cell.tensorqtl-in.30.covariates.tsv
│   │   │   ├── T_cell.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── Whole_blood.tensorqtl-in.100.covariates.tsv
│   │   │   └── Whole_blood.tensorqtl-in.phenotypes.bed.gz
│   │   ├── cis-sqtl
│   │   │   ├── Lung.tensorqtl-in.10.covariates.tsv
│   │   │   ├── Lung.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── Monocyte.tensorqtl-in.10.covariates.tsv
│   │   │   ├── Monocyte.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── Nasal_epithelial.tensorqtl-in.10.covariates.tsv
│   │   │   ├── Nasal_epithelial.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── PBMC.tensorqtl-in.10.covariates.tsv
│   │   │   ├── PBMC.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── T_cell.tensorqtl-in.10.covariates.tsv
│   │   │   ├── T_cell.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── Whole_blood.tensorqtl-in.10.covariates.tsv
│   │   │   └── Whole_blood.tensorqtl-in.phenotypes.bed.gz
│   │   ├── trans-eqtl
│   │   │   ├── Lung.tensorqtl-in.75.covariates.tsv
│   │   │   ├── Lung.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── Monocyte.tensorqtl-in.30.covariates.tsv
│   │   │   ├── Monocyte.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── Nasal_epithelial.tensorqtl-in.30.covariates.tsv
│   │   │   ├── Nasal_epithelial.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── PBMC.tensorqtl-in.30.covariates.tsv
│   │   │   ├── PBMC.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── T_cell.tensorqtl-in.30.covariates.tsv
│   │   │   ├── T_cell.tensorqtl-in.phenotypes.bed.gz
│   │   │   ├── Whole_blood___EUR.tensorqtl-in.40.covariates.tsv
│   │   │   ├── Whole_blood.tensorqtl-in.100.covariates.tsv
│   │   │   └── Whole_blood.tensorqtl-in.phenotypes.bed.gz
│   │   └── trans-sqtl
│   │       ├── Lung.tensorqtl-in.10.covariates.tsv
│   │       ├── Lung.tensorqtl-in.phenotypes.bed.gz
│   │       ├── Monocyte.tensorqtl-in.10.covariates.tsv
│   │       ├── Monocyte.tensorqtl-in.phenotypes.bed.gz
│   │       ├── Nasal_epithelial.tensorqtl-in.10.covariates.tsv
│   │       ├── Nasal_epithelial.tensorqtl-in.phenotypes.bed.gz
│   │       ├── PBMC.tensorqtl-in.10.covariates.tsv
│   │       ├── PBMC.tensorqtl-in.phenotypes.bed.gz
│   │       ├── T_cell.tensorqtl-in.10.covariates.tsv
│   │       ├── T_cell.tensorqtl-in.phenotypes.bed.gz
│   │       ├── Whole_blood.tensorqtl-in.10.covariates.tsv
│   │       └── Whole_blood.tensorqtl-in.phenotypes.bed.gz
│   └── saturation
│       ├── cis-eqtl
│       │   ├── Whole_blood_1000.tensorqtl-in.0.covariates.tsv
│       │   ├── Whole_blood_1000.tensorqtl-in.100.covariates.tsv
│       │   ├── Whole_blood_1000.tensorqtl-in.10.covariates.tsv
│       │   ├── Whole_blood_1000.tensorqtl-in.15.covariates.tsv
│       │   ├── Whole_blood_1000.tensorqtl-in.20.covariates.tsv
│       │   ├── Whole_blood_1000.tensorqtl-in.25.covariates.tsv
│       │   ├── Whole_blood_1000.tensorqtl-in.30.covariates.tsv
│       │   ├── Whole_blood_1000.tensorqtl-in.35.covariates.tsv
│       │   ├── Whole_blood_1000.tensorqtl-in.40.covariates.tsv
│       │   ├── Whole_blood_1000.tensorqtl-in.45.covariates.tsv
│       │   ├── Whole_blood_1000.tensorqtl-in.50.covariates.tsv
│       │   ├── Whole_blood_1000.tensorqtl-in.55.covariates.tsv
│       │   ├── Whole_blood_1000.tensorqtl-in.5.covariates.tsv
│       │   ├── Whole_blood_1000.tensorqtl-in.60.covariates.tsv
│       │   ├── Whole_blood_1000.tensorqtl-in.65.covariates.tsv
│       │   ├── Whole_blood_1000.tensorqtl-in.70.covariates.tsv
│       │   ├── Whole_blood_1000.tensorqtl-in.75.covariates.tsv
│       │   ├── Whole_blood_1000.tensorqtl-in.80.covariates.tsv
│       │   ├── Whole_blood_1000.tensorqtl-in.85.covariates.tsv
│       │   ├── Whole_blood_1000.tensorqtl-in.90.covariates.tsv
│       │   ├── Whole_blood_1000.tensorqtl-in.95.covariates.tsv
│       │   ├── Whole_blood_1000.tensorqtl-in.phenotypes.bed.gz
│       │   ├── Whole_blood_1500.tensorqtl-in.0.covariates.tsv
│       │   ├── Whole_blood_1500.tensorqtl-in.100.covariates.tsv
│       │   ├── Whole_blood_1500.tensorqtl-in.10.covariates.tsv
│       │   ├── Whole_blood_1500.tensorqtl-in.15.covariates.tsv
│       │   ├── Whole_blood_1500.tensorqtl-in.20.covariates.tsv
│       │   ├── Whole_blood_1500.tensorqtl-in.25.covariates.tsv
│       │   ├── Whole_blood_1500.tensorqtl-in.30.covariates.tsv
│       │   ├── Whole_blood_1500.tensorqtl-in.35.covariates.tsv
│       │   ├── Whole_blood_1500.tensorqtl-in.40.covariates.tsv
│       │   ├── Whole_blood_1500.tensorqtl-in.45.covariates.tsv
│       │   ├── Whole_blood_1500.tensorqtl-in.50.covariates.tsv
│       │   ├── Whole_blood_1500.tensorqtl-in.55.covariates.tsv
│       │   ├── Whole_blood_1500.tensorqtl-in.5.covariates.tsv
│       │   ├── Whole_blood_1500.tensorqtl-in.60.covariates.tsv
│       │   ├── Whole_blood_1500.tensorqtl-in.65.covariates.tsv
│       │   ├── Whole_blood_1500.tensorqtl-in.70.covariates.tsv
│       │   ├── Whole_blood_1500.tensorqtl-in.75.covariates.tsv
│       │   ├── Whole_blood_1500.tensorqtl-in.80.covariates.tsv
│       │   ├── Whole_blood_1500.tensorqtl-in.85.covariates.tsv
│       │   ├── Whole_blood_1500.tensorqtl-in.90.covariates.tsv
│       │   ├── Whole_blood_1500.tensorqtl-in.95.covariates.tsv
│       │   ├── Whole_blood_1500.tensorqtl-in.phenotypes.bed.gz
│       │   ├── Whole_blood_2000.tensorqtl-in.0.covariates.tsv
│       │   ├── Whole_blood_2000.tensorqtl-in.100.covariates.tsv
│       │   ├── Whole_blood_2000.tensorqtl-in.10.covariates.tsv
│       │   ├── Whole_blood_2000.tensorqtl-in.15.covariates.tsv
│       │   ├── Whole_blood_2000.tensorqtl-in.20.covariates.tsv
│       │   ├── Whole_blood_2000.tensorqtl-in.25.covariates.tsv
│       │   ├── Whole_blood_2000.tensorqtl-in.30.covariates.tsv
│       │   ├── Whole_blood_2000.tensorqtl-in.35.covariates.tsv
│       │   ├── Whole_blood_2000.tensorqtl-in.40.covariates.tsv
│       │   ├── Whole_blood_2000.tensorqtl-in.45.covariates.tsv
│       │   ├── Whole_blood_2000.tensorqtl-in.50.covariates.tsv
│       │   ├── Whole_blood_2000.tensorqtl-in.55.covariates.tsv
│       │   ├── Whole_blood_2000.tensorqtl-in.5.covariates.tsv
│       │   ├── Whole_blood_2000.tensorqtl-in.60.covariates.tsv
│       │   ├── Whole_blood_2000.tensorqtl-in.65.covariates.tsv
│       │   ├── Whole_blood_2000.tensorqtl-in.70.covariates.tsv
│       │   ├── Whole_blood_2000.tensorqtl-in.75.covariates.tsv
│       │   ├── Whole_blood_2000.tensorqtl-in.80.covariates.tsv
│       │   ├── Whole_blood_2000.tensorqtl-in.85.covariates.tsv
│       │   ├── Whole_blood_2000.tensorqtl-in.90.covariates.tsv
│       │   ├── Whole_blood_2000.tensorqtl-in.95.covariates.tsv
│       │   ├── Whole_blood_2000.tensorqtl-in.phenotypes.bed.gz
│       │   ├── Whole_blood_2500.tensorqtl-in.0.covariates.tsv
│       │   ├── Whole_blood_2500.tensorqtl-in.100.covariates.tsv
│       │   ├── Whole_blood_2500.tensorqtl-in.10.covariates.tsv
│       │   ├── Whole_blood_2500.tensorqtl-in.15.covariates.tsv
│       │   ├── Whole_blood_2500.tensorqtl-in.20.covariates.tsv
│       │   ├── Whole_blood_2500.tensorqtl-in.25.covariates.tsv
│       │   ├── Whole_blood_2500.tensorqtl-in.30.covariates.tsv
│       │   ├── Whole_blood_2500.tensorqtl-in.35.covariates.tsv
│       │   ├── Whole_blood_2500.tensorqtl-in.40.covariates.tsv
│       │   ├── Whole_blood_2500.tensorqtl-in.45.covariates.tsv
│       │   ├── Whole_blood_2500.tensorqtl-in.50.covariates.tsv
│       │   ├── Whole_blood_2500.tensorqtl-in.55.covariates.tsv
│       │   ├── Whole_blood_2500.tensorqtl-in.5.covariates.tsv
│       │   ├── Whole_blood_2500.tensorqtl-in.60.covariates.tsv
│       │   ├── Whole_blood_2500.tensorqtl-in.65.covariates.tsv
│       │   ├── Whole_blood_2500.tensorqtl-in.70.covariates.tsv
│       │   ├── Whole_blood_2500.tensorqtl-in.75.covariates.tsv
│       │   ├── Whole_blood_2500.tensorqtl-in.80.covariates.tsv
│       │   ├── Whole_blood_2500.tensorqtl-in.85.covariates.tsv
│       │   ├── Whole_blood_2500.tensorqtl-in.90.covariates.tsv
│       │   ├── Whole_blood_2500.tensorqtl-in.95.covariates.tsv
│       │   ├── Whole_blood_2500.tensorqtl-in.phenotypes.bed.gz
│       │   ├── Whole_blood_3000.tensorqtl-in.0.covariates.tsv
│       │   ├── Whole_blood_3000.tensorqtl-in.100.covariates.tsv
│       │   ├── Whole_blood_3000.tensorqtl-in.10.covariates.tsv
│       │   ├── Whole_blood_3000.tensorqtl-in.15.covariates.tsv
│       │   ├── Whole_blood_3000.tensorqtl-in.20.covariates.tsv
│       │   ├── Whole_blood_3000.tensorqtl-in.25.covariates.tsv
│       │   ├── Whole_blood_3000.tensorqtl-in.30.covariates.tsv
│       │   ├── Whole_blood_3000.tensorqtl-in.35.covariates.tsv
│       │   ├── Whole_blood_3000.tensorqtl-in.40.covariates.tsv
│       │   ├── Whole_blood_3000.tensorqtl-in.45.covariates.tsv
│       │   ├── Whole_blood_3000.tensorqtl-in.50.covariates.tsv
│       │   ├── Whole_blood_3000.tensorqtl-in.55.covariates.tsv
│       │   ├── Whole_blood_3000.tensorqtl-in.5.covariates.tsv
│       │   ├── Whole_blood_3000.tensorqtl-in.60.covariates.tsv
│       │   ├── Whole_blood_3000.tensorqtl-in.65.covariates.tsv
│       │   ├── Whole_blood_3000.tensorqtl-in.70.covariates.tsv
│       │   ├── Whole_blood_3000.tensorqtl-in.75.covariates.tsv
│       │   ├── Whole_blood_3000.tensorqtl-in.80.covariates.tsv
│       │   ├── Whole_blood_3000.tensorqtl-in.85.covariates.tsv
│       │   ├── Whole_blood_3000.tensorqtl-in.90.covariates.tsv
│       │   ├── Whole_blood_3000.tensorqtl-in.95.covariates.tsv
│       │   ├── Whole_blood_3000.tensorqtl-in.phenotypes.bed.gz
│       │   ├── Whole_blood_3500.tensorqtl-in.0.covariates.tsv
│       │   ├── Whole_blood_3500.tensorqtl-in.100.covariates.tsv
│       │   ├── Whole_blood_3500.tensorqtl-in.10.covariates.tsv
│       │   ├── Whole_blood_3500.tensorqtl-in.15.covariates.tsv
│       │   ├── Whole_blood_3500.tensorqtl-in.20.covariates.tsv
│       │   ├── Whole_blood_3500.tensorqtl-in.25.covariates.tsv
│       │   ├── Whole_blood_3500.tensorqtl-in.30.covariates.tsv
│       │   ├── Whole_blood_3500.tensorqtl-in.35.covariates.tsv
│       │   ├── Whole_blood_3500.tensorqtl-in.40.covariates.tsv
│       │   ├── Whole_blood_3500.tensorqtl-in.45.covariates.tsv
│       │   ├── Whole_blood_3500.tensorqtl-in.50.covariates.tsv
│       │   ├── Whole_blood_3500.tensorqtl-in.55.covariates.tsv
│       │   ├── Whole_blood_3500.tensorqtl-in.5.covariates.tsv
│       │   ├── Whole_blood_3500.tensorqtl-in.60.covariates.tsv
│       │   ├── Whole_blood_3500.tensorqtl-in.65.covariates.tsv
│       │   ├── Whole_blood_3500.tensorqtl-in.70.covariates.tsv
│       │   ├── Whole_blood_3500.tensorqtl-in.75.covariates.tsv
│       │   ├── Whole_blood_3500.tensorqtl-in.80.covariates.tsv
│       │   ├── Whole_blood_3500.tensorqtl-in.85.covariates.tsv
│       │   ├── Whole_blood_3500.tensorqtl-in.90.covariates.tsv
│       │   ├── Whole_blood_3500.tensorqtl-in.95.covariates.tsv
│       │   ├── Whole_blood_3500.tensorqtl-in.phenotypes.bed.gz
│       │   ├── Whole_blood_4000.tensorqtl-in.0.covariates.tsv
│       │   ├── Whole_blood_4000.tensorqtl-in.100.covariates.tsv
│       │   ├── Whole_blood_4000.tensorqtl-in.10.covariates.tsv
│       │   ├── Whole_blood_4000.tensorqtl-in.15.covariates.tsv
│       │   ├── Whole_blood_4000.tensorqtl-in.20.covariates.tsv
│       │   ├── Whole_blood_4000.tensorqtl-in.25.covariates.tsv
│       │   ├── Whole_blood_4000.tensorqtl-in.30.covariates.tsv
│       │   ├── Whole_blood_4000.tensorqtl-in.35.covariates.tsv
│       │   ├── Whole_blood_4000.tensorqtl-in.40.covariates.tsv
│       │   ├── Whole_blood_4000.tensorqtl-in.45.covariates.tsv
│       │   ├── Whole_blood_4000.tensorqtl-in.50.covariates.tsv
│       │   ├── Whole_blood_4000.tensorqtl-in.55.covariates.tsv
│       │   ├── Whole_blood_4000.tensorqtl-in.5.covariates.tsv
│       │   ├── Whole_blood_4000.tensorqtl-in.60.covariates.tsv
│       │   ├── Whole_blood_4000.tensorqtl-in.65.covariates.tsv
│       │   ├── Whole_blood_4000.tensorqtl-in.70.covariates.tsv
│       │   ├── Whole_blood_4000.tensorqtl-in.75.covariates.tsv
│       │   ├── Whole_blood_4000.tensorqtl-in.80.covariates.tsv
│       │   ├── Whole_blood_4000.tensorqtl-in.85.covariates.tsv
│       │   ├── Whole_blood_4000.tensorqtl-in.90.covariates.tsv
│       │   ├── Whole_blood_4000.tensorqtl-in.95.covariates.tsv
│       │   ├── Whole_blood_4000.tensorqtl-in.phenotypes.bed.gz
│       │   ├── Whole_blood_4500.tensorqtl-in.0.covariates.tsv
│       │   ├── Whole_blood_4500.tensorqtl-in.100.covariates.tsv
│       │   ├── Whole_blood_4500.tensorqtl-in.10.covariates.tsv
│       │   ├── Whole_blood_4500.tensorqtl-in.15.covariates.tsv
│       │   ├── Whole_blood_4500.tensorqtl-in.20.covariates.tsv
│       │   ├── Whole_blood_4500.tensorqtl-in.25.covariates.tsv
│       │   ├── Whole_blood_4500.tensorqtl-in.30.covariates.tsv
│       │   ├── Whole_blood_4500.tensorqtl-in.35.covariates.tsv
│       │   ├── Whole_blood_4500.tensorqtl-in.40.covariates.tsv
│       │   ├── Whole_blood_4500.tensorqtl-in.45.covariates.tsv
│       │   ├── Whole_blood_4500.tensorqtl-in.50.covariates.tsv
│       │   ├── Whole_blood_4500.tensorqtl-in.55.covariates.tsv
│       │   ├── Whole_blood_4500.tensorqtl-in.5.covariates.tsv
│       │   ├── Whole_blood_4500.tensorqtl-in.60.covariates.tsv
│       │   ├── Whole_blood_4500.tensorqtl-in.65.covariates.tsv
│       │   ├── Whole_blood_4500.tensorqtl-in.70.covariates.tsv
│       │   ├── Whole_blood_4500.tensorqtl-in.75.covariates.tsv
│       │   ├── Whole_blood_4500.tensorqtl-in.80.covariates.tsv
│       │   ├── Whole_blood_4500.tensorqtl-in.85.covariates.tsv
│       │   ├── Whole_blood_4500.tensorqtl-in.90.covariates.tsv
│       │   ├── Whole_blood_4500.tensorqtl-in.95.covariates.tsv
│       │   ├── Whole_blood_4500.tensorqtl-in.phenotypes.bed.gz
│       │   ├── Whole_blood_5000.tensorqtl-in.0.covariates.tsv
│       │   ├── Whole_blood_5000.tensorqtl-in.100.covariates.tsv
│       │   ├── Whole_blood_5000.tensorqtl-in.10.covariates.tsv
│       │   ├── Whole_blood_5000.tensorqtl-in.15.covariates.tsv
│       │   ├── Whole_blood_5000.tensorqtl-in.20.covariates.tsv
│       │   ├── Whole_blood_5000.tensorqtl-in.25.covariates.tsv
│       │   ├── Whole_blood_5000.tensorqtl-in.30.covariates.tsv
│       │   ├── Whole_blood_5000.tensorqtl-in.35.covariates.tsv
│       │   ├── Whole_blood_5000.tensorqtl-in.40.covariates.tsv
│       │   ├── Whole_blood_5000.tensorqtl-in.45.covariates.tsv
│       │   ├── Whole_blood_5000.tensorqtl-in.50.covariates.tsv
│       │   ├── Whole_blood_5000.tensorqtl-in.55.covariates.tsv
│       │   ├── Whole_blood_5000.tensorqtl-in.5.covariates.tsv
│       │   ├── Whole_blood_5000.tensorqtl-in.60.covariates.tsv
│       │   ├── Whole_blood_5000.tensorqtl-in.65.covariates.tsv
│       │   ├── Whole_blood_5000.tensorqtl-in.70.covariates.tsv
│       │   ├── Whole_blood_5000.tensorqtl-in.75.covariates.tsv
│       │   ├── Whole_blood_5000.tensorqtl-in.80.covariates.tsv
│       │   ├── Whole_blood_5000.tensorqtl-in.85.covariates.tsv
│       │   ├── Whole_blood_5000.tensorqtl-in.90.covariates.tsv
│       │   ├── Whole_blood_5000.tensorqtl-in.95.covariates.tsv
│       │   ├── Whole_blood_5000.tensorqtl-in.phenotypes.bed.gz
│       │   ├── Whole_blood_500.tensorqtl-in.0.covariates.tsv
│       │   ├── Whole_blood_500.tensorqtl-in.100.covariates.tsv
│       │   ├── Whole_blood_500.tensorqtl-in.10.covariates.tsv
│       │   ├── Whole_blood_500.tensorqtl-in.15.covariates.tsv
│       │   ├── Whole_blood_500.tensorqtl-in.20.covariates.tsv
│       │   ├── Whole_blood_500.tensorqtl-in.25.covariates.tsv
│       │   ├── Whole_blood_500.tensorqtl-in.30.covariates.tsv
│       │   ├── Whole_blood_500.tensorqtl-in.35.covariates.tsv
│       │   ├── Whole_blood_500.tensorqtl-in.40.covariates.tsv
│       │   ├── Whole_blood_500.tensorqtl-in.45.covariates.tsv
│       │   ├── Whole_blood_500.tensorqtl-in.50.covariates.tsv
│       │   ├── Whole_blood_500.tensorqtl-in.55.covariates.tsv
│       │   ├── Whole_blood_500.tensorqtl-in.5.covariates.tsv
│       │   ├── Whole_blood_500.tensorqtl-in.60.covariates.tsv
│       │   ├── Whole_blood_500.tensorqtl-in.65.covariates.tsv
│       │   ├── Whole_blood_500.tensorqtl-in.70.covariates.tsv
│       │   ├── Whole_blood_500.tensorqtl-in.75.covariates.tsv
│       │   ├── Whole_blood_500.tensorqtl-in.80.covariates.tsv
│       │   ├── Whole_blood_500.tensorqtl-in.85.covariates.tsv
│       │   ├── Whole_blood_500.tensorqtl-in.90.covariates.tsv
│       │   ├── Whole_blood_500.tensorqtl-in.95.covariates.tsv
│       │   ├── Whole_blood_500.tensorqtl-in.phenotypes.bed.gz
│       │   ├── Whole_blood_5500.tensorqtl-in.0.covariates.tsv
│       │   ├── Whole_blood_5500.tensorqtl-in.100.covariates.tsv
│       │   ├── Whole_blood_5500.tensorqtl-in.10.covariates.tsv
│       │   ├── Whole_blood_5500.tensorqtl-in.15.covariates.tsv
│       │   ├── Whole_blood_5500.tensorqtl-in.20.covariates.tsv
│       │   ├── Whole_blood_5500.tensorqtl-in.25.covariates.tsv
│       │   ├── Whole_blood_5500.tensorqtl-in.30.covariates.tsv
│       │   ├── Whole_blood_5500.tensorqtl-in.35.covariates.tsv
│       │   ├── Whole_blood_5500.tensorqtl-in.40.covariates.tsv
│       │   ├── Whole_blood_5500.tensorqtl-in.45.covariates.tsv
│       │   ├── Whole_blood_5500.tensorqtl-in.50.covariates.tsv
│       │   ├── Whole_blood_5500.tensorqtl-in.55.covariates.tsv
│       │   ├── Whole_blood_5500.tensorqtl-in.5.covariates.tsv
│       │   ├── Whole_blood_5500.tensorqtl-in.60.covariates.tsv
│       │   ├── Whole_blood_5500.tensorqtl-in.65.covariates.tsv
│       │   ├── Whole_blood_5500.tensorqtl-in.70.covariates.tsv
│       │   ├── Whole_blood_5500.tensorqtl-in.75.covariates.tsv
│       │   ├── Whole_blood_5500.tensorqtl-in.80.covariates.tsv
│       │   ├── Whole_blood_5500.tensorqtl-in.85.covariates.tsv
│       │   ├── Whole_blood_5500.tensorqtl-in.90.covariates.tsv
│       │   ├── Whole_blood_5500.tensorqtl-in.95.covariates.tsv
│       │   ├── Whole_blood_5500.tensorqtl-in.phenotypes.bed.gz
│       │   ├── Whole_blood_6000.tensorqtl-in.0.covariates.tsv
│       │   ├── Whole_blood_6000.tensorqtl-in.100.covariates.tsv
│       │   ├── Whole_blood_6000.tensorqtl-in.10.covariates.tsv
│       │   ├── Whole_blood_6000.tensorqtl-in.15.covariates.tsv
│       │   ├── Whole_blood_6000.tensorqtl-in.20.covariates.tsv
│       │   ├── Whole_blood_6000.tensorqtl-in.25.covariates.tsv
│       │   ├── Whole_blood_6000.tensorqtl-in.30.covariates.tsv
│       │   ├── Whole_blood_6000.tensorqtl-in.35.covariates.tsv
│       │   ├── Whole_blood_6000.tensorqtl-in.40.covariates.tsv
│       │   ├── Whole_blood_6000.tensorqtl-in.45.covariates.tsv
│       │   ├── Whole_blood_6000.tensorqtl-in.50.covariates.tsv
│       │   ├── Whole_blood_6000.tensorqtl-in.55.covariates.tsv
│       │   ├── Whole_blood_6000.tensorqtl-in.5.covariates.tsv
│       │   ├── Whole_blood_6000.tensorqtl-in.60.covariates.tsv
│       │   ├── Whole_blood_6000.tensorqtl-in.65.covariates.tsv
│       │   ├── Whole_blood_6000.tensorqtl-in.70.covariates.tsv
│       │   ├── Whole_blood_6000.tensorqtl-in.75.covariates.tsv
│       │   ├── Whole_blood_6000.tensorqtl-in.80.covariates.tsv
│       │   ├── Whole_blood_6000.tensorqtl-in.85.covariates.tsv
│       │   ├── Whole_blood_6000.tensorqtl-in.90.covariates.tsv
│       │   ├── Whole_blood_6000.tensorqtl-in.95.covariates.tsv
│       │   └── Whole_blood_6000.tensorqtl-in.phenotypes.bed.gz
│       ├── cis-sqtl
│       └── trans-eqtl
```


## Running

Wait for each step to complete before running the next one.

1. Subset genotypes to scan variants: `make scan-variant-vcf-files`
2. Subset genotypes to CS variants: `make cs-variant-vcf-files`
3. Clump trans variants: `make clump-trans-variants`
4. Build joint models: `make joint-model-with-cs-variants`
5. Allelic fold change (cis-eQTLs): `make allelic-fold-change-cis-eqtl`
6. Allelic fold change (trans-eQTLs): `make allelic-fold-change-trans-eqtl`
7. Calculate statistics used for creating control credible sets: `make precalculate-matching-statistics-cis`
8. Create control credible sets: `make control-credible-sets-cis`
9. Create trans control variants: `make control-snps-trans`
10. Compute ancestry allele counts for all variants in cis-eQTL scans: `make compute-ancestry-allele-counts-for-all-variants-in-cis-eqtl-scans-75-AMR-50`
11. Compute ancestry allele counts for all CS variants for top trans variants: `make compute-ancestry-allele-counts-for-cs-or-top-variants-75-AMR-50`
12. Make variant annotation matrices: `make variant-annotation-matrices`
13. Blood as proxy: `make blood-as-proxy`
15. Summarize cis-e/sQTL rare signals `make cis-rare-variants-summary`
16. Rank signals: `make rank-signals`
17. Attempt to replicate eQTLGen trans-eQTL signals: `make replicate-trans-eqtlgen`; when finished, run `make replication-of-eqtlgen-trans-eqtl` to plot results
18. Attempt to replicate DIRECT trans-eQTL signals: `make replicate-trans-direct`; when finished, run `make replication-of-direct-trans-eqtl` to plot results
19. Attempt to replicate GTEx trans-eQTL signals: `make replicate-trans-gtex`; when finished, run `make replication-of-gtex-trans-qtl` to examine results
20. Run the coloc commands (below)
14. Summarize cis-e/sQTL colocs: `make summarize-cis-eqtl-coloc-with-cis-sqtl`
15. Summarize PanUKBB colocs: `make summarize-panukbb-coloc`
15. Summarize ancestry counts: `make sample-ancestry-assignment`
20. Run functional enrichments: `make functional-enrichments`
21. trans-eQTL GO enrichments: `make trans-go-enrichment`
22. Whole blood trans enrichments in cis signals: `make trans-enrichment-in-cis`
23. Get LD buddies for eQTLGen and DIRECT signals: `make eqtlgen-and-direct-ld-buddies`

## Coloc

1. tensorQTL output (python objects) need to be converted to R objects: `make tensorqtl-to-coloc-in`
2. `make susie-json-joint susie-json-panukbb susie-json-ancestry-specific susie-json-saturation`
3. Stage the files: `make coloc-stage-joint coloc-stage-panukbb coloc-stage-ancestry-specific coloc-stage-saturation`
3. Launch TOPMed xQTL - PanUKBB GWAS colocs: `make coloc-panukbb-joint`
4. Postprocess TOPMed xQTL - PanUKBB GWAS colocs: `make postprocess-coloc-panukbb-joint`
5. Launch ancestry-specific TOPMed xQTL - PanUKBB GWAS colocs: `make coloc-panukbb-ancestry-specific`
6. Postprocess TOPMed xQTL - PanUKBB GWAS colocs: `make postprocess-coloc-panukbb-ancestry-specific`
7. Launch TOPMed saturation cis-eQTL - PanUKBB GWAS colocs: `make coloc-panukbb-saturation`
8. Postprocess TOPMed saturation cis-eQTL - PanUKBB GWAS colocs: `make postprocess-coloc-panukbb-saturation`
9. Launch xQTL-xQTL colocs (e.g., cis-trans colocs, cis-e vs cis-sQTL colocs): `make coloc-xqtl-joint`
10. Postprocess xQTL-xQTL colocs (e.g., cis-trans colocs, cis-e vs cis-sQTL colocs): `make postprocess-coloc-xqtl-joint`
