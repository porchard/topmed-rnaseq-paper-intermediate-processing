ROOT=/net/topmed11/working/porchard/rnaseq-2024-10-16-sample-update/manuscript-intermediate-processing
DATA=$(ROOT)/data
WORK=$(ROOT)/work
BIN=$(ROOT)/bin

SIF=singularity exec --bind /net /net/topmed3/working/porchard/rnaseq-paper-figures/general_20240516.sif

ANALYSIS=$(WORK)/$@

.PHONY: all

define NL


endef

singularity:
	singularity pull general.sif docker://porchard/general_jl:20231026105312

data: fasta-hg38 tensorqtl-in tensorqtl-in-new tensorqtl-out saturation-analysis ancestry-specific freeze-1.1RNA metadata gene-counts genotypes roadmap variant-sensitive-motif-scanning crossmap-gencode-v30 panukbb-finemapping direct ensembl-regulatory-build
	
fasta-hg38:
	mkdir -p $(DATA)/fasta/hg38
	cp /net/topmed10/working/porchard/rnaseq/data/fasta/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta $(DATA)/fasta/hg38/

roadmap:
	mkdir -p $(DATA)/$@
	$(foreach i,E029 E034 E062 E096 E114,wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/$(i)_15_coreMarks_hg38lift_mnemonics.bed.gz --directory-prefix $(DATA)/$@/$(NL))

ensembl-regulatory-build:
	mkdir -p $(DATA)/$@ && cd $(DATA)/$@ && wget https://ftp.ensembl.org/pub/current_regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20221007.gff.gz

crossmap-gencode-v30:
	mkdir -p $(DATA)/$@
	ln -s /net/topmed10/working/porchard/rnaseq/work/crossmap-gencode-v30/results/snp-mappability/snp_mappability_100mer_2mismatch.bed.gz $(DATA)/$@/
	ln -s /net/topmed10/working/porchard/rnaseq/work/crossmap-gencode-v30/results/snp-mappability/snp_mappability_100mer_2mismatch.bed.gz.tbi $(DATA)/$@/

direct:
	mkdir -p $(DATA)/$@
	# cd $(DATA)/$@ && wget https://zenodo.org/api/records/7521410/files-archive
	# description of tables: https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-023-40569-3/MediaObjects/41467_2023_40569_MOESM3_ESM.pdf
	# cd $(DATA)/$@ && wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-023-40569-3/MediaObjects/41467_2023_40569_MOESM4_ESM.txt -O Table-S1.txt
	# cd $(DATA)/$@ && wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-023-40569-3/MediaObjects/41467_2023_40569_MOESM6_ESM.txt -O Table-S3.txt
	# cd $(DATA)/$@ && wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-023-40569-3/MediaObjects/41467_2023_40569_MOESM10_ESM.txt -O Table-S7.txt
	cd $(DATA)/$@ && wget https://zenodo.org/records/7521410/files/Pvalues_nominal_trans_eQTLs_10e4_Genes_DIRECT.txt.gz
	

# TODO: add phenotype groups for sQTL scans
tensorqtl-in:
	mkdir -p $(DATA)/$@/{joint,saturation,ancestry-specific}/{cis-eqtl,cis-sqtl,trans-eqtl,trans-sqtl}
	rmdir $(DATA)/$@/saturation/trans-sqtl # did no such analysis
	# joint, cis-eQTL
	cp /net/topmed11/working/porchard/rnaseq-2024-10-16-sample-update/qtl-scan-input/work/tensorqtl-in/eqtl/results/tensorqtl-in/{Whole_blood,Lung,Nasal_epithelial}*.phenotypes.bed.gz $(DATA)/$@/joint/cis-eqtl/
	cp /net/topmed11/working/porchard/rnaseq-paper-intermediate-processing/data/tensorqtl-in-new/cis-eqtl/{PBMC,T_cell,Monocyte}*.phenotypes.bed.gz $(DATA)/$@/joint/cis-eqtl/
	$(foreach t,PBMC T_cell Monocyte,cp /net/topmed11/working/porchard/rnaseq-paper-intermediate-processing/data/tensorqtl-in-new/cis-eqtl/$(t).tensorqtl-in.30.covariates.tsv $(DATA)/$@/joint/cis-eqtl/$(NL))
	cp /net/topmed11/working/porchard/rnaseq-2024-10-16-sample-update/qtl-scan-input/work/tensorqtl-in/eqtl/results/tensorqtl-in/Nasal_epithelial.tensorqtl-in.30.covariates.tsv $(DATA)/$@/joint/cis-eqtl/
	cp /net/topmed11/working/porchard/rnaseq-2024-10-16-sample-update/qtl-scan-input/work/tensorqtl-in/eqtl/results/tensorqtl-in/Lung.tensorqtl-in.75.covariates.tsv $(DATA)/$@/joint/cis-eqtl/
	cp /net/topmed11/working/porchard/rnaseq-2024-10-16-sample-update/qtl-scan-input/work/tensorqtl-in/eqtl/results/tensorqtl-in/Whole_blood.tensorqtl-in.100.covariates.tsv $(DATA)/$@/joint/cis-eqtl/
	# joint, cis-sQTL
	cp /net/topmed11/working/porchard/rnaseq-2024-10-16-sample-update/qtl-scan-input/work/tensorqtl-in/sqtl/results/tensorqtl-in/{Whole_blood,Lung,Nasal_epithelial}*.phenotypes.bed.gz $(DATA)/$@/joint/cis-sqtl/
	cp /net/topmed11/working/porchard/rnaseq-2024-10-16-sample-update/qtl-scan-input/work/tensorqtl-in/sqtl/results/tensorqtl-in/{Whole_blood,Lung,Nasal_epithelial}.tensorqtl-in.10.covariates.tsv $(DATA)/$@/joint/cis-sqtl/
	cp /net/topmed11/working/porchard/rnaseq-paper-intermediate-processing/data/tensorqtl-in-new/cis-sqtl/{PBMC,T_cell,Monocyte}*.phenotypes.bed.gz $(DATA)/$@/joint/cis-sqtl/
	cp /net/topmed11/working/porchard/rnaseq-paper-intermediate-processing/data/tensorqtl-in-new/cis-sqtl/{PBMC,T_cell,Monocyte}.tensorqtl-in.10.covariates.tsv $(DATA)/$@/joint/cis-sqtl/
	# joint, trans-eQTL
	cp $(DATA)/$@/joint/cis-eqtl/* $(DATA)/$@/joint/trans-eqtl/
	grep $(foreach DROP,25 $(shell seq 51 100),-w -v -e phenotype_PC$(DROP)) /net/topmed11/working/porchard/rnaseq-2024-10-16-sample-update/qtl-scan-input/work/tensorqtl-in/eqtl/results/tensorqtl-in/Whole_blood.tensorqtl-in.100.covariates.tsv > $(DATA)/$@/joint/trans-eqtl/Whole_blood.tensorqtl-in.100.covariates.tsv
	# joint, trans-sQTL
	cp $(DATA)/$@/joint/cis-sqtl/* $(DATA)/$@/joint/trans-sqtl/
	# ancestry-specific, cis-eQTL
	cp /net/topmed11/working/porchard/rnaseq-2024-10-16-sample-update/qtl-scan-input/work/ancestry-specific/tensorqtl-in/eqtl/results/tensorqtl-in/{Whole_blood,Lung,Nasal_epithelial}*.phenotypes.bed.gz $(DATA)/$@/ancestry-specific/cis-eqtl/
	cp /net/topmed11/working/porchard/rnaseq-paper-intermediate-processing/data/ancestry-specific/tensorqtl-in/eqtl/results/tensorqtl-in/{PBMC,Monocyte,T_cell}*.phenotypes.bed.gz $(DATA)/$@/ancestry-specific/cis-eqtl/
	$(foreach t,PBMC___EAS,cp /net/topmed11/working/porchard/rnaseq-paper-intermediate-processing/data/ancestry-specific/tensorqtl-in/eqtl/results/tensorqtl-in/$(t).tensorqtl-in.10.covariates.tsv $(DATA)/$@/ancestry-specific/cis-eqtl/$(NL))
	$(foreach t,Monocyte___EUR T_cell___EUR,cp /net/topmed11/working/porchard/rnaseq-paper-intermediate-processing/data/ancestry-specific/tensorqtl-in/eqtl/results/tensorqtl-in/$(t).tensorqtl-in.15.covariates.tsv $(DATA)/$@/ancestry-specific/cis-eqtl/$(NL))
	$(foreach t,PBMC___AFR PBMC___EUR,cp /net/topmed11/working/porchard/rnaseq-paper-intermediate-processing/data/ancestry-specific/tensorqtl-in/eqtl/results/tensorqtl-in/$(t).tensorqtl-in.20.covariates.tsv $(DATA)/$@/ancestry-specific/cis-eqtl/$(NL))
	$(foreach t,Nasal_epithelial___EUR Lung___EUR,cp /net/topmed11/working/porchard/rnaseq-2024-10-16-sample-update/qtl-scan-input/work/ancestry-specific/tensorqtl-in/eqtl/results/tensorqtl-in/$(t).tensorqtl-in.20.covariates.tsv $(DATA)/$@/ancestry-specific/cis-eqtl/$(NL))
	$(foreach t,Whole_blood___AFR Whole_blood___EUR,cp /net/topmed11/working/porchard/rnaseq-2024-10-16-sample-update/qtl-scan-input/work/ancestry-specific/tensorqtl-in/eqtl/results/tensorqtl-in/$(t).tensorqtl-in.40.covariates.tsv $(DATA)/$@/ancestry-specific/cis-eqtl/$(NL))
	# ancestry-specific, cis-sQTL
	cp /net/topmed11/working/porchard/rnaseq-2024-10-16-sample-update/qtl-scan-input/work/ancestry-specific/tensorqtl-in/sqtl/results/tensorqtl-in/{Whole_blood,Lung,Nasal_epithelial}*.phenotypes.bed.gz $(DATA)/$@/ancestry-specific/cis-sqtl/
	cp /net/topmed11/working/porchard/rnaseq-paper-intermediate-processing/data/ancestry-specific/tensorqtl-in/sqtl/results/tensorqtl-in/{PBMC,Monocyte,T_cell}*.phenotypes.bed.gz $(DATA)/$@/ancestry-specific/cis-sqtl/
	$(foreach t,PBMC___AFR PBMC___EAS Monocyte___EUR T_cell___EUR,cp /net/topmed11/working/porchard/rnaseq-paper-intermediate-processing/data/ancestry-specific/tensorqtl-in/sqtl/results/tensorqtl-in/$(t).tensorqtl-in.5.covariates.tsv $(DATA)/$@/ancestry-specific/cis-sqtl/$(NL))
	$(foreach t,PBMC___EUR,cp /net/topmed11/working/porchard/rnaseq-paper-intermediate-processing/data/ancestry-specific/tensorqtl-in/sqtl/results/tensorqtl-in/$(t).tensorqtl-in.10.covariates.tsv $(DATA)/$@/ancestry-specific/cis-sqtl/$(NL))
	$(foreach t,Lung___EUR Nasal_epithelial___EUR Whole_blood___AFR Whole_blood___EUR,cp /net/topmed11/working/porchard/rnaseq-2024-10-16-sample-update/qtl-scan-input/work/ancestry-specific/tensorqtl-in/sqtl/results/tensorqtl-in/$(t).tensorqtl-in.10.covariates.tsv $(DATA)/$@/ancestry-specific/cis-sqtl/$(NL))
	# ancestry-specific, trans-eQTL
	cp $(DATA)/$@/ancestry-specific/cis-eqtl/* $(DATA)/$@/ancestry-specific/trans-eqtl/
	grep $(foreach DROP,19 22 27,-w -v -e phenotype_PC$(DROP)) /net/topmed11/working/porchard/rnaseq-2024-10-16-sample-update/qtl-scan-input/work/ancestry-specific/tensorqtl-in/eqtl/results/tensorqtl-in/Whole_blood___EUR.tensorqtl-in.40.covariates.tsv > $(DATA)/$@/joint/trans-eqtl/Whole_blood___EUR.tensorqtl-in.40.covariates.tsv
	# ancestry-specific, trans-sQTL
	cp $(DATA)/$@/ancestry-specific/cis-sqtl/* $(DATA)/$@/ancestry-specific/trans-sqtl/
	# saturation, cis-eQTL
	cp /net/topmed11/working/porchard/rnaseq-2024-10-16-sample-update/qtl-scan-input/work/saturation/tensorqtl-in/eqtl/results/tensorqtl-in/* $(DATA)/$@/saturation/cis-eqtl/

scan-results:
	mkdir -p $(DATA)
	ln -s /net/topmed11/working/porchard/rnaseq-2024-10-16-sample-update/qtl-scan-output $(DATA)/$@

ancestry:
	mkdir -p $(DATA)/$@
	cp /net/snowwhite/home/porchard/github/topmed-manuscript/topmed_manuscript/data/ancestries.txt $(DATA)/$@/

metadata:
	mkdir -p $(DATA)/$@
	cp /net/topmed11/working/porchard/rnaseq-2024-10-16-sample-update/qtl-scan-input/data/metadata/metadata.txt $(DATA)/$@/metadata.big.txt
	cut-name cohort,sequencing_center,tor,tissue,inferred_sex,wgs,used_for_scan $(DATA)/$@/metadata.big.txt > $(DATA)/$@/metadata.tm.txt

scan-samples:
	mkdir -p $(DATA)/$@/joint && cp /net/topmed11/working/porchard/rnaseq-2024-10-16-sample-update/qtl-scan-input/data/scan-samples/* $(DATA)/$@/joint/

gene-counts:
	mkdir -p $(DATA)/$@
	cp /net/topmed10/working/porchard/rnaseq/data/wolverine/all-cohorts.gene_counts.hdf5 $(DATA)/$@/
	cp /net/topmed10/working/porchard/rnaseq/data/wolverine/all-cohorts.gene_tpm.hdf5 $(DATA)/$@/

gtf:
	mkdir -p $(DATA)/$@
	cp /net/topmed10/working/porchard/rnaseq/data/gtf/gencode.v30.GRCh38.ERCC.genes.collapsed_only.gtf.gz $(DATA)/$@/

genotypes:
	mkdir -p $(DATA)/$@/
	ln -s /net/topmed10/working/porchard/rnaseq/work/subset-topmed-bcf/freeze-alpha/results/vcfs-updated-ids-pass-filter $(DATA)/$@/
	ln -s /net/topmed10/working/porchard/rnaseq/work/subset-topmed-bcf/freeze-alpha/results/plink-pass-filter $(DATA)/$@/

genotype-pca:
	mkdir -p $(DATA)/$@
	cp /net/topmed10/working/porchard/rnaseq/work/genotype-pca/freeze-beta/unrelated/results/pca/PC-scores.txt $(DATA)/$@/
	cp /net/topmed10/working/porchard/rnaseq/work/genotype-pca/freeze-beta/unrelated/results/pca/PC-variance-explained.txt $(DATA)/$@/

sites:
	mkdir -p $(DATA)/$@
	ln -s /net/topmed2/working/gt-release/exchange-area/freeze.9b/sites/* $(DATA)/$@/

variant-sensitive-motif-scanning:
	mkdir -p $(DATA)
	ln -s /net/topmed11/working/porchard/variant-sensitive-motif-scanning $(DATA)/


panukbb-finemapping:
	ln -s /net/topmed11/working/porchard/panukbb-finemapping $(DATA)/$@


scan-variant-vcf-files: ANALYSIS=$(WORK)/genotype-subsets/scan-variants
scan-variant-vcf-files: TISSUES=Lung Nasal_epithelial T_cell PBMC Monocyte Whole_blood
scan-variant-vcf-files:
	mkdir -p $(ANALYSIS)/data
	$(foreach t,$(TISSUES),zcat $(DATA)/tensorqtl-in/joint/cis-eqtl/$(t).tensorqtl-in.phenotypes.bed.gz | awk 'NR==1' | perl -pe 's/\t/\n/g' | grep NWD > $(ANALYSIS)/data/$(t).samples.txt$(NL))
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --samples_glob '$(ANALYSIS)/data/*.samples.txt' --vcf_glob '$(DATA)/genotypes/vcfs-updated-ids-pass-filter/*.vcf.gz' --vcf_index_glob '$(DATA)/genotypes/vcfs-updated-ids-pass-filter/*.vcf.gz.tbi' $(ROOT)/scan-variant-vcf-files.nf &

cs-variant-vcf-files: ANALYSIS=$(WORK)/genotype-subsets/cs-variants
cs-variant-vcf-files:
	mkdir -p $(ANALYSIS)
	cat $(DATA)/scan-results/joint/cis-eqtl/susie/maf*/*.cs.txt | cut -f2 | grep -v variant_id | sort | uniq > $(ANALYSIS)/variants.txt
	cat $(DATA)/scan-results/joint/cis-sqtl/susie/maf*/postprocessed/*cs.txt | cut -f2 | grep -v variant_id | sort | uniq >> $(ANALYSIS)/variants.txt
	cat $(DATA)/scan-results/joint/trans-eqtl/maf005/trans-susie/*.cs.txt | cut -f2 | grep -v variant_id | sort | uniq >> $(ANALYSIS)/variants.txt
	cat $(DATA)/scan-results/joint/trans-sqtl/maf005/trans-susie/*.cs.txt | cut -f2 | grep -v variant_id | sort | uniq >> $(ANALYSIS)/variants.txt
	cat $(ANALYSIS)/variants.txt | sort | uniq | perl -pe 's/_/\t/g' | awk '{print($$1, $$2-1, $$2)}' | perl -pe 's/ /\t/g' > $(ANALYSIS)/variants.bed
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --regions $(ANALYSIS)/variants.bed --vcf_glob '$(DATA)/genotypes/vcfs-updated-ids-pass-filter/*.vcf.gz' --vcf_index_glob '$(DATA)/genotypes/vcfs-updated-ids-pass-filter/*.vcf.gz.tbi' $(ROOT)/cs-variant-vcf-files.nf &

# TODO: not in README
scan-variant-info:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --scan_variant_vcf_glob '$(WORK)/genotype-subsets/scan-variants/results/vcfs-by-chrom/*.vcf.gz' --sites_vcf_glob '/net/topmed2/working/gt-release/exchange-area/freeze.9b/sites/*bcf' --sites_vcf_index_glob '/net/topmed2/working/gt-release/exchange-area/freeze.9b/sites/*bcf.csi' $(ROOT)/scan-variant-site-files.nf &

clump-trans-variants:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && singularity exec --bind /net docker://porchard/general:20241111 python $(BIN)/clump-trans-signals.py '$(DATA)/scan-results/joint/trans-eqtl/maf005/trans-top/*.top.txt' '$(DATA)/scan-results/joint/trans-sqtl/maf005/trans-top/*.top.txt' $(DATA)/genotypes/vcfs-updated-ids-pass-filter $(DATA)/metadata/metadata.tm.txt

joint-model-with-cs-variants: ANALYSIS=$(WORK)/joint-model-with-cs-variants
joint-model-with-cs-variants: TISSUES=Lung Nasal_epithelial T_cell PBMC Monocyte Whole_blood
joint-model-with-cs-variants:
	mkdir -p $(ANALYSIS)/data/{covariates,phenotypes,cs}
	$(foreach t,$(TISSUES),ln -s $(DATA)/tensorqtl-in/joint/cis-sqtl/$(t).tensorqtl-in.phenotypes.bed.gz $(ANALYSIS)/data/phenotypes/$(t).cissqtl.phenotypes.bed.gz$(NL))
	$(foreach t,$(TISSUES),ln -s $(DATA)/tensorqtl-in/joint/cis-sqtl/$(t).tensorqtl-in.10.covariates.tsv $(ANALYSIS)/data/covariates/$(t).cissqtl.covariates.tsv$(NL))
	$(foreach t,$(TISSUES),ln -s $(DATA)/tensorqtl-in/joint/trans-sqtl/$(t).tensorqtl-in.phenotypes.bed.gz $(ANALYSIS)/data/phenotypes/$(t).transsqtl.phenotypes.bed.gz$(NL))
	$(foreach t,$(TISSUES),ln -s $(DATA)/tensorqtl-in/joint/trans-sqtl/$(t).tensorqtl-in.10.covariates.tsv $(ANALYSIS)/data/covariates/$(t).transsqtl.covariates.tsv$(NL))
	$(foreach t,$(TISSUES),ln -s $(DATA)/tensorqtl-in/joint/cis-eqtl/$(t).tensorqtl-in.phenotypes.bed.gz $(ANALYSIS)/data/phenotypes/$(t).ciseqtl.phenotypes.bed.gz$(NL))
	$(foreach t,$(TISSUES),ln -s $(DATA)/tensorqtl-in/joint/cis-eqtl/$(t).tensorqtl-in.*.covariates.tsv $(ANALYSIS)/data/covariates/$(t).ciseqtl.covariates.tsv$(NL))
	$(foreach t,$(TISSUES),ln -s $(DATA)/tensorqtl-in/joint/trans-eqtl/$(t).tensorqtl-in.phenotypes.bed.gz $(ANALYSIS)/data/phenotypes/$(t).transeqtl.phenotypes.bed.gz$(NL))
	$(foreach t,$(TISSUES),ln -s $(DATA)/tensorqtl-in/joint/trans-eqtl/$(t).tensorqtl-in.*.covariates.tsv $(ANALYSIS)/data/covariates/$(t).transeqtl.covariates.tsv$(NL))
	$(foreach t,$(TISSUES),ln -s $(DATA)/scan-results/joint/cis-eqtl/susie/maf001/$(t).*.cs.txt $(ANALYSIS)/data/cs/$(t).maf001.ciseqtl.cs.txt$(NL))
	ln -s $(DATA)/scan-results/joint/cis-eqtl/susie/maf0001/Whole_blood.100.cs.txt $(ANALYSIS)/data/cs/Whole_blood.maf0001.ciseqtl.cs.txt
	$(foreach t,$(TISSUES),ln -s $(DATA)/scan-results/joint/cis-sqtl/susie/maf001/postprocessed/$(t).by-gene.per-intron-cs.txt $(ANALYSIS)/data/cs/$(t).maf001.cissqtl.cs.txt$(NL))
	ln -s $(DATA)/scan-results/joint/cis-sqtl/susie/maf0001/postprocessed/Whole_blood.by-gene.per-intron-cs.txt $(ANALYSIS)/data/cs/Whole_blood.maf0001.cissqtl.cs.txt
	$(foreach t,Lung PBMC Whole_blood,ln -s $(DATA)/scan-results/joint/trans-sqtl/maf005/trans-susie/$(t).cs.txt $(ANALYSIS)/data/cs/$(t).maf005.transsqtl.cs.txt$(NL))
	$(foreach t,$(TISSUES),ln -s $(DATA)/scan-results/joint/trans-eqtl/maf005/trans-susie/$(t).cs.txt $(ANALYSIS)/data/cs/$(t).maf005.transeqtl.cs.txt$(NL))
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --vcf_path '$(WORK)/genotype-subsets/cs-variants/results/vcfs-by-chrom' --covariates_glob '$(ANALYSIS)/data/covariates/*' --phenotypes_glob '$(ANALYSIS)/data/phenotypes/*' --credible_sets_glob '$(ANALYSIS)/data/cs/*' $(ROOT)/joint-models-using-cs-variants.nf &

allelic-fold-change-cis-eqtl: ANALYSIS=$(WORK)/allelic-fold-change/cis-eqtl
allelic-fold-change-cis-eqtl: TISSUES=Lung Nasal_epithelial T_cell PBMC Monocyte Whole_blood
allelic-fold-change-cis-eqtl:
	mkdir -p $(ANALYSIS)/data/{covariates,samples,susie}
	$(foreach t,$(TISSUES),ln -s $(DATA)/scan-results/joint/cis-eqtl/susie/maf001/$(t).*.cs.txt $(ANALYSIS)/data/susie/$(t).maf001.cs.txt$(NL))
	ln -s $(DATA)/scan-results/joint/cis-eqtl/susie/maf0001/Whole_blood.100.cs.txt $(ANALYSIS)/data/susie/Whole_blood.maf0001.cs.txt
	cp $(DATA)/tensorqtl-in/joint/cis-eqtl/*.tensorqtl-in.*.covariates.tsv $(ANALYSIS)/data/covariates/
	$(foreach t,$(TISSUES),cp $(DATA)/scan-samples/samples-to-use-for-scan.$(t).tsv $(ANALYSIS)/data/samples/$(t).samples.txt$(NL))
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --plink_glob '$(DATA)/genotypes/plink-pass-filter/*' --metadata $(DATA)/metadata/metadata.tm.txt --rna_counts $(DATA)/gene-counts/all-cohorts.gene_counts.hdf5 --samples_glob '$(ANALYSIS)/data/samples/*' --covariates_glob '$(ANALYSIS)/data/covariates/*' --credible_sets_glob '$(ANALYSIS)/data/susie/*' $(ROOT)/allelic-fold-change.nf &

allelic-fold-change-trans-eqtl: ANALYSIS=$(WORK)/allelic-fold-change/trans-eqtl
allelic-fold-change-trans-eqtl: TISSUES=Lung Nasal_epithelial T_cell PBMC Monocyte Whole_blood
allelic-fold-change-trans-eqtl:
	mkdir -p $(ANALYSIS)/data/{covariates,samples,gene-variant-pairs}
	$(foreach t,$(TISSUES),ln -s $(DATA)/tensorqtl-in/joint/trans-eqtl/$(t).tensorqtl-in.*.covariates.tsv $(ANALYSIS)/data/covariates/$(t).transeqtl.covariates.tsv$(NL))
	$(foreach t,$(TISSUES),cp $(DATA)/scan-samples/samples-to-use-for-scan.$(t).tsv $(ANALYSIS)/data/samples/$(t).samples.txt$(NL))
	$(foreach t,$(TISSUES),grep -e $(t) -e variant_id $(WORK)/clump-trans-variants/clump-trans-signals.significant-trans-eqtl-clumped.tsv | cut -f1,2 | awk '{print($$2, $$1)}' | perl -pe 's/ /\t/' > $(ANALYSIS)/data/gene-variant-pairs/$(t).gene-variant-pairs.txt$(NL))
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --plink_glob '$(DATA)/genotypes/plink-pass-filter/*' --metadata $(DATA)/metadata/metadata.tm.txt --rna_counts $(DATA)/gene-counts/all-cohorts.gene_counts.hdf5 --samples_glob '$(ANALYSIS)/data/samples/*' --covariates_glob '$(ANALYSIS)/data/covariates/*' --gene_variant_pairs_glob '$(ANALYSIS)/data/gene-variant-pairs/*' $(ROOT)/allelic-fold-change-trans.nf &

precalculate-matching-statistics-cis: precalculate-matching-statistics-cis-maf001 precalculate-matching-statistics-cis-maf0001

precalculate-matching-statistics-cis-maf001: ANALYSIS=$(WORK)/control-credible-sets/matching-statistics/cis/maf001
precalculate-matching-statistics-cis-maf001: TISSUES=Lung Nasal_epithelial T_cell PBMC Monocyte Whole_blood
precalculate-matching-statistics-cis-maf001:
	mkdir -p $(ANALYSIS)/samples-wgs
	$(foreach t,$(TISSUES),$(BIN)/cut-name wgs,tor $(DATA)/metadata/metadata.tm.txt | grep -f $(DATA)/scan-samples/samples-to-use-for-scan.$(t).tsv | cut -f1 > $(ANALYSIS)/samples-wgs/$(t).samples.txt$(NL))
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --min_maf '0.01' --plink_glob '$(DATA)/genotypes/plink-pass-filter/*' --sample_glob '$(ANALYSIS)/samples-wgs/*.samples.txt' $(ROOT)/precalculate-matching-statistics-cis.nf &

precalculate-matching-statistics-cis-maf0001: ANALYSIS=$(WORK)/control-credible-sets/matching-statistics/cis/maf0001
precalculate-matching-statistics-cis-maf0001: TISSUES=Whole_blood
precalculate-matching-statistics-cis-maf0001:
	mkdir -p $(ANALYSIS)/samples-wgs
	$(foreach t,$(TISSUES),$(BIN)/cut-name wgs,tor $(DATA)/metadata/metadata.tm.txt | grep -f $(DATA)/scan-samples/samples-to-use-for-scan.$(t).tsv | cut -f1 > $(ANALYSIS)/samples-wgs/$(t).samples.txt$(NL))
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --min_maf '0.001' --plink_glob '$(DATA)/genotypes/plink-pass-filter/*' --sample_glob '$(ANALYSIS)/samples-wgs/*.samples.txt' $(ROOT)/precalculate-matching-statistics-cis.nf &

control-credible-sets-cis: control-credible-sets-cis-eqtl-maf001 control-credible-sets-cis-eqtl-maf0001 control-credible-sets-cis-sqtl-maf001 control-credible-sets-cis-sqtl-maf0001

control-credible-sets-cis-eqtl-maf001: ANALYSIS=$(WORK)/control-credible-sets/cis-eqtl/maf001
control-credible-sets-cis-eqtl-maf001:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --susie_glob '$(DATA)/scan-results/joint/cis-eqtl/susie/maf001/*.cs.txt' --ld_glob '$(WORK)/control-credible-sets/matching-statistics/cis/maf001/results/maf-and-ld/*.ld.txt' --maf_glob '$(WORK)/control-credible-sets/matching-statistics/cis/maf001/results/maf-and-ld/*.maf.txt' --phenotype_glob '$(DATA)/tensorqtl-in/joint/cis-eqtl/*.phenotypes.bed.gz' $(ROOT)/control-credible-sets-cis.nf &

control-credible-sets-cis-eqtl-maf0001: ANALYSIS=$(WORK)/control-credible-sets/cis-eqtl/maf0001
control-credible-sets-cis-eqtl-maf0001:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --susie_glob '$(DATA)/scan-results/joint/cis-eqtl/susie/maf0001/*.cs.txt' --ld_glob '$(WORK)/control-credible-sets/matching-statistics/cis/maf0001/results/maf-and-ld/*.ld.txt' --maf_glob '$(WORK)/control-credible-sets/matching-statistics/cis/maf0001/results/maf-and-ld/*.maf.txt' --phenotype_glob '$(DATA)/tensorqtl-in/joint/cis-eqtl/*.phenotypes.bed.gz' $(ROOT)/control-credible-sets-cis.nf &

control-credible-sets-cis-sqtl-maf001: ANALYSIS=$(WORK)/control-credible-sets/cis-sqtl/maf001
control-credible-sets-cis-sqtl-maf001:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --susie_glob '$(DATA)/scan-results/joint/cis-sqtl/susie/maf001/postprocessed/*.by-gene.cs.txt' --ld_glob '$(WORK)/control-credible-sets/matching-statistics/cis/maf001/results/maf-and-ld/*.ld.txt' --maf_glob '$(WORK)/control-credible-sets/matching-statistics/cis/maf001/results/maf-and-ld/*.maf.txt' --phenotype_glob '$(DATA)/tensorqtl-in/joint/cis-sqtl/*.phenotypes.bed.gz' $(ROOT)/control-credible-sets-cis.nf &

control-credible-sets-cis-sqtl-maf0001: ANALYSIS=$(WORK)/control-credible-sets/cis-sqtl/maf0001
control-credible-sets-cis-sqtl-maf0001:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --susie_glob '$(DATA)/scan-results/joint/cis-sqtl/susie/maf0001/postprocessed/*.by-gene.cs.txt' --ld_glob '$(WORK)/control-credible-sets/matching-statistics/cis/maf0001/results/maf-and-ld/*.ld.txt' --maf_glob '$(WORK)/control-credible-sets/matching-statistics/cis/maf0001/results/maf-and-ld/*.maf.txt' --phenotype_glob '$(DATA)/tensorqtl-in/joint/cis-sqtl/*.phenotypes.bed.gz' $(ROOT)/control-credible-sets-cis.nf &

compute-ancestry-allele-counts-for-all-variants-in-cis-eqtl-scans-75-AMR-50: ANALYSIS=$(WORK)/ancestry-allele-counts-75-AMR-50/cis-eqtl-scans
compute-ancestry-allele-counts-for-all-variants-in-cis-eqtl-scans-75-AMR-50: TISSUES=Lung Nasal_epithelial T_cell PBMC Monocyte Whole_blood
compute-ancestry-allele-counts-for-all-variants-in-cis-eqtl-scans-75-AMR-50:
	mkdir -p $(ANALYSIS)/data/nominal
	$(foreach t,$(TISSUES),ln -s $(DATA)/scan-results/joint/cis-eqtl/nominal/maf001/$(t).tsv.gz $(ANALYSIS)/data/nominal/$(t).maf001.tsv.gz$(NL))
	ln -s $(DATA)/scan-results/joint/cis-eqtl/nominal/maf0001/Whole_blood.tsv.gz $(ANALYSIS)/data/nominal/Whole_blood.maf0001.tsv.gz
	mkdir -p $(ANALYSIS)/data/samples/50
	cd $(ANALYSIS)/data/samples && python $(BIN)/ancestry-frac-to-assignment.py --threshold 0.75 $(DATA)/metadata/metadata.tm.txt $(DATA)/ancestry/ancestries.txt
	cd $(ANALYSIS)/data/samples/50 && python $(BIN)/ancestry-frac-to-assignment.py --threshold 0.50 $(DATA)/metadata/metadata.tm.txt $(DATA)/ancestry/ancestries.txt
	mv $(ANALYSIS)/data/samples/50/*.AMR.samples.txt $(ANALYSIS)/data/samples/
	rm -r $(ANALYSIS)/data/samples/50
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --nominal_glob '$(ANALYSIS)/data/nominal/*.tsv.gz' --samples_glob '$(ANALYSIS)/data/samples/*' --plink_glob '$(DATA)/genotypes/plink-pass-filter/*' $(ROOT)/compute-ancestry-allele-counts-plink.nf &

compute-ancestry-allele-counts-for-cs-or-top-variants-75-AMR-50: ANALYSIS=$(WORK)/ancestry-allele-counts-75-AMR-50/cs-variants
compute-ancestry-allele-counts-for-cs-or-top-variants-75-AMR-50:
	# use all CS variants, or trans-e/sVariants
	mkdir -p $(ANALYSIS)/data/samples/50
	cd $(ANALYSIS)/data/samples && python $(BIN)/ancestry-frac-to-assignment.py --threshold 0.75 $(DATA)/metadata/metadata.tm.txt $(DATA)/ancestry/ancestries.txt
	cd $(ANALYSIS)/data/samples/50 && python $(BIN)/ancestry-frac-to-assignment.py --threshold 0.50 $(DATA)/metadata/metadata.tm.txt $(DATA)/ancestry/ancestries.txt
	mv $(ANALYSIS)/data/samples/50/*.AMR.samples.txt $(ANALYSIS)/data/samples/
	rm -r $(ANALYSIS)/data/samples/50
	cat $(DATA)/scan-results/joint/cis-eqtl/susie/maf00*/*.cs.txt | cut -f2 | grep -v variant_id > $(ANALYSIS)/variants.txt
	cat $(DATA)/scan-results/joint/cis-sqtl/susie/maf00*/postprocessed/*.per-intron-cs.txt | cut -f2 | grep -v variant_id >> $(ANALYSIS)/variants.txt
	cat $(DATA)/scan-results/joint/trans-*qtl/maf005/trans-top/* | cut -f1 | grep -v variant_id >> $(ANALYSIS)/variants.txt
	cat $(DATA)/scan-results/joint/trans-*qtl/maf005/trans-susie/*.cs.txt | cut -f2 | grep -v variant_id >> $(ANALYSIS)/variants.txt
	cat $(ANALYSIS)/variants.txt | sort | uniq > $(ANALYSIS)/variants.uniq.txt
	mv $(ANALYSIS)/variants.uniq.txt $(ANALYSIS)/variants.txt
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --variants $(ANALYSIS)/variants.txt --samples_glob '$(ANALYSIS)/data/samples/*' --plink_glob '$(DATA)/genotypes/plink-pass-filter/*' $(ROOT)/compute-ancestry-allele-counts-plink-selected-variants.nf &

control-snps-trans: control-snps-trans-eqtl control-snps-trans-sqtl

control-snps-trans-eqtl: ANALYSIS=$(WORK)/control-snps/trans-eqtl/maf005
control-snps-trans-eqtl: TISSUES=Lung Nasal_epithelial T_cell PBMC Monocyte Whole_blood
control-snps-trans-eqtl:
	mkdir -p $(ANALYSIS)/samples-wgs
	mkdir -p $(ANALYSIS)/top-variants
	$(foreach t,$(TISSUES),$(BIN)/cut-name wgs,tor $(DATA)/metadata/metadata.tm.txt | grep -f $(DATA)/scan-samples/samples-to-use-for-scan.$(t).tsv | cut -f1 > $(ANALYSIS)/samples-wgs/$(t).samples.txt$(NL))
	$(foreach t,$(TISSUES),$(BIN)/cut-name tissue,clumped_variant_id $(WORK)/clump-trans-variants/clump-trans-signals.significant-trans-eqtl-clumped.tsv | grep $(t) | cut -f2 | sort | uniq > $(ANALYSIS)/top-variants/$(t).variants.txt$(NL))
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --snp_mappability $(DATA)/crossmap-gencode-v30/snp_mappability_100mer_2mismatch.bed.gz --min_maf '0.05' --plink_glob '$(DATA)/genotypes/plink-pass-filter/*' --top_variants_glob '$(ANALYSIS)/top-variants/*.variants.txt' --sample_glob '$(ANALYSIS)/samples-wgs/*.samples.txt' $(ROOT)/get-maf-matched-snps-trans.nf &

# no hits for several cell types, so omitting them
control-snps-trans-sqtl: ANALYSIS=$(WORK)/control-snps/trans-sqtl/maf005
control-snps-trans-sqtl: TISSUES=Lung PBMC Whole_blood
control-snps-trans-sqtl:
	mkdir -p $(ANALYSIS)/samples-wgs
	mkdir -p $(ANALYSIS)/top-variants
	$(foreach t,$(TISSUES),$(BIN)/cut-name wgs,tor $(DATA)/metadata/metadata.tm.txt | grep -f $(DATA)/scan-samples/samples-to-use-for-scan.$(t).tsv | cut -f1 > $(ANALYSIS)/samples-wgs/$(t).samples.txt$(NL))
	$(foreach t,$(TISSUES),$(BIN)/cut-name tissue,clumped_variant_id $(WORK)/clump-trans-variants/clump-trans-signals.significant-trans-sqtl-clumped.tsv | grep $(t) | cut -f2 | sort | uniq > $(ANALYSIS)/top-variants/$(t).variants.txt$(NL))
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --snp_mappability $(DATA)/crossmap-gencode-v30/snp_mappability_100mer_2mismatch.bed.gz --min_maf '0.05' --plink_glob '$(DATA)/genotypes/plink-pass-filter/*'  --top_variants_glob '$(ANALYSIS)/top-variants/*.variants.txt' --sample_glob '$(ANALYSIS)/samples-wgs/*.samples.txt' $(ROOT)/get-maf-matched-snps-trans.nf &

variant-annotation-matrices: ANALYSIS=$(WORK)/variant-annotation-matrices
variant-annotation-matrices:
	mkdir -p $(ANALYSIS)
	cat $(WORK)/control-credible-sets/*/*/results/controls/*.control-credible-sets.tsv | grep -v "^phenotype" | cut -f2 | uniq | sort --parallel=10 | uniq > $(ANALYSIS)/variants.txt
	cat $(DATA)/scan-results/joint/cis-eqtl/susie/maf*/*.cs.txt | grep -v "^phenotype" | cut -f2 | uniq | sort --parallel=10 | uniq >> $(ANALYSIS)/variants.txt
	cat $(DATA)/scan-results/joint/cis-sqtl/susie/maf*/postprocessed/*.by-gene.cs.txt | grep -v "^phenotype" | cut -f2 | uniq | sort --parallel=10 | uniq >> $(ANALYSIS)/variants.txt
	cat $(DATA)/scan-results/joint/trans-*qtl/maf005/trans-susie/*.cs.txt | cut -f2 | grep -v -w variant_id >> $(ANALYSIS)/variants.txt
	cat $(WORK)/control-snps/trans-*qtl/maf005/results/controls/* | grep -v "^variant" | cut -f1 | uniq | sort --parallel=10 | uniq >> $(ANALYSIS)/variants.txt
	cat $(WORK)/control-snps/trans-*qtl/maf005/results/controls/* | grep -v "^variant" | cut -f2 | sort --parallel=10 | uniq >> $(ANALYSIS)/variants.txt
	cat $(ANALYSIS)/variants.txt | sort --parallel=10 | uniq > $(ANALYSIS)/variants.tmp && mv $(ANALYSIS)/variants.tmp $(ANALYSIS)/variants.txt
	cd $(ANALYSIS) && sbatch --mem-per-cpu=120G --partition=topmed-working --exclude=topmed,topmed[2-10] --time=1-00:00:00 --job-name=vam --wrap="singularity exec --bind /net docker://porchard/general:20241111 python $(BIN)/make-variant-annotation-matrices.py $(ANALYSIS)/variants.txt $(DATA)/sites $(DATA)/roadmap $(DATA)/ensembl-regulatory-build/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20221007.gff.gz annotations."

blood-as-proxy:
	cd $(ROOT)/notebooks && singularity exec --bind /net docker://porchard/general:20241111 jupyter nbconvert --to notebook --execute blood-as-proxy.ipynb --output blood-as-proxy.ipynb

summarize-cis-eqtl-coloc-with-cis-sqtl:
	cd $(ROOT)/notebooks && singularity exec --bind /net $(ROOT)/general.sif jupyter nbconvert --to notebook --execute $@.ipynb --output $@.ipynb

summarize-panukbb-coloc:
	cd $(ROOT)/notebooks && singularity exec --bind /net $(ROOT)/general.sif jupyter nbconvert --to notebook --execute $@.ipynb --output $@.ipynb

sample-ancestry-assignment:
	cd $(ROOT)/notebooks && singularity exec --bind /net $(ROOT)/general.sif jupyter nbconvert --to notebook --execute $@.ipynb --output $@.ipynb

rank-signals:
	cd $(ROOT)/notebooks && singularity exec --bind /net $(ROOT)/general.sif jupyter nbconvert --to notebook --execute $@.ipynb --output $@.ipynb

cis-rare-variants-summary:
	cd $(ROOT)/notebooks && singularity exec --bind /net docker://porchard/general_jl:20231026105312 jupyter nbconvert --to notebook --execute $@.ipynb --output $@.ipynb

functional-enrichments:
	cd $(ROOT)/notebooks && singularity exec --bind /net $(ROOT)/general.sif jupyter nbconvert --to notebook --execute calculate-functional-enrichments-set-level-logistic-regression.ipynb --output calculate-functional-enrichments-set-level-logistic-regression.ipynb

blood-cell-type-open-chromatin:
	cd $(ROOT)/notebooks && singularity exec --bind /net $(ROOT)/general.sif jupyter nbconvert --to notebook --execute $@.ipynb --output $@.ipynb

assign-signal-cell-type-specificity:
	cd $(ROOT)/notebooks && singularity exec --bind /net $(ROOT)/general.sif jupyter nbconvert --to notebook --execute $@.ipynb --output $@.ipynb

trans-go-enrichment:
	cd $(ROOT)/notebooks && jupyter nbconvert --to notebook --execute $@.ipynb --output $@.ipynb

trans-enrichment-in-cis:
	cd $(ROOT)/notebooks && jupyter nbconvert --to notebook --execute $@.ipynb --output $@.ipynb

replicate-trans-eqtlgen: ANALYSIS=$(WORK)/replicate-trans/eqtlgen
replicate-trans-eqtlgen:
	mkdir -p $(ANALYSIS)
	zcat /net/topmed11/working/porchard/eqtlgen-preprocessing/work/lift-and-tabix/trans-significant/results/tabixed/eqtlgen.txt.gz | cut -f5,7 | awk 'NR>1' > $(ANALYSIS)/pairs.txt
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --phenotypes $(DATA)/tensorqtl-in/joint/cis-eqtl/Whole_blood.tensorqtl-in.phenotypes.bed.gz --covariates $(DATA)/tensorqtl-in/joint/trans-eqtl/Whole_blood.tensorqtl-in.100.covariates.tsv --pairs $(ANALYSIS)/pairs.txt --vcf_dir $(DATA)/genotypes/vcfs-updated-ids-pass-filter $(ROOT)/replicate-trans.nf &

replicate-trans-direct: ANALYSIS=$(WORK)/replicate-trans/direct
replicate-trans-direct:
	mkdir -p $(ANALYSIS)
	zcat /net/topmed11/working/porchard/direct-preprocessing/work/lift-and-tabix/trans-significant/results/tabixed/direct.txt.gz | cut -f4,6 > $(ANALYSIS)/gene-variant.txt
	cut -f1 $(ANALYSIS)/gene-variant.txt > $(ANALYSIS)/gene.txt
	cut -f2 $(ANALYSIS)/gene-variant.txt > $(ANALYSIS)/variant.txt
	paste $(ANALYSIS)/variant.txt $(ANALYSIS)/gene.txt | grep -v "SNPid" | grep -v GL000 > $(ANALYSIS)/pairs.txt
	rm $(ANALYSIS)/gene.txt $(ANALYSIS)/gene-variant.txt $(ANALYSIS)/variant.txt
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --phenotypes $(DATA)/tensorqtl-in/joint/cis-eqtl/Whole_blood.tensorqtl-in.phenotypes.bed.gz --covariates $(DATA)/tensorqtl-in/joint/trans-eqtl/Whole_blood.tensorqtl-in.100.covariates.tsv --pairs $(ANALYSIS)/pairs.txt --vcf_dir $(DATA)/genotypes/vcfs-updated-ids-pass-filter $(ROOT)/replicate-trans.nf &

replicate-trans-gtex: replicate-trans-gtex-lung replicate-trans-gtex-whole-blood

replicate-trans-gtex-lung: ANALYSIS=$(WORK)/replicate-trans/gtex-lung
replicate-trans-gtex-lung:
	mkdir -p $(ANALYSIS)
	grep Lung /net/topmed11/working/porchard/gtex-preprocessing/data/gtex/GTEx_Analysis_v8_trans_eGenes_fdr05.txt | cut -f2,7 | perl -pe 's/_b38//' > $(ANALYSIS)/gene-variant.txt
	cut -f1 $(ANALYSIS)/gene-variant.txt > $(ANALYSIS)/gene.txt
	cut -f2 $(ANALYSIS)/gene-variant.txt > $(ANALYSIS)/variant.txt
	paste $(ANALYSIS)/variant.txt $(ANALYSIS)/gene.txt > $(ANALYSIS)/pairs.txt
	rm $(ANALYSIS)/gene.txt $(ANALYSIS)/gene-variant.txt $(ANALYSIS)/variant.txt
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --phenotypes $(DATA)/tensorqtl-in/joint/cis-eqtl/Lung.tensorqtl-in.phenotypes.bed.gz --covariates $(DATA)/tensorqtl-in/joint/trans-eqtl/Lung.tensorqtl-in.75.covariates.tsv --pairs $(ANALYSIS)/pairs.txt --vcf_dir $(DATA)/genotypes/vcfs-updated-ids-pass-filter $(ROOT)/replicate-trans.nf &

replicate-trans-gtex-whole-blood: ANALYSIS=$(WORK)/replicate-trans/gtex-whole-blood
replicate-trans-gtex-whole-blood:
	mkdir -p $(ANALYSIS)
	grep Whole_Blood /net/topmed11/working/porchard/gtex-preprocessing/data/gtex/GTEx_Analysis_v8_trans_eGenes_fdr05.txt | cut -f2,7 | perl -pe 's/_b38//' > $(ANALYSIS)/gene-variant.txt
	cut -f1 $(ANALYSIS)/gene-variant.txt > $(ANALYSIS)/gene.txt
	cut -f2 $(ANALYSIS)/gene-variant.txt > $(ANALYSIS)/variant.txt
	paste $(ANALYSIS)/variant.txt $(ANALYSIS)/gene.txt > $(ANALYSIS)/pairs.txt
	rm $(ANALYSIS)/gene.txt $(ANALYSIS)/gene-variant.txt $(ANALYSIS)/variant.txt
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --phenotypes $(DATA)/tensorqtl-in/joint/cis-eqtl/Whole_blood.tensorqtl-in.phenotypes.bed.gz --covariates $(DATA)/tensorqtl-in/joint/trans-eqtl/Whole_blood.tensorqtl-in.100.covariates.tsv --pairs $(ANALYSIS)/pairs.txt --vcf_dir $(DATA)/genotypes/vcfs-updated-ids-pass-filter $(ROOT)/replicate-trans.nf &

replication-of-gtex-trans-qtl:
	#cd $(ROOT)/notebooks && singularity exec --bind /net $(ROOT)/general.sif jupyter nbconvert --to notebook --execute $@.ipynb --output $@.ipynb
	cd $(ROOT)/notebooks && jupyter nbconvert --to notebook --execute $@.ipynb --output $@.ipynb

replication-of-eqtlgen-trans-eqtl:
	#cd $(ROOT)/notebooks && singularity exec --bind /net docker://porchard/general:20241111 jupyter nbconvert --to notebook --execute $@.ipynb --output $@.ipynb
	cd $(ROOT)/notebooks && jupyter nbconvert --to notebook --execute $@.ipynb --output $@.ipynb

replication-of-direct-trans-eqtl:
	#cd $(ROOT)/notebooks && singularity exec --bind /net docker://porchard/general_jl:20231026105312 jupyter nbconvert --to notebook --execute $@.ipynb --output $@.ipynb
	cd $(ROOT)/notebooks && jupyter nbconvert --to notebook --execute $@.ipynb --output $@.ipynb


tensorqtl-to-coloc-in:  tensorqtl-to-coloc-in-joint tensorqtl-to-coloc-in-ancestry-specific tensorqtl-to-coloc-in-saturation
tensorqtl-to-coloc-in-joint: tensorqtl-to-coloc-in-joint-cis-eqtl-maf001 tensorqtl-to-coloc-in-joint-cis-eqtl-maf0001 tensorqtl-to-coloc-in-joint-cis-sqtl-maf001 tensorqtl-to-coloc-in-joint-cis-sqtl-maf0001 tensorqtl-to-coloc-in-joint-trans-eqtl-maf005 tensorqtl-to-coloc-in-joint-trans-sqtl-maf005
tensorqtl-to-coloc-in-ancestry-specific: tensorqtl-to-coloc-in-ancestry-specific-cis-eqtl-maf001 tensorqtl-to-coloc-in-ancestry-specific-cis-sqtl-maf001 tensorqtl-to-coloc-in-ancestry-specific-trans-eqtl-maf005 tensorqtl-to-coloc-in-ancestry-specific-trans-sqtl-maf005
tensorqtl-to-coloc-in-saturation: tensorqtl-to-coloc-in-saturation-cis-eqtl-maf001

tensorqtl-to-coloc-in-joint-cis-eqtl-maf001: ANALYSIS=$(WORK)/tensorqtl-to-coloc-in/joint/cis-eqtl/maf001
tensorqtl-to-coloc-in-joint-cis-eqtl-maf001:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --bin $(BIN) --container 'docker.io/porchard/tensorqtl_dev:5ea048f_pyarrow_11' --susie_glob '$(DATA)/scan-results/joint/cis-eqtl/susie/maf001/*.pickle' --cs_glob '$(DATA)/scan-results/joint/cis-eqtl/susie/maf001/*.cs.txt' $(ROOT)/tensorqtl-susie-to-susieR.nf &

tensorqtl-to-coloc-in-joint-cis-eqtl-maf0001: ANALYSIS=$(WORK)/tensorqtl-to-coloc-in/joint/cis-eqtl/maf0001
tensorqtl-to-coloc-in-joint-cis-eqtl-maf0001:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --bin $(BIN) --container 'docker.io/porchard/tensorqtl_dev:5ea048f_pyarrow_11' --susie_glob '$(DATA)/scan-results/joint/cis-eqtl/susie/maf0001/*.pickle' --cs_glob '$(DATA)/scan-results/joint/cis-eqtl/susie/maf0001/*.cs.txt' $(ROOT)/tensorqtl-susie-to-susieR.nf &

tensorqtl-to-coloc-in-joint-cis-sqtl-maf001: ANALYSIS=$(WORK)/tensorqtl-to-coloc-in/joint/cis-sqtl/maf001
tensorqtl-to-coloc-in-joint-cis-sqtl-maf001:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --bin $(BIN) --container 'docker.io/porchard/tensorqtl_dev:5ea048f_pyarrow_11' --susie_glob '$(DATA)/scan-results/joint/cis-sqtl/susie/maf001/postprocessed/*.by-gene.pickle' --cs_glob '$(DATA)/scan-results/joint/cis-sqtl/susie/maf001/postprocessed/*.by-gene.cs.txt' $(ROOT)/tensorqtl-susie-to-susieR.nf &

tensorqtl-to-coloc-in-joint-cis-sqtl-maf0001: ANALYSIS=$(WORK)/tensorqtl-to-coloc-in/joint/cis-sqtl/maf0001
tensorqtl-to-coloc-in-joint-cis-sqtl-maf0001:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --bin $(BIN) --container 'docker.io/porchard/tensorqtl_dev:5ea048f_pyarrow_11' --susie_glob '$(DATA)/scan-results/joint/cis-sqtl/susie/maf0001/postprocessed/*.by-gene.pickle' --cs_glob '$(DATA)/scan-results/joint/cis-sqtl/susie/maf0001/postprocessed/*.by-gene.cs.txt' $(ROOT)/tensorqtl-susie-to-susieR.nf &

tensorqtl-to-coloc-in-joint-trans-eqtl-maf005: ANALYSIS=$(WORK)/tensorqtl-to-coloc-in/joint/trans-eqtl/maf005
tensorqtl-to-coloc-in-joint-trans-eqtl-maf005:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --bin $(BIN) --container 'docker.io/porchard/tensorqtl_dev:5ea048f_pyarrow_11' --susie_glob '$(DATA)/scan-results/joint/trans-eqtl/maf005/trans-susie/*.pickle' --cs_glob '$(DATA)/scan-results/joint/trans-eqtl/maf005/trans-susie/*.cs.txt' $(ROOT)/tensorqtl-susie-to-susieR.nf &

tensorqtl-to-coloc-in-joint-trans-sqtl-maf005: ANALYSIS=$(WORK)/tensorqtl-to-coloc-in/joint/trans-sqtl/maf005
tensorqtl-to-coloc-in-joint-trans-sqtl-maf005:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --bin $(BIN) --container 'docker.io/porchard/tensorqtl_dev:5ea048f_pyarrow_11' --susie_glob '$(DATA)/scan-results/joint/trans-sqtl/maf005/trans-susie/*.pickle' --cs_glob '$(DATA)/scan-results/joint/trans-sqtl/maf005/trans-susie/*.cs.txt' $(ROOT)/tensorqtl-susie-to-susieR.nf &

tensorqtl-to-coloc-in-ancestry-specific-cis-eqtl-maf001: ANALYSIS=$(WORK)/tensorqtl-to-coloc-in/ancestry-specific/cis-eqtl/maf001
tensorqtl-to-coloc-in-ancestry-specific-cis-eqtl-maf001:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --bin $(BIN) --container 'docker.io/porchard/tensorqtl_dev:5ea048f_pyarrow_11' --susie_glob '$(DATA)/scan-results/ancestry-specific/cis-eqtl/susie/maf001/*.pickle' --cs_glob '$(DATA)/scan-results/ancestry-specific/cis-eqtl/susie/maf001/*.cs.txt' $(ROOT)/tensorqtl-susie-to-susieR.nf &

tensorqtl-to-coloc-in-ancestry-specific-cis-sqtl-maf001: ANALYSIS=$(WORK)/tensorqtl-to-coloc-in/ancestry-specific/cis-sqtl/maf001
tensorqtl-to-coloc-in-ancestry-specific-cis-sqtl-maf001:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --bin $(BIN) --container 'docker.io/porchard/tensorqtl_dev:5ea048f_pyarrow_11' --susie_glob '$(DATA)/scan-results/ancestry-specific/cis-sqtl/susie/maf001/postprocessed/*.by-gene.pickle' --cs_glob '$(DATA)/scan-results/ancestry-specific/cis-sqtl/susie/maf001/postprocessed/*.by-gene.cs.txt' $(ROOT)/tensorqtl-susie-to-susieR.nf &

tensorqtl-to-coloc-in-ancestry-specific-trans-eqtl-maf005: ANALYSIS=$(WORK)/tensorqtl-to-coloc-in/ancestry-specific/trans-eqtl/maf005
tensorqtl-to-coloc-in-ancestry-specific-trans-eqtl-maf005:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --bin $(BIN) --container 'docker.io/porchard/tensorqtl_dev:5ea048f_pyarrow_11' --susie_glob '$(DATA)/scan-results/ancestry-specific/trans-eqtl/maf005/trans-susie/*.pickle' --cs_glob '$(DATA)/scan-results/ancestry-specific/trans-eqtl/maf005/trans-susie/*.cs.txt' $(ROOT)/tensorqtl-susie-to-susieR.nf &

tensorqtl-to-coloc-in-ancestry-specific-trans-sqtl-maf005: ANALYSIS=$(WORK)/tensorqtl-to-coloc-in/ancestry-specific/trans-sqtl/maf005
tensorqtl-to-coloc-in-ancestry-specific-trans-sqtl-maf005:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --bin $(BIN) --container 'docker.io/porchard/tensorqtl_dev:5ea048f_pyarrow_11' --susie_glob '$(DATA)/scan-results/ancestry-specific/trans-sqtl/maf005/trans-susie/*.pickle' --cs_glob '$(DATA)/scan-results/ancestry-specific/trans-sqtl/maf005/trans-susie/*.cs.txt' $(ROOT)/tensorqtl-susie-to-susieR.nf &

tensorqtl-to-coloc-in-saturation-cis-eqtl-maf001: ANALYSIS=$(WORK)/tensorqtl-to-coloc-in/saturation/cis-eqtl/maf001
tensorqtl-to-coloc-in-saturation-cis-eqtl-maf001:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --bin $(BIN) --container 'docker.io/porchard/tensorqtl_dev:5ea048f_pyarrow_11' --susie_glob '$(DATA)/scan-results/saturation/cis-eqtl/susie/maf001/*.100.susie.pickle' --cs_glob '$(DATA)/scan-results/saturation/cis-eqtl/susie/maf001/*.100.cs.txt' $(ROOT)/tensorqtl-susie-to-susieR.nf &

susie-json-joint: ANALYSIS=$(WORK)/coloc/json
susie-json-joint:
	mkdir -p $(ANALYSIS)
	singularity exec --bind /net docker://porchard/general:20241111 python $(BIN)/make-susie-json-joint.py > $(ANALYSIS)/joint.json

susie-json-ancestry-specific: ANALYSIS=$(WORK)/coloc/json
susie-json-ancestry-specific:
	mkdir -p $(ANALYSIS)
	singularity exec --bind /net docker://porchard/general:20241111 python $(BIN)/make-susie-json-ancestry-specific.py > $(ANALYSIS)/ancestry-specific.json

susie-json-saturation: ANALYSIS=$(WORK)/coloc/json
susie-json-saturation:
	mkdir -p $(ANALYSIS)
	singularity exec --bind /net docker://porchard/general:20241111 python $(BIN)/make-susie-json-saturation.py > $(ANALYSIS)/saturation.json

susie-json-panukbb: ANALYSIS=$(WORK)/coloc/json
susie-json-panukbb:
	mkdir -p $(ANALYSIS)
	singularity exec --bind /net docker://porchard/general:20241111 python $(BIN)/make-susie-json-panukbb.py > $(ANALYSIS)/panukbb.json

coloc-stage-joint: INPUT_JSON=$(WORK)/coloc/json/joint.json
coloc-stage-joint: ANALYSIS=/tmp/rnaseq-2024-10-16-sample-update/stage/coloc/joint
coloc-stage-joint: OUTPUT_JSON=$(ANALYSIS)/joint.json
coloc-stage-joint:
	ssh topmed11 'mkdir -p $(ANALYSIS) && cd $(ANALYSIS) && sbatch --time=10:00:00 --mail-user=porchard@umich.edu --mail-type=END,FAIL --partition=topmed-working --exclude=topmed,topmed[2-10] --wrap="python $(BIN)/stage.py --stage-using modality ancestry tissue maf --json $(INPUT_JSON) > $(OUTPUT_JSON)"'

coloc-stage-ancestry-specific: INPUT_JSON=$(WORK)/coloc/json/ancestry-specific.json
coloc-stage-ancestry-specific: ANALYSIS=/tmp/rnaseq-2024-10-16-sample-update/stage/coloc/ancestry-specific
coloc-stage-ancestry-specific: OUTPUT_JSON=$(ANALYSIS)/ancestry-specific.json
coloc-stage-ancestry-specific:
	ssh topmed11 'mkdir -p $(ANALYSIS) && cd $(ANALYSIS) && sbatch --time=10:00:00 --mail-user=porchard@umich.edu --mail-type=END,FAIL --partition=topmed-working --exclude=topmed,topmed[2-10] --wrap="python $(BIN)/stage.py --stage-using modality ancestry tissue maf --json $(INPUT_JSON) > $(OUTPUT_JSON)"'

coloc-stage-saturation: INPUT_JSON=$(WORK)/coloc/json/saturation.json
coloc-stage-saturation: ANALYSIS=/tmp/rnaseq-2024-10-16-sample-update/stage/coloc/saturation
coloc-stage-saturation: OUTPUT_JSON=$(ANALYSIS)/saturation.json
coloc-stage-saturation:
	ssh topmed11 'mkdir -p $(ANALYSIS) && cd $(ANALYSIS) && sbatch --time=10:00:00 --mail-user=porchard@umich.edu --mail-type=END,FAIL --partition=topmed-working --exclude=topmed,topmed[2-10] --wrap="python $(BIN)/stage.py --stage-using modality ancestry tissue maf --json $(INPUT_JSON) > $(OUTPUT_JSON)"'

coloc-stage-panukbb: INPUT_JSON=$(WORK)/coloc/json/panukbb.json
coloc-stage-panukbb: ANALYSIS=/tmp/rnaseq-2024-10-16-sample-update/stage/coloc/panukbb
coloc-stage-panukbb: OUTPUT_JSON=$(ANALYSIS)/panukbb.json
coloc-stage-panukbb:
	ssh topmed11 'mkdir -p $(ANALYSIS) && cd $(ANALYSIS) && sbatch --time=10:00:00 --mail-user=porchard@umich.edu --mail-type=END,FAIL --partition=topmed-working --exclude=topmed,topmed[2-10] --wrap="python $(BIN)/stage.py --stage-using modality ancestry tissue maf --json $(INPUT_JSON) > $(OUTPUT_JSON)"'

coloc-panukbb-joint: ANALYSIS=$(WORK)/coloc/panukbb/joint
coloc-panukbb-joint: SIF=singularity exec --bind /net docker://porchard/general_jl:20231026105312
coloc-panukbb-joint:
	mkdir -p $(ANALYSIS)
	ssh topmed11 '$(SIF) python $(BIN)/filter.py --modality ciseqtl cissqtl transeqtl transsqtl --ancestry joint --input /tmp/rnaseq-2024-10-16-sample-update/stage/coloc/joint/joint.json > $(ANALYSIS)/xqtl.json'
	ssh topmed11 '$(SIF) python $(BIN)/filter.py --modality gwas --input /tmp/rnaseq-2024-10-16-sample-update/stage/coloc/panukbb/panukbb.json > $(ANALYSIS)/gwas.json'
	$(SIF) python $(BIN)/get-coloc-pairs.py $(ANALYSIS)/xqtl.json $(ANALYSIS)/gwas.json > $(ANALYSIS)/pairs.txt
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --label_1 xqtl --label_2 gwas --pairs $(ANALYSIS)/pairs.txt $(ROOT)/run-coloc-batch-preload.nf &

postprocess-coloc-panukbb-joint: ANALYSIS=$(WORK)/coloc/panukbb/joint
postprocess-coloc-panukbb-joint:
	$(SIF) python $(BIN)/postprocess-panukbb-coloc.py --colocs-only --panukbb-root /net/topmed11/working/porchard/panukbb-finemapping --coloc-out $(ANALYSIS)/results/txt/coloc.txt --susie-json $(ANALYSIS)/*.json > $(ANALYSIS)/postprocessed.txt
	$(SIF) python $(BIN)/postprocess-panukbb-coloc.py --panukbb-root /net/topmed11/working/porchard/panukbb-finemapping --coloc-out $(ANALYSIS)/results/txt/coloc.txt --susie-json $(ANALYSIS)/*.json > $(ANALYSIS)/all-tested-pairs.txt

coloc-panukbb-ancestry-specific: ANALYSIS=$(WORK)/coloc/panukbb/ancestry-specific
coloc-panukbb-ancestry-specific: SIF=singularity exec --bind /net docker://porchard/general_jl:20231026105312
coloc-panukbb-ancestry-specific:
	mkdir -p $(ANALYSIS)
	ssh topmed11 '$(SIF) python $(BIN)/filter.py --modality ciseqtl cissqtl transeqtl transsqtl --ancestry AFR --input /tmp/rnaseq-2024-10-16-sample-update/stage/coloc/ancestry-specific/ancestry-specific.json > $(ANALYSIS)/xqtl.AFR.json'
	ssh topmed11 '$(SIF) python $(BIN)/filter.py --modality ciseqtl cissqtl transeqtl transsqtl --ancestry EUR --input /tmp/rnaseq-2024-10-16-sample-update/stage/coloc/ancestry-specific/ancestry-specific.json > $(ANALYSIS)/xqtl.EUR.json'
	ssh topmed11 '$(SIF) python $(BIN)/filter.py --modality gwas --ancestry AFR --input /tmp/rnaseq-2024-10-16-sample-update/stage/coloc/panukbb/panukbb.json > $(ANALYSIS)/gwas.AFR.json'
	ssh topmed11 '$(SIF) python $(BIN)/filter.py --modality gwas --ancestry EUR --input /tmp/rnaseq-2024-10-16-sample-update/stage/coloc/panukbb/panukbb.json > $(ANALYSIS)/gwas.EUR.json'
	$(SIF) python $(BIN)/get-coloc-pairs.py $(ANALYSIS)/xqtl.AFR.json $(ANALYSIS)/gwas.AFR.json > $(ANALYSIS)/pairs.txt
	$(SIF) python $(BIN)/get-coloc-pairs.py $(ANALYSIS)/xqtl.EUR.json $(ANALYSIS)/gwas.EUR.json >> $(ANALYSIS)/pairs.txt
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --label_1 xqtl --label_2 gwas --pairs $(ANALYSIS)/pairs.txt $(ROOT)/run-coloc-batch-preload.nf &

postprocess-coloc-panukbb-ancestry-specific: ANALYSIS=$(WORK)/coloc/panukbb/ancestry-specific
postprocess-coloc-panukbb-ancestry-specific:
	$(SIF) python $(BIN)/postprocess-panukbb-coloc.py --colocs-only --panukbb-root /net/topmed11/working/porchard/panukbb-finemapping --coloc-out $(ANALYSIS)/results/txt/coloc.txt --susie-json $(ANALYSIS)/*.json > $(ANALYSIS)/postprocessed.txt
	$(SIF) python $(BIN)/postprocess-panukbb-coloc.py --panukbb-root /net/topmed11/working/porchard/panukbb-finemapping --coloc-out $(ANALYSIS)/results/txt/coloc.txt --susie-json $(ANALYSIS)/*.json > $(ANALYSIS)/all-tested-pairs.txt

coloc-panukbb-saturation: ANALYSIS=$(WORK)/coloc/panukbb/saturation
coloc-panukbb-saturation: SIF=singularity exec --bind /net docker://porchard/general_jl:20231026105312
coloc-panukbb-saturation:
	mkdir -p $(ANALYSIS)
	ssh topmed11 '$(SIF) python $(BIN)/filter.py --modality ciseqtl cissqtl transeqtl transsqtl --ancestry joint --input /tmp/rnaseq-2024-10-16-sample-update/stage/coloc/saturation/saturation.json > $(ANALYSIS)/xqtl.json'
	ssh topmed11 '$(SIF) python $(BIN)/filter.py --modality gwas --input /tmp/rnaseq-2024-10-16-sample-update/stage/coloc/panukbb/panukbb.json > $(ANALYSIS)/gwas.json'
	$(SIF) python $(BIN)/get-coloc-pairs.py $(ANALYSIS)/xqtl.json $(ANALYSIS)/gwas.json > $(ANALYSIS)/pairs.txt
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --label_1 xqtl --label_2 gwas --pairs $(ANALYSIS)/pairs.txt $(ROOT)/run-coloc-batch-preload.nf &

postprocess-coloc-panukbb-saturation: ANALYSIS=$(WORK)/coloc/panukbb/saturation
postprocess-coloc-panukbb-saturation:
	$(SIF) python $(BIN)/postprocess-panukbb-coloc.py --colocs-only --panukbb-root /net/topmed11/working/porchard/panukbb-finemapping --coloc-out $(ANALYSIS)/results/txt/coloc.txt --susie-json $(ANALYSIS)/*.json > $(ANALYSIS)/postprocessed.txt
	$(SIF) python $(BIN)/postprocess-panukbb-coloc.py --panukbb-root /net/topmed11/working/porchard/panukbb-finemapping --coloc-out $(ANALYSIS)/results/txt/coloc.txt --susie-json $(ANALYSIS)/*.json > $(ANALYSIS)/all-tested-pairs.txt

coloc-xqtl-joint: TISSUES=Lung Nasal_epithelial T_cell PBMC Monocyte Whole_blood
coloc-xqtl-joint: TRANS_SQTL_TISSUES=Lung PBMC Whole_blood
coloc-xqtl-joint: JSON=/tmp/rnaseq-2024-10-16-sample-update/stage/coloc/joint/joint.json
coloc-xqtl-joint: ANALYSIS=$(WORK)/coloc/xqtl/joint
coloc-xqtl-joint:
	mkdir -p $(ANALYSIS)
	rm -rf $(ANALYSIS)/pairs.txt
	$(foreach t,$(TISSUES),ssh topmed11 '$(SIF) python $(BIN)/filter.py --tissue $(t) --modality ciseqtl --maf 1% --ancestry joint --input $(JSON) > $(ANALYSIS)/$(t).ciseqtl.json'$(NL))
	$(foreach t,$(TISSUES),ssh topmed11 '$(SIF) python $(BIN)/filter.py --tissue $(t) --modality cissqtl --maf 1% --ancestry joint --input $(JSON) > $(ANALYSIS)/$(t).cissqtl.json'$(NL))
	$(foreach t,$(TISSUES),ssh topmed11 '$(SIF) python $(BIN)/filter.py --tissue $(t) --modality transeqtl --maf 5% --ancestry joint --input $(JSON) > $(ANALYSIS)/$(t).transeqtl.json'$(NL))
	$(foreach t,$(TRANS_SQTL_TISSUES),ssh topmed11 '$(SIF) python $(BIN)/filter.py --tissue $(t) --modality transsqtl --maf 5% --ancestry joint --input $(JSON) > $(ANALYSIS)/$(t).transsqtl.json'$(NL))
	$(foreach t,$(TISSUES),$(SIF) python $(BIN)/get-coloc-pairs.py $(ANALYSIS)/$(t).ciseqtl.json $(ANALYSIS)/$(t).cissqtl.json >> $(ANALYSIS)/pairs.txt$(NL))
	$(foreach t,$(TISSUES),$(SIF) python $(BIN)/get-coloc-pairs.py $(ANALYSIS)/$(t).ciseqtl.json $(ANALYSIS)/$(t).transeqtl.json >> $(ANALYSIS)/pairs.txt$(NL))
	$(foreach t,$(TRANS_SQTL_TISSUES),$(SIF) python $(BIN)/get-coloc-pairs.py $(ANALYSIS)/$(t).ciseqtl.json $(ANALYSIS)/$(t).transsqtl.json >> $(ANALYSIS)/pairs.txt$(NL))
	$(foreach t,$(TISSUES),$(SIF) python $(BIN)/get-coloc-pairs.py $(ANALYSIS)/$(t).cissqtl.json $(ANALYSIS)/$(t).transeqtl.json >> $(ANALYSIS)/pairs.txt$(NL))
	$(foreach t,$(TRANS_SQTL_TISSUES),$(SIF) python $(BIN)/get-coloc-pairs.py $(ANALYSIS)/$(t).cissqtl.json $(ANALYSIS)/$(t).transsqtl.json >> $(ANALYSIS)/pairs.txt$(NL))
	$(foreach t,$(TRANS_SQTL_TISSUES),$(SIF) python $(BIN)/get-coloc-pairs.py $(ANALYSIS)/$(t).transeqtl.json $(ANALYSIS)/$(t).transsqtl.json >> $(ANALYSIS)/pairs.txt$(NL))
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --label_1 xqtl1 --label_2 xqtl2 --pairs $(ANALYSIS)/pairs.txt $(ROOT)/run-coloc-batch-preload.nf &

postprocess-coloc-xqtl-joint: ANALYSIS=$(WORK)/coloc/xqtl/joint
postprocess-coloc-xqtl-joint:
	$(SIF) python $(BIN)/postprocess-xqtl-coloc.py --colocs-only --coloc-out $(ANALYSIS)/results/txt/coloc.txt --susie-json $(ANALYSIS)/*.json > $(ANALYSIS)/postprocessed.txt
	$(SIF) python $(BIN)/postprocess-xqtl-coloc.py --coloc-out $(ANALYSIS)/results/txt/coloc.txt --susie-json $(ANALYSIS)/*.json > $(ANALYSIS)/all-tested-pairs.txt


eqtlgen-and-direct-ld-buddies: ANALYSIS=$(WORK)/eqtlgen-and-direct-ld-buddies
eqtlgen-and-direct-ld-buddies: TISSUES=Whole_blood
eqtlgen-and-direct-ld-buddies:
	#mkdir -p $(ANALYSIS)/samples-wgs
	#cut -f3 /net/topmed11/working/porchard/eqtlgen-preprocessing/work/top-hit-per-gene-after-lifting/top-per-gene.txt | awk 'NR>1' | sort | uniq > $(ANALYSIS)/eqtlgen-variants.txt
	#zcat /net/topmed11/working/porchard/direct-preprocessing/work/lift-and-tabix/cis-eqtl-significant/results/tabixed/direct.txt.gz | cut -f7 | awk 'NR>1' | sort | uniq > $(ANALYSIS)/direct-variants.txt
	#cat $(ANALYSIS)/eqtlgen-variants.txt $(ANALYSIS)/direct-variants.txt | sort | uniq > $(ANALYSIS)/variants.txt
	$(foreach t,$(TISSUES),# $(BIN)/cut-name wgs,tor $(DATA)/metadata/metadata.tm.txt | grep -f $(DATA)/scan-samples/samples-to-use-for-scan.$(t).tsv | cut -f1 > $(ANALYSIS)/samples-wgs/$(t).samples.txt$(NL))
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --variants $(ANALYSIS)/variants.txt --plink_glob '$(DATA)/genotypes/plink-pass-filter/*' --sample_glob '$(ANALYSIS)/samples-wgs/*.samples.txt' $(ROOT)/get-ld-buddies.nf &
