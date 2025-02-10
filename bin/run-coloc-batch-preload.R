#!/usr/bin/env Rscript
# coding: utf-8

# library(optparse)

# option_list <- list(
#   make_option(c("--input-list"), action = "store", type = "character", default = "", help = "[Required] File w/ two columns of RDA paths"),
#   make_option(c("--label-1"), action = "store", type = "character", default = "gwas", help = "[Optional] Label for modality represented by first column of --input-list (default: gwas)"),
#   make_option(c("--label-2"), action = "store", type = "character", default = "eqtl", help = "[Optional] Label for modality represented by second column of --input-list (default: eqtl)"),
#   make_option(c("--output"), action = "store", type = "character", default = "", help = "[Required] File to write"),
#   make_option(c("--fix-regionwide-priors"), action = "store_true", type = "logical", default = F, help = "[Optional] Select per-SNP priors as to keep a fixed regionwide prior")
# )

# option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
# opts <- parse_args(option_parser)

library(argparser)

# Create a parser
p <- arg_parser("Run coloc on a list of RDA file pairs")

# Add command line arguments
p <- add_argument(p, "--label-1", default='gwas', help="[Optional] Label for modality represented by first column of --input-list (default: gwas)")
p <- add_argument(p, "--label-2", default='eqtl', help="[Optional] Label for modality represented by second column of --input-list (default: eqtl)")
p <- add_argument(p, "--fix-regionwide-priors", type='logical', flag=T, help=" Select per-SNP priors as to keep a fixed regionwide prior")
p <- add_argument(p, "input_list", help="[Required] File w/ two columns of RDA paths")
p <- add_argument(p, "output", help="[Required] File to write")

# Parse the command line arguments
opts <- parse_args(p)


library(coloc)
library(glue)
library(dplyr)

args <- commandArgs(T)
INPUT_FILES <- opts$input_list
FIX_REGIONWIDE_PRIORS <- opts$fix_regionwide_priors
OUT <- opts$output
LABEL_1 <- opts$label_1
LABEL_2 <- opts$label_2

calculate_regionwide_priors <- function(p1=1e-4, p2=1e-4, p12=5e-6, nsnp=1000) {
    tmp <- c(nsnp * p1,
                nsnp * p2,
                nsnp * (nsnp-1) * p1 * p2,
                nsnp * p12)
    tmp <- c(1-sum(tmp), tmp)
    names(tmp) <- paste0('H', seq(0, 4))
    return(tmp)
}

calculate_priors <- function(H0, H1, H2, H3, H4, nsnps) {
    p1 <- H1/nsnps
    p2 <- H2/nsnps
    p12 <- H4/nsnps
    tmp <- c(p1, p2, p12)
    names(tmp) <- c('p1', 'p2', 'p12')
    return(tmp)
}

defaults <- calculate_regionwide_priors()

coloc_df <- read.table(INPUT_FILES, sep='\t', head=F, as.is=T, col.names=c('rda_1', 'rda_2'))


# load up all the RDA files
message(Sys.time())
message('Loading R objects')
load_rda <- function(rda)
{
    rda_name <- load(rda)
    rda_obj <- get(rda_name)
    return(rda_obj)
}


rda_files <- unique(c(coloc_df$rda_1, coloc_df$rda_2))
rda_list <- lapply(rda_files, load_rda)
names(rda_list) <- rda_files

rda_from_tensorqtl <- lapply(rda_list, function(x){!'susie' %in% class(x)})
message(Sys.time())
message('Finished loading R objects')




run_coloc <- function(rda_1, rda_2, fix_regionwide_priors, rda_1_label, rda_2_label, rda_1_from_tensorqtl, rda_2_from_tensorqtl, rda_1_file, rda_2_file) {

    # make coloc.susie represent these as complete susie objects
    class(rda_1) <- c(class(rda_1), 'susie')
    class(rda_2) <- c(class(rda_2), 'susie')

    # strip _b38 from variant names (for GTEx)
    colnames(rda_1$lbf_variable) <- gsub('_b38', '', colnames(rda_1$lbf_variable))
    colnames(rda_2$lbf_variable) <- gsub('_b38', '', colnames(rda_2$lbf_variable))

    if (fix_regionwide_priors) {
        number_intersecting_snps <- length(intersect(colnames(rda_1$lbf_variable), colnames(rda_2$lbf_variable)))
        new_priors <- calculate_priors(defaults['H0'], defaults['H1'], defaults['H2'], defaults['H3'], defaults['H4'], number_intersecting_snps)
        names(new_priors) <- NULL
        coloc_out <- coloc.susie(rda_1, rda_2, p1=new_priors[1], p2=new_priors[2], p12=new_priors[3])
    } else {
        coloc_out <- coloc.susie(rda_1, rda_2)
    }

    # if (is.data.frame(coloc_out$summary)) {
    #     #sens <- lapply(1:nrow(coloc_out$summary), function(x){sensitivity(coloc_out, rule="H4 > 0.5", doplot=F, plot.manhattans=F, row=x)})
    #     sens <- NULL
    # } else {
    #     sens <- NULL
    # }

    # add some additional metadata about SNPs present
    #coloc_out$nsnps <- c(ncol(rda_1$lbf_variable), ncol(rda_2$lbf_variable))
    #names(coloc_out$nsnps) <- c(rda_1_label, rda_2_label)
    #coloc_out[[glue('{rda_1_label}_snps')]] <- colnames(rda_1$lbf_variable)
    #coloc_out[[glue('{rda_2_label}_snps')]] <- colnames(rda_2$lbf_variable)
    # coloc_out$sensitivity <- sens

    snp_positions <- as.numeric(sapply(strsplit(coloc_out$results$snp, '_'), function(x){x[2]}))
    region_start <- min(snp_positions) - 1
    region_end <- max(snp_positions)
    s <- coloc_out$summary
    s[[glue('{rda_1_label}_nsnps')]] <- coloc_out$nsnps[rda_1_label]
    s[[glue('{rda_2_label}_nsnps')]] <- coloc_out$nsnps[rda_2_label]
    s$chrom <- strsplit(coloc_out$results$snp[1], '_')[[1]][1]
    s$region_start <- region_start
    s$region_end <- region_end 
    s[[glue('{rda_1_label}_cs')]] <- if (!rda_1_from_tensorqtl) {paste0('L', s$idx1)} else {names(rda_1$sets$cs)[s$idx1]}
    s[[glue('{rda_2_label}_cs')]] <- if (!rda_2_from_tensorqtl) {paste0('L', s$idx2)} else {names(rda_2$sets$cs)[s$idx2]}
    s$p1 <- coloc_out$priors['p1']
    s$p2 <- coloc_out$priors['p2']
    s$p12 <- coloc_out$priors['p12']
    s[[glue('{rda_1_label}_file')]] <- rda_1_file
    s[[glue('{rda_2_label}_file')]] <- rda_2_file
    return(s)
}

results <- list()
lines_to_process <- nrow(coloc_df)

for(i in 1:nrow(coloc_df)) {
    if(i %% 1000 == 0) {
        message(Sys.time())
        message(glue('Processed {i} of {lines_to_process} lines'))
    }
    tryCatch(
        {
            tmp <- run_coloc(rda_list[[coloc_df$rda_1[i]]], rda_list[[coloc_df$rda_2[i]]], FIX_REGIONWIDE_PRIORS, LABEL_1, LABEL_2, rda_from_tensorqtl[[coloc_df$rda_1[i]]], rda_from_tensorqtl[[coloc_df$rda_2[i]]], coloc_df$rda_1[i], coloc_df$rda_2[i])
            results[[length(results)+1]] <- tmp
        }, error=function(cond) {
            M <- paste('Failed', coloc_df$rda_1[i], coloc_df$rda_2[i], sep=' - ')
            #message(M)
            #message(cond)
        },
        warning=function(cond) {
            M <- paste('Failed', coloc_df$rda_1[i], coloc_df$rda_2[i], sep=' - ')
            #message(M)
            #message(cond)
        }
    )
}

results <- bind_rows(results)
write.table(results, OUT, sep='\t', quote=F, row.names=F)