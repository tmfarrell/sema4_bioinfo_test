## 
## vcf_to_mutation_types_table.R
## 
## Tim Farrell (tfarrell01@gmail.com) 
## 20190319
## 

# useful functions 
import <- function(package) {
    if (!is.element(package, installed.packages()[,1]))
        install.packages(package, dep=TRUE, repos="http://cran.us.r-project.org")
    require(package, character.only=TRUE)
}

# import libraries
import("vcfR")
import("dplyr")
import("ggplot2")
import("argparse")

# build cmd-line parser 
parser = ArgumentParser(description='Print VCF "mutation types", in tabular format.')
parser$add_argument('--vcf', required=TRUE, help='Path to vcf.')
parser$add_argument('--mutation_types_table', required=TRUE, help='Path to tabular file of mutation types.')

# parse cmd-line args 
args = parser$parse_args(commandArgs(trailingOnly=TRUE))

# read mutation types file 
mut_types = read.table(args$mutation_types_table, header=TRUE)

# read VCF and get each variant's ref/alt alleles 
vcf = read.vcfR(args$vcf, verbose=FALSE)
ref_alt_alleles = as.data.frame(getFIX(vcf)[,c("REF","ALT")])
colnames(ref_alt_alleles) = c("ref","alt")

# add the mutation type column to ref/ alt alleles 
ref_alt_alleles = left_join(ref_alt_alleles, mut_types, by=c("ref","alt"))

# count # of each mutation type 
mut_type_counts = as.data.frame(table(ref_alt_alleles["mutation_type"]))
colnames(mut_type_counts) = c("mutation_type", "count")

# plot and save 
p = ggplot(data=mut_type_counts, aes(x=mutation_type, y=count)) + geom_bar(stat="identity")
ggsave(paste0(strsplit(basename(args$vcf), "\\.")[[1]][1], ".mutation_type_counts.png"), plot=p)

