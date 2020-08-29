#########################################################################################################
# Author: Roula Tsonaka (s.tsonaka@lumc.nl)
#         Leiden University Medical Center (LUMC)
#########################################################################################################
# Title: Rcode_NBmixed.R
# Aim: R code to analyse longitudinal RNAseq data using 
#      Negative Binomial mixed models. 
# Notes: 
#      This code can be used to analyse longitudinal RNAseq data 
#      via the GLMMadaptive R package and the boostrap method to 
#      empirically derive the p-values per gene.
#      The main function is longRNAseq(.) with arguments:
#       - formula.nbmixed = a two-sided linear formula object to be used in mixed_models(.) 
#                        describing the fixed-effects part of the model under the alternative 
#                        hypothesis, with the response on the left of a ~ operator and the terms, 
#                        separated by + operators, on the right. 
#                        Regarding the random-effects terms only random-intercepts are used at the moment.
#       - formula.nbmixed.0 = a two-sided linear formula object to be used in mixed_models(.) 
#                        describing the fixed-effects part of the model under the null 
#                        hypothesis.
#       - formula.random = a two-sided linear formula object to be used in mixed_models(.) 
#                        describing the random-effects part of the model.
#       - data = dataframe containing the variables named in formula and the RNAseq data in each column per genomic feature.
#       - nboots = numeric scalar denoting the number of bootstrap samples.
#       - test.coefs = a character string with the the names of the fixed-effects parameters to be tested.
#       - n.genes = numeric scalar denoting the number of genes.
#       - gene.nams = a character string with the gene names in data.
#       - nAGQ = numeric scalar denoting the number of quadrature points.
#
# Date: 27AUGUST2020
###########################################################################################################


# Load main function
source("./MainFunction.R")

# Load data
data. <- read.table("./Dataset.txt")


#####################

# Differential expression testing
diff.expr <- longRNAseq(formula.nbmixed = "time. + group + time.:group",
                        formula.nbmixed.0 = "time.",
                        formula.random = "~1|ID", 
                        data = data.,
                        nboots = 2,
                        n.genes = 5,
                        gene.nams = paste("X", 1:5, sep = ""),
                        test.coefs = c("group", "time.:group"),
                        nAGQ = 11)

# Results
diff.expr$plrt.boot # pvalues of LRT per gene using the bootstrap method
diff.expr$plrt # pvalues of LRT per gene (using the asymptotic chisq_2)
diff.expr$var.b # estimates of random effects variance per gene
diff.expr$disp # estimates of dispersion per gene