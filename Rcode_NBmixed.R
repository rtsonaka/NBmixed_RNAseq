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
#       - formula.lmer = a two-sided linear formula object to be used in lmer(.) describing both the 
#                        fixed-effects and random-effects part of the model, 
#                        with the response on the left of a ~ operator and the terms, 
#                        separated by + operators, on the right. 
#                        Regarding the random-effects terms only random-intercepts are used at the moment.
#       - data = dataframe containing the variables named in formula and the RNAseq data in each column per genomic feature.
#       - Time = a character string with the variable name in data for the followup time.
#       - ID = a character string with the variable name in data for the subject indicator.
#       - gene.nams = a character string with the gene names in data.
#       - dmat.voom = design matrix to be used in limma.
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
                        data = data.,
                        nboots = 2,
                        n.genes = 5,
                        gene.nams = as.character(1:5),
                        test.coefs = c("group", "time.:group"),
                        nAGQ = 11)

# Results
diff.expr$plrt.boot # pvalues of LRT per gene using the bootstrap method
diff.expr$plrt # pvalues of LRT per gene (using the asymptotic chisq_2)
diff.expr$var.b # estimates of random effects variance per gene
diff.expr$disp # estimates of dispersion per gene