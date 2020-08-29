# NBmixed_RNAseq
Modelling longitudinal RNAseq data using Negative Binomial mixed models.


* Simulate_Data.R: contains code  to simulate longitudinal RNAseq gene expression data as described in Tsonaka \& Spitali (2020). 
* MainFunction.R: Main function to model longitudinal RNAseq data using Negative Binomial mixed models via the GLMMadaptive R package and the bootstrap method to empirically derive the p-values per gene.  The output contains a list with elements: (1) gene names, (2) estimates of random effects variance per gene (3) estimates of dispersion per gene and (4) p-values of the LRT per gene using the asymptotic chisq distribution and (5) p-values of the LRT per gene using the bootstrap method. 
* Rcode_NBmixed.R: Example to model longitudinal RNAseq data using Negative Binomial mixed models and bootstrap sampling.
* Dataset.txt: Contains a toy dataset with longitudinal RNAseq collected on 10 subjects at 5 occasions: baseline and at 0.5, 1.0, 1.5, 2.0 weeks, afterwards. Subjects are assumed to be assigned to 2 treatments groups A and B and our interest is to test group differences over time. The data have been simulated from a Negative Binomial mixed effects model and the first half of the genes are assumed to be differentially expressed between the 2 groups at baseline and over time.

This code has been used to run the simulation study descibed in Tsonaka 
& Spitali (2020). 

Reference: Tsonaka, R. and Spitali, P. (2020). Negative Binomial mixed models under the maximum likelihood method can be used for longitudinal RNAseq data. Under review.

###
