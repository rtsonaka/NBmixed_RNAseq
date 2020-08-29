simulate_data <- function(n, n.i, ngenes = 1, betas, phi, sigmab, seed.){
  # phi: dispersion parameter, e.g. phi= 4 corresponds to 0.25 dispersion in edgeR parameterization
  # n: sample size
  # n.i: repeated measurements
  # ngenes: number of genes
  # betas: $n_betas \times ngenes$ matrix with coefs
  # sigmab: random effects sd
  # This parameterization implies that size is the dispersion parameter 
  # and the variance is given by mu + mu^2/size.
  # So big dispersion/size implies poisson.
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
    runif(1)
  if (is.null(seed.)) 
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed.)
    RNGstate <- structure(seed., kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  
  time. <- seq(0, 2, length.out = n.i)
  group <- rep(1:2, each = n/2)-1
  data. <- data.frame(ID = rep(1:n, each = n.i), 
                      group = rep(group, each = n.i), 
                      time. = rep(time., n))
  X <- model.matrix(~ time. + group + time.:group, data = data.)
  n.g <- ngenes 
  DE.genes <- 1:(n.g/2)#indicator for DE genes
  DEstatus <- rep(0, ngenes)
  DEstatus[DE.genes] <- 1
  ln.lib.size <- rnorm(n*n.i, 0, 0.125)# increase to allow higher total counts
  betas.0 <- runif(ngenes, 0, 2) 
  phi <- 4 # this corresponds to 0.25 dispersion in edgeR parameterization
  betas.1 <- rep(1, ngenes)
  betas.2 <- rep(0, ngenes)
  betas.2[DE.genes] <- 2.4# set to 0 to be under the null hypothesis
  betas.3 <- rep(0, ngenes)
  betas.3[DE.genes] <- 2.4# set to 0 to be under the null hypothesis
  betas <- rbind(betas.0, betas.1, betas.2, betas.3)
  bi <- rnorm(n, 0, sigmab)
  eta <- X %*% betas + bi[data.$ID]
  mu <- exp(eta)
  counts <- matrix(rnbinom(ngenes*n*n.i, size = phi, 
                           mu = mu*exp(ln.lib.size)), 
                   ncol = ngenes, nrow = n*n.i)
  data.$ln.libsize <- log(rowSums(counts))
  data.$ln.lib.size <- ln.lib.size
  # This parameterization implies that size is the dispersion parameter 
  # and the variance is given by mu + mu^2/size.
  # So big dispersion/size implies poisson.
  data. <- cbind(data., counts) 
  data.
}

ngenes <- 10
data. <- simulate_data(n = 10, n.i = 5, ngenes = ngenes, 
                       betas = rbind(4, 1, 0, 0), 
                       phi = 4, sigmab = 2,
                       seed. = 1234)


file. <- "C:/RoulaTsonaka/Supervisions/MirkoSignorelli/Manuscripts/Paper2/BRIBIO_LaTex_Template-Jan2020/Revision1_August2020/github/"
write.table(data., paste(file., "Dataset.txt"))


