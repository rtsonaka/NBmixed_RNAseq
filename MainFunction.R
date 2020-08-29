if(FALSE){
formula.nbmixed <- "time. + group + time.:group" 
formula.nbmixed.0 <- "time." 
data <- data.
nboots <- 5
n.genes <- 10
gene.nams <- as.character(1:10)
nAGQ <- nAGQ. <- 21
test.coefs <- c("group", "time.:group")
}

longRNAseq <- function(formula.nbmixed, formula.nbmixed.0, 
                       data, nboots, test.coefs,
                       n.genes, gene.nams, nAGQ){
  
  library(GLMMadaptive)
  distr.genes <- numeric(n.genes) 
  var.genes <- disp.genes <- numeric(n.genes) 
  lrt.genes <- plrt.genes <- numeric(n.genes)
  plrtboot.genes <- numeric(n.genes)

  for(k in 1:n.genes){
    
    lrt.boot <- numeric(nboots)
    
    data$gene <- data[, as.character(gene.nams[k])] 
    form <- as.formula(paste("gene ~ offset(ln.lib.size)", 
                             formula.nbmixed, sep = "+"))
    form.0 <- as.formula(paste("gene ~ offset(ln.lib.size)", 
                             formula.nbmixed.0, sep = "+"))
    
    model.1 <- try(mixed_model(form, 
                                  data, random = ~ 1 | ID, 
                                  family = GLMMadaptive::negative.binomial(), 
                                  nAGQ = nAGQ))
    distr.genes[k] <- "NB"
    n.betas <- length(fixef(model.1))
    
    if(inherits(model.1, "try-error") || 
       if(!inherits(model.1, "try-error")){
         exp(model.1$phis) > 1000
         }){
      model.1 <- try(mixed_model(fixed = form, 
                                 random = ~ 1 | ID, data = data, 
                                 family = poisson(), nAGQ = nAGQ))
      model.0 <- try(mixed_model(fixed = form.0, 
                                 random = ~ 1 | ID, data = data, 
                                 family = poisson(), nAGQ = nAGQ))
      distr.genes[k] <- "P"
      disp.genes[k] <- NA
    }else{
      
      model.0 <- try(mixed_model(fixed = form.0, 
                                 random = ~ 1 | ID, data = data, 
                                 family = GLMMadaptive::negative.binomial(), 
                                 nAGQ = nAGQ))
      disp.genes[k] <- exp(model.1$phis)
      names(disp.genes[k]) <- ""
    }
    
    if(!inherits(model.1, "try-error")){
      s <- try(summary(model.1))
      if(!inherits(s, "try-error")){
        coefs. <- fixef(model.1)
        var.genes[k] <- model.1$D[1]
        se. <- s$ coef_table[ , "Std.Err"]
        gr. <- try(anova(model.0, model.1))
        
        if(!inherits(gr., "try-error")){
          lrt.genes[k] <-  gr.$ LRT
          plrt.genes[k] <-  gr.$ p.value
        } else{
          lrt.genes[k] <- NA
          plrt.genes[k] <-  NA
        }
        #bootstrap
        model.1.boot <- model.1
        model.1.boot$coefficients[test.coefs] <- c(0,0)#change!
        
        form.boot <- as.formula(paste("Resp ~ offset(ln.lib.size)", 
                                 formula.nbmixed, sep = "+"))
        form.boot.0 <- as.formula(paste("Resp ~ offset(ln.lib.size)", 
                                   formula.nbmixed.0, sep = "+"))
        
        for(seed_boot in 1:nboots){      
          data$Resp <- simulate(model.1.boot, nsim = 1, 
                                seed = 12345 + seed_boot, 
                                new_RE = TRUE)
          if(distr.genes[k] == "P"){
            model.1boot <- try(mixed_model(fixed = form.boot, 
                                           random = ~ 1 | ID, data = data, 
                                           family = poisson(), nAGQ = nAGQ))
            model.0boot <- try(mixed_model(fixed = form.boot.0, 
                                           random = ~ 1 | ID, data = data, 
                                           family = poisson(), nAGQ = nAGQ))
          }else{
            model.1boot <- try(mixed_model(fixed = form.boot, 
                                           random = ~ 1 | ID, data = data, 
                                           family = GLMMadaptive::negative.binomial(), 
                                           nAGQ = nAGQ))
            model.0boot <- try(mixed_model(fixed = form.boot.0, 
                                           random = ~ 1 | ID, data = data, 
                                           family = GLMMadaptive::negative.binomial(), 
                                           nAGQ = nAGQ))
          }
          
          gr.boot <- try(anova(model.0boot, model.1boot))
          if(!inherits(gr.boot, "try-error")){
            lrt.boot[seed_boot] <-  gr.boot$ LRT
          } else{
            lrt.boot[seed_boot] <- NA
          }
          plrtboot.genes[k] <- if(all(is.na(lrt.boot))){
          NA
        }else{
          (sum(lrt.boot >= lrt.genes[k], na.rm = TRUE) + 1)/ (1 + sum(!is.na(lrt.boot)))
        }
      }
    }
  }
  }    
  list(Gene = gene.nams, 
       var.b = var.genes, 
       disp = disp.genes,
       plrt = plrt.genes,
       plrt.boot = plrtboot.genes)
}

