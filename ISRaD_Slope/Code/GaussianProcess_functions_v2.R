library(mlegp)
library(ggplot2)
#library(ggtree)
library(tidyr)
library(dplyr)
library(factoextra)
library(pbmcapply)
library(Rfast)

# define some functions for the analysis
# Calculate the variance covariance matrix for two sets of x-coordinates
calcK <- function(X1,X2,beta) {
  K <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1))
  for (i in 1:nrow(K)) {
    for (j in 1:ncol(K)) {
      K[i,j] <- exp(-beta*((X1[i]-X2[j])^2))
    }
  }
  return(K)
}

# Function to sample from an IGP object. The function is a rewrite of the Gaupro sample function. 
sampleIGP = function(n=1, model, x.star) {
  
  Mxx = matrix(NA, ncol=length(x.star), nrow=length(x.star))
  for (i in 1:length(x.star)) {
    for (j in i:length(x.star)) {
      v = exp(-model$theta()*(x.star[i] - x.star[j])^2)
      Mxx[i,j] = v
      Mxx[j,i] = v
    }
    Mxx[i,i = 1]
  }
  
  kxx <- Mxx + mean(model$nugget())
  
  Mx.xx = matrix(NA, nrow =length(model$X), ncol =length(x.star))
  for (i in 1:length(model$X)) {
    for (j in 1:length(x.star)) {
      Mx.xx[i,j] = exp(-sum(model$theta()*(model$X[i]-x.star[j])^2))
    }
  }
  kx.xx <- Mx.xx
  
  s2_hat = model$s2()
  
  K = matrix(NA, nrow =length(model$X), ncol =length(model$X))
  for (i in 1:length(model$X)) {
    for (j in 1:length(model$X)) {
      v = exp(-model$theta()*(model$X[i] - model$X[j])^2)
      K[i,j] = v
      K[j,i] = v    
    }
    K[i,i = 1]
  }
  K = K + diag(model$nugget(), length(model$X))
  Kinv = chol2inv(chol(K))
  
  covmat = s2_hat * (kxx - t(kx.xx)%*%Kinv%*%kx.xx)
  
  
  y = MASS::mvrnorm(n=n, mu = model$predict(matrix(x.star, ncol=1)), Sigma = covmat)
  
  return(t(y))
  
}

# Calculate realization of GP model, using x.xtar as independant variable
sampleGP = function(n=1, model, x.star) {
  # Calculate the covariance matrices
  # using the same x.star values as above
  x <- model$X
  z <- model$Z
  mu.hat <- model$mu
  k.xx <- calcK(x,x, model$beta)
  k.xxs <- calcK(x,x.star, model$beta)
  k.xsx <- calcK(x.star,x, model$beta)
  k.xsxs <- calcK(x.star,x.star, model$beta)
  
  z.s = mlegp::predict.gp(model, newData = as.data.frame(x=x.star))
  zbar = mlegp::predict.gp(model)
  
  sigma2 = diag(as.numeric(model$nugget^2), ncol(k.xx)) #model$sig2
  
  f.star.bar <-  z.s + k.xsx%*%solve(k.xx + sigma2)%*%(z-zbar) #
  cov.f.star <-  k.xsxs - k.xsx%*%solve(k.xx + sigma2)%*%k.xxs #solve(k.xx + sigma2)
 
  values <- try(mvtnorm::rmvnorm(n,  f.star.bar, cov.f.star,method = "svd")) #MASS::mvrnorm(n, f.star.bar, cov.f.star)
  
  if (inherits(values, "try-error")){
    out = values
  } else {
    out =  cbind(x=x.star,as.data.frame(t(values)))
  }
  return(
    out 
  )
}

# Function to get the signed maximum in a vector. 
signed.max=function(x) {
  if (!all(is.na(x))) {
    sign(x[which.max(abs(x))])*max(abs(x))
  } else {
    NA
  }
}


# Main function. takes a long table with columns containing, x, y and standard error values. The name of the columns are given using arguments idcol, 
compare_GP = function(longdata, xcol = "x", ycol = "y", errcol = "se", idcol = "id", nsamples = 100, verbose=1, engine ="mlegp", parallel = T, ncores = 1, position = F, keep_simulations = F) {
  require(parallel)
  print("new version")
  if (engine=="mlegp") {
    require(mlegp)
  }
  if (engine=="GauPro") {
    require(GauPro)
  }
  if (engine=="laGP_gaupro") {
    require(IGP)
    require(laGP)
    require(GauPro)
  }
  if (engine=="laGP") {
    require(laGP)
  }
  if (!(engine %in% c("mlegp", "GauPro","laGP", "laGP_gaupro"))) {
    stop("Engine must be mlegp, GauPro, laGP or laGP_gaupro")
  }
  
  ## Fitting GPs to data
  print("Fitting GP models...")
  
  models = lapply( #This fits GPs to all profile ids.
    split(longdata, longdata[,idcol]), #splitting data table on id -> makes a list with each element containing a single profile in each element
    function(data) {
      if (any(is.na(data[,errcol]))) { #Test whether a  nugget is available or not and fit the corresponding model
        if (engine == "mlegp") {
          out = try(mlegp::mlegp(data[,xcol, drop=T], data[,ycol, drop=T], nugget.known = 0, verbose=verbose))
        }
        if (engine == "GauPro") {
          out =try(GauPro(data[,xcol, drop=T], data[,ycol, drop=T], useGrad=T)$update(restarts=50))
        }
        if (engine == "laGP_gaupro") {
          out =try(IGP::IGP(package='laGP', data[,xcol, drop=T], data[,ycol, drop=T]))
        }
        if (engine == "laGP") {
          out =try({
            x.star = longdata %>% ungroup %>% tidyr::expand(!!as.name(xcol))
            d <- darg(list(mle=TRUE), data[,xcol])
            g <- garg(list(mle=TRUE), data[,ycol, drop = T])
            gpi <- newGP(data[,xcol], data[,ycol, drop=T], d=d$start, g=g$start, dK=TRUE)
            mle = mleGP(gpi, param = c("d"), tmax = 100000)
            updateGP(gpi, X, Z,verb = T)
            model <- predGP(gpi, matrix(pull(x.star), ncol = 1))
            deleteGP(gpi)
            return(model)
            })
          }
      } else {
        if (engine == "mlegp") {
          out =try(mlegp::mlegp(data[,xcol, drop=T], data[,ycol, drop=T], nugget= data[,errcol, drop=T], nugget.known = 1, verbose=verbose))
        }
        if (engine == "GauPro") {
          out =try(GauPro(data[,xcol, drop=T], data[,ycol, drop=T], nug= mean(data[,errcol, drop=T]), nug.est=F, verbose=verbose, useGrad=T,useC=T)$update(restarts=50))
        }
        if (engine == "laGP_gaupro") {
          out =try(IGP::IGP(package='laGP_gaupro', data[,xcol, drop=T], data[,ycol, drop=T], nug= data[,errcol, drop=T],estimate.nugget = F, lite=F))
        }
        if (engine == "laGP") {
          # THIS ONE DOES NOT WORK
          out =try({
            x.star = longdata %>% ungroup %>% tidyr::expand(!!as.name(xcol))
            d <- darg(list(mle=TRUE),data[,xcol] )
           # g <- garg(list(mle=T), data[,ycol, drop = T]*data[,errcol, drop = T])
            gpi <- newGP(data[,xcol], data[,ycol, drop=T], d=d$start, g=0.001, dK=TRUE)
            #mle <- jmleGP(gpi, drange = c(d$min,d$max), grange=c(g$min, g$max))
            mle = mleGP(gpi, param = c("d"), tmax = 100000)
            updateGP(gpi, X, Z,verb = T)
            model <- predGP(gpi, matrix(pull(x.star), ncol = 1))
            deleteGP(gpi)
            return(model)
          })
        }
      }
      return(out)
    }
    )
  #Data for which GP did not fit inherits the "try-error" class
  
  ## Simulating new data from the fitted GPs
  # Here we create a vector of all x values that are in the original dataset. Simulations will be calculated for all these values and data will be compared to the simulations at the corresponding x values. 
  x.star = longdata %>% ungroup %>% tidyr::expand(!!as.name(xcol))
  
  # for each model, we simulate new data at each x-value
  print("Simulating new data...")
  
  simuls = lapply(
    models, 
    function(model) {
      if(!inherits(model, "try-error")) { #check that model fitted
        if (engine == "mlegp") {
          out = try(sampleGP(nsamples, model, pull(x.star)))
        }
        if (engine == "GauPro") {
          out = try(cbind.data.frame(x=x.star, t(model$sample(pull(x.star), nsamples))))
        }
        if (engine == "laGP_gaupro") {
          out = try(cbind.data.frame(x=x.star,sampleIGP(n=nsamples, model, pull(x.star))))
          if(!inherits(out, "try-error")) {
          names(out) = c(xcol, 1:nsamples)
          }
        }
        if (engine == "laGP") {
          out = try({
            cbind.data.frame( x=x.star, t(Rfast::rmvnorm(nsamples, model$mean, model$Sigma)))
            })
          # if(!inherits(out, "try-error")) {
          #   names(out) = c(xcol, 1:nsamples)
          # }
        }
      } else {
        out = model
      }
      return(out)
    }
    )
  
  ## Compare original to simulated data  
  print("Comparing data to simulated data...")
  print("Calculating d")
  expanded_tab  = full_join(expand.grid(id = pull(unique(longdata[,idcol ])), Profondeur_num = pull(x.star)), longdata) %>% arrange(!!as.name(idcol), !!as.name(xcol))
  results = pbmclapply(names(models),
                function(p) {
                    measures = expanded_tab %>% filter(!!as.name(idcol)==p) %>% pull(!!as.name(ycol))
                    if (any(!is.na(measures))) {
                    GPfit = lapply(names(simuls), function(ns) {
                      s = simuls[[ns]]
                      if(!inherits(s, "try-error") ) {
                        diff = abs(eachrow(t(as.matrix(s[,-1])), measures,"-"))
                        d = mean(apply(diff, 1, max,na.rm=T), na.rm=T)
                        p = mean(s[apply(diff, 1, which.max),1], na.rm=T)
                      } else {
                        d=NA
                        p = NA
                      }
                      return(c(d=d,p=p))
                      }
                  )
                    GPfit = do.call(rbind, GPfit)
                    GPfit = data.frame(data = rep(p,length(simuls)), model = names(simuls), d = GPfit[,'d'], p = GPfit[,'p'])
                    } else {
                      GPfit = data.frame(
                        data = rep(p,length(simuls)), 
                        model = names(simuls),
                        d=NA, 
                        p=NA)
                    }
                    return(GPfit)
                    }, mc.cores = min(ncores,detectCores()-1))

   results = do.call(rbind.data.frame, results)
   #names(d)  = names(models)
   #rownames(d) = names(models)
  # p = do.call(rbind.data.frame, lapply(results, function(x) ifelse(!inherits(x, "try-error"),x$p,NA)))
   #names(p)  = names(models)
   #rownames(p) = names(models)
  
   if (keep_simulations) {
     output =  list(models = models,
                    simuls = simuls,
                    results = results)
   } else {
     output =  list(models = models,
                    results = results)
   }
   
  return(output)
  
}

bestGauPro = function(X, Z, ..., ntries = 1) {
  tmp = list()
  for (i in 1:ntries) {
    tmp[[i]] = try(GauPro(X, Z, ...))
  }
  best = which.max(unlist(lapply(tmp, function(m) m$loglikelihood())))
  return(tmp[[best]])
}
