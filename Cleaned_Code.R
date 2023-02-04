#' Code to run analysis for manuscript: Individual-level biotic interactions and species distribution models
#' Constants, data and intial values are stored in .txt files
#'########################################
#' The following describes all parameters seen in the NIMBLE model, in order of appearance:
#'  alpha0.a     intercept for detection function for species a; probability of detecting an indiviudal whose activity center is directly over trap/net location  
#'  sigma.a      scale parameter of detection function for species a
#'  alpha0.b     intercept for detection function for species b
#'  sigma.b      scale parameter of detection function for species b
#'  sig.inhib    strength of inhibition between two species
#'  n.years      years of data
#'  B0.a[t]      intercept in year t in log-linear function of yearly abundance for species a
#'  B0.b[t]      intercept in year t in log-linear function of yearly abundance for species b
#'  B1.a         relationship between environmental variables and abundance for species a
#'  B1.b         relationship between environmental variables and abundance for species b
#'  n.sites      how many sites are present in data
#'  FF           how many "mini-pixels" are present inside each site
#'  M.a          total number of individuals, real + agumented being analyzed
#'  exp.dist     expected distance between each activity center of species a and each mini-pixel center
#'  owned.ind    expected area not "owned" by each individual in species a in each mini pixel
#'  owned.f      expected available proportion for species b in each mini pixel
#'  Avail        expected proportion of each site availabe for species b
#'  lambda.a     expected abundance of species a at site q in time t
#'  lambda.b     expected abundance of species b at site q in time t
#'  psi          proportion of individuals in M.a or M.b at each site that are real 
#'  obs.a        number of individuals in species a seen during sample 
#'  s.a[,1:2,,]  activity center x and y coordinates for each individual in species a
#'  z.a          binary indicator (1 = individual is real, 0 = not part of population)
#'  like.za      for WAIC calculation (likelihood of z's for species a)
#'  d.a          distance between traps and activity centers of species a
#'  p.a          probability of catching each individual at each trap 
#'  y.a          binary indicator of capture for each individual at each trap (1  = yes, 0 = no)
#'  likeya       for WAIC calculation (likelihood of detections for species a)
#'  zeros.trick2 computationally efficient way of dealing with all the augmented non-detections
#'  s.b[,1:2,,]  activity center x and y coordinates for each individual in species b
#'  ones.b       forcing JAGS/NIMBLE to respect soft-core strauss process
#'  like.onesb   for WAIC calculation (likelihood of one's being true given inhibition)
#'  dist.b       distance between all species a activity centers and activity centers of species b
#'  inhib        inhibition strength based on all activity centers of real individuals of species a 
#'  neighbors    product of all inhibition of all individuals of species a 
#'  z.b          binary indicator (1 = individual is real, 0 = not part of population)
#'  like.zb      for WAIC calculation (likelihood of z's for species b)
#'  d.b          distance between traps and activity centers of species b
#'  p.b          probability of catching each individual at each trap 
#'  y.b          binary indicator of capture for each individual at each trap (1  = yes, 0 = no)
#'  likeyb       for WAIC calculation (likelihood of detections for species b)
#'  zeros.trick2.b computationally efficient way of dealing with all the augmented non-detections
#'  N            Abundance in each year
#'  N.a          Abundance of species a in each year
#'  N.b          Abundance of species b in each year
#'  a            Abundance of species a in each year at each site
#'  b            Abundance of species b in each year at each site
#'  like         for WAIC calculation (total log-likelihood for this iteration) 


library(nimble)
Test_nimble <- nimbleCode( {
  
  alpha0.a ~dnorm(0,.01)
  sigma.a ~dunif(0,200)
  alpha0.b ~dnorm(0,.01) 
  sigma.b ~dunif(0,200)
  sig.inhib ~ dunif(0, 50)
  
  for (t in 1:n.years){
    B0.a[t] ~ dunif(-10,10)
    B0.b[t] ~ dunif(-10,10)
  }
  
  B1.a ~ dunif(-10,10)
  B1.b ~ dunif(-10,10)
  
  for (t in 1:n.years){ #years of data
    for (q in 1:n.sites){ #run through each site
      for (f in 1:FF){ #mini pixels to estimate availability
        for (a in 1:M.a){ #M.a will be the same size in every year
          exp.dist[f,a,q,t] <- pow(pow(pix[f,1,q]-s.a[a,1,q,t],2) + pow(pix[f,2,q]-s.a[a,2,q,t],2),0.5)
          owned.ind[f,a,q,t] <- 1-exp((-exp.dist[f,a,q,t]*exp.dist[f,a,q,t])/(sig.inhib*sig.inhib))*z.a[a,q,t]
        } #end a
        
        owned.f[f,q,t] <- prod(owned.ind[f,1:M.a,q,t])
      } #end f
      
      Avail[q,t] <- sum(owned.f[1:FF,q,t])/FF #proportion available in this grid cell this year
      
      lambda.a[q,t] <- exp(B0.a[t] + B1.a*pca[q,t])
      lambda.b[q,t] <- exp(B0.b[t] + B1.b*pca[q,t])*Avail[q,t]
      
      psi[1,q,t] <- lambda.a[q,t]/M.a 
      psi[2,q,t] <- lambda.b[q,t]/M.b
      
      for (i in 1:obs.a){ #the "real" birds of species a 
        s.a[i,1,q,t] ~dunif(xmins[q],xmaxs[q])  #x coord
        s.a[i,2,q,t] ~dunif(ymins[q],ymaxs[q]) #y coord
        
        
        z.a[i,q,t]  ~dbern(psi[1,q,t])
        like.za[i,q,t] <- dbinom(z.a[i,q,t],1, psi[1,q,t], log = T)
        
        for (j in 1:20){ #traps at this site
          d.a[i,j,q,t] <- pow(pow(s.a[i,1,q,t]-trap[j,1,q],2) + pow(s.a[i,2,q,t]-trap[j,2,q],2),0.5)  #distance between traps and locations
          for (k in 1:K){ #capture occasions
            p.a[i,j,k,q,t] <- alpha0.a*exp(-d.a[i,j,q,t]*d.a[i,j,q,t]/(sigma.a^2))*z.a[i,q,t]*open[j,k,q,t]  #detection prob
            y.a[i,j,k,q,t] ~ dbern(p.a[i,j,k,q,t])   #captured? based on probability of detection and # occasions
            like.ya[i,j,k,q,t] <- dbinom(y.a[i,j,k,q,t],1,  p.a[i,j,k,q,t], log = T)
          } #end k
          likeya_1[i,j,q,t] <- sum(like.ya[i,j,1:K,q,t])
        } #end j
        likeya[i,q,t] <- sum(likeya_1[i,1:20,q,t]) # loglike of y across traps and days
      } #end i (real a)
      
      for(i in ((obs.a+1):M.a)){ #augmented birds species A
        s.a[i,1,q,t] ~dunif(xmins[q],xmaxs[q])  #x coord of the activity center
        s.a[i,2,q,t] ~dunif(ymins[q],ymaxs[q])  # y coord
        
        for (k in 1:K){
          zeros.trick2[i,k,q,t] ~ dbern(1-(prod(1-p.a[i,1:20,k,q,t])))
        }
        
        z.a[i,q,t]  ~dbern(psi[1,q,t])
        like.za[i,q,t] <- dbinom(z.a[i,q,t],1,psi[1,q,t], log = T)
        
        for (j in 1:20){
          d.a[i,j,q,t] <- pow(pow(s.a[i,1,q,t]-trap[j,1,q],2) + pow(s.a[i,2,q,t]-trap[j,2,q],2),0.5)  #distance between traps and locations
          for (k in 1:K){
            p.a[i,j,k,q,t] <- alpha0.a*exp(-d.a[i,j,q,t]*d.a[i,j,q,t]/(sigma.a^2))*z.a[i,q,t]*open[j,k,q,t]
          }}
      } #end augmented species A
      
      
      for (i in 1:obs.b){   #real b 
        s.b[i,1,q,t] ~dunif(xmins[q],xmaxs[q])  #coords boxed in by site bounds
        s.b[i,2,q,t] ~dunif(ymins[q],ymaxs[q])
        
        ones.b[i,q,t] ~dbern(neighbors[i,q,t])
        like.onesb[i,q,t] <- dbinom(ones.b[i,q,t],1, neighbors[i,q,t], log = T)
        
        for (a in 1:M.a){ #species A at this site
          dist.b[i,a,q,t] <- pow(pow(s.b[i,1,q,t]-s.a[a,1,q,t],2) + pow(s.b[i,2,q,t]-s.a[a,2,q,t],2),0.5) #distance to all locations species a
          #inhib[i,a,q,t] <- 1 #for removing interaction
          inhib[i,a,q,t] <- 1-exp((-dist.b[i,a,q,t]*dist.b[i,a,q,t])/(sig.inhib*sig.inhib))*z.a[a,q,t]*z.b[i,q,t] #half normal inhibition (~= soft-core strauss)
        }
        
        neighbors[i,q,t] <- prod(inhib[i,1:M.a,q,t]) #product of probabilities of all neighbors
        z.b[i,q,t]  ~dbern(psi[2,q,t])
        like.zb[i,q,t] <- dbinom(z.b[i,q,t],1, psi[2,q,t], log = T)
        
        for (j in 1:20){ #traps at this site
          d.b[i,j,q,t] <- pow(pow(s.b[i,1,q,t]-trap[j,1,q],2) + pow(s.b[i,2,q,t]-trap[j,2,q],2),0.5)  #distance between traps and locations
          for (k in 1:K){
            p.b[i,j,k,q,t] <- alpha0.b*exp(-d.b[i,j,q,t]*d.b[i,j,q,t]/(sigma.b^2))*z.b[i,q,t]*open[j,k,q,t]  #detection prob
            y.b[i,j,k,q,t] ~ dbern(p.b[i,j,k,q,t])   #captured? based on probability of detection and # occasions
            like.yb[i,j,k,q,t] <- dbinom(y.b[i,j,k,q,t], 1, p.b[i,j,k,q,t], log = T)
          } #end k
          likeyb_1[i,j,q,t] <- sum(like.yb[i,j,1:K,q,t])
        } #end j
        likeyb[i,q,t] <- sum(likeyb_1[i,1:20,q,t]) # loglike of y across traps and days
        
      }
      
      
      for (i in (obs.b+1):M.b){ #augmented species b 
        s.b[i,1,q,t] ~dunif(xmins[q],xmaxs[q])  #x coord of the activity center
        s.b[i,2,q,t] ~dunif(ymins[q],ymaxs[q])  # y coord
        
        ones.b[i,q,t] ~ dbern(neighbors[i,q,t])
        like.onesb[i,q,t] <- dbinom(ones.b[i,q,t], 1, neighbors[i,q,t], log = T) 
        
        for (a in 1:M.a){ #species A at this site
          dist.b[i,a,q,t] <- pow(pow(s.b[i,1,q,t]-s.a[a,1,q,t],2) + pow(s.b[i,2,q,t]-s.a[a,2,q,t],2),0.5) #distance to all locations species a
          #inhib[i,a,q,t] <- 1 #for removing interaction
          inhib[i,a,q,t] <- 1-exp((-dist.b[i,a,q,t]*dist.b[i,a,q,t])/(sig.inhib*sig.inhib))*z.a[a,q,t]*z.b[i,q,t] #half normal inhibition (~= soft-core strauss)
        }
        
        neighbors[i,q,t] <- prod(inhib[i,1:M.a,q,t]) #product of probabilities of all neighbors
        z.b[i,q,t]  ~dbern(psi[2,q,t])
        like.zb[i,q,t] <- dbinom(z.b[i,q,t],1, psi[2,q,t], log = T)
        
        for (k in 1:K){
          zeros.trick2.b[i,k,q,t] ~ dbern(1-(prod(1-p.b[i,1:20,k,q,t])))
        }
        
        for (j in 1:20){ #traps at this site
          d.b[i,j,q,t] <- pow(pow(s.b[i,1,q,t]-trap[j,1,q],2) + pow(s.b[i,2,q,t]-trap[j,2,q],2),0.5)  #distance between traps and locations
          for (k in 1:K){
            p.b[i,j,k,q,t] <- alpha0.b*exp(-d.b[i,j,q,t]*d.b[i,j,q,t]/(sigma.b^2))*z.b[i,q,t]*open[j,k,q,t] 
          } # end K
        } #end J
      } # end aug species B
      
    } # end site loop
    
    N[t] <- N.a[t] + N.b[t] #all species together
    N.b[t] <- sum(z.b[1:M.b,1:n.sites,t]) #all alive for species 2
    N.a[t] <- sum(z.a[1:M.a, 1:n.sites,t])
    
    for (q in 1:n.sites){
      a[q,t] <- sum(z.a[1:M.a,q,t])
      b[q,t] <- sum(z.b[1:M.b,q,t])
    } #end q
  }#end t
  
  ## For WAIC likelihood estimation 
  for (q in 1:n.sites){
    for(t in 1:n.years){
      like[q,t] <- sum(like.za[1:M.a,q,t])+sum(likeya[1:obs.a, q,t])+
        sum(like.zb[1:M.b,q,t])+sum(likeyb[1:obs.b, q,t]) 
    }
  }
  
  
})


params <- c("B0.b",
            "B0.a", "B1.a", "B1.b",
            "alpha0.a", "sigma.a",
            "alpha0.b", "sigma.b",
            "N.b", "N.a","N",
            "sig.inhib","a", "b", "psi"
)

constants <- dget("Strauss_constants2_minpix.txt")
data.test <- dget("Strauss_data2.txt")
inits <- dget("Strauss_inits2.txt")

### If running WAIC calculations, change params to:
#params <- c("like.za", "likeya", "like.zb", "likeyb")

options(scipen=999)
library(parallel)
cl <- makeCluster(3)
clusterExport(cl = cl, varlist = c("constants", "data.test", "inits", "params", 'Test_nimble'))
system.time(outs <- clusterEvalQ(cl = cl, {
  library(nimble)
  library(coda)
  testy.mod <- nimbleModel(code = Test_nimble,
                           constants = constants,
                           data = data.test,
                           inits = inits)
  testy.mod$initializeInfo()
  mcmcSCR<-configureMCMC(testy.mod,monitors=params, print = T, warnNoSamplerAssigned = T)
  SCRMCMC <- buildMCMC(mcmcSCR)
  Cmodel <- compileNimble(testy.mod)
  CompSCRMCMC <- compileNimble(SCRMCMC, project = testy.mod)
  CompSCRMCMC$run(niter = 50000, nburnin = 30000, thin = 10)
  #if you run this in your console it will say "null". 
  return(as.mcmc(as.matrix(CompSCRMCMC$mvSamples)))
})
)

#20K iterations w/ WAIC took 2.3 hr
#including WAIC but w/ density surface monitoring took ~ 10 hours


library(coda)
birds.mod.nimble <- mcmc.list(outs)
mod <- summary(birds.mod.nimble)

### WAIC calc #####
### This only works if you ran the parameters commented out in line 227
### With interaction 
w_int <- as.matrix(mod)
like <- array(NA, dim = c(nrow(w_int), 19, 4))
for(t in 1:4){
  for(q in 1:19){
    pat_za <- paste("like.za[", 1:40, ", ", q, ", ", t, "]", sep = "")
    pat_zb <- paste("like.zb[", 1:20, ", ", q, ", ", t, "]", sep = "")
    pat_ya <- paste("like.ya[", 1:18, ", ", q, ", ", t, "]", sep = "")
    pat_yb <- paste("like.yb[", 1:9, ", ", q, ", ", t, "]", sep = "")
    like[,q,t] <- rowSums(w_int[,colnames(w_int) %in% pat_za]) + rowSums(w_int[,colnames(w_int) %in% pat_ya]) + rowSums(w_int[,colnames(w_int) %in% pat_zb])+ rowSums(w_int[,colnames(w_int) %in% pat_yb])
  }
}
fbar <- colMeans(exp(like))
Pw <- sum(apply(like, 2, var))
WAIC_ish <- -2*sum(log(fbar))+2*Pw
## Pw = 606.821, WAIC = 4376.523

### Without interaction 
## This only works if you ran the parameters commented out in line 227 and turned on lines 139 and 170 (which turns off inhibition in model)
w_noint <- as.matrix(mod)
like_no <- array(NA, dim = c(nrow(w_noint), 19, 4))
for(t in 1:4){
  for(q in 1:19){
    pat_za <- paste("like.za[", 1:40, ", ", q, ", ", t, "]", sep = "")
    pat_zb <- paste("like.zb[", 1:20, ", ", q, ", ", t, "]", sep = "")
    pat_ya <- paste("like.ya[", 1:18, ", ", q, ", ", t, "]", sep = "")
    pat_yb <- paste("like.yb[", 1:9, ", ", q, ", ", t, "]", sep = "")
    like_no[,q,t] <- rowSums(w_noint[,colnames(w_noint) %in% pat_za]) + rowSums(w_noint[,colnames(w_noint) %in% pat_ya]) + rowSums(w_noint[,colnames(w_noint) %in% pat_zb])+ rowSums(w_noint[,colnames(w_noint) %in% pat_yb])
  }
}
fbar_no <- colMeans(exp(like_no))
Pw_no <- sum(apply(like_no, 2, var))
WAIC_ish_no <- -2*sum(log(fbar_no))+2*Pw_no
## Pw = 599.06, WAIC = 4362.348
