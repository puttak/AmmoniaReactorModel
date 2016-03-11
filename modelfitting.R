fitdata <- read.csv("freshcatalyst.csv")
#divide concentrations by 100
normalize <- function(x){
    xsum <- sum(x)
    return(x/xsum)
    }
conc_in <-t(apply(fitdata[c("D401_in_hydrogen",
                           "D401_in_argon","D401_in_nitrogen",
                           "D401_in_methane","D401_in_ammonia")],
                 1,
                 normalize))
conc_out <-t(apply(fitdata[c("D401_out_hydrogen",
                    "D401_out_argon","D401_out_nitrogen",
                    "D401_out_methane","D401_out_ammonia")],
          1,
          normalize))
processconditions <- data.frame(reactor_fr = fitdata[["D401_total_feed"]]-
                                    fitdata[["X1st_bed_quench"]]-
                                    fitdata[["X2nd_bed_quench"]],
                                quench1 = fitdata[["X1st_bed_quench"]],
                                quench2 = fitdata[["X2nd_bed_quench"]],
                                pressure = fitdata[["pressure_after_J401"]],
                                temperature = fitdata[["Feed_temperature"]]
                                )
fitdata <- cbind(processconditions, conc_in, conc_out)
rm(conc_in, conc_out, processconditions)
#inlets
inletstreams <- apply(fitdata,1,
                      function(x){
                          new("Stream", t = x[["temperature"]],
                              p = x[["pressure"]],
                              flowrate = x[["reactor_fr"]]*1000,
                              composition = c(
                                  hydrogen = x[["D401_in_hydrogen"]],
                                  nitrogen = x[["D401_in_nitrogen"]],
                                  ammonia = x[["D401_in_ammonia"]],
                                  methane = x[["D401_in_methane"]],
                                  argon = x[["D401_in_argon"]]
                              )
                              )
                      })
#2 bed quenches
quenches2 <- apply(fitdata,1,
                      function(x){
                          new("Stream", t = x[["temperature"]],
                              p = x[["pressure"]],
                              flowrate = x[["quench2"]]*1000,
                              composition = c(
                                  hydrogen = x[["D401_in_hydrogen"]],
                                  nitrogen = x[["D401_in_nitrogen"]],
                                  ammonia = x[["D401_in_ammonia"]],
                                  methane = x[["D401_in_methane"]],
                                  argon = x[["D401_in_argon"]]
                              )
                          )
                      })
#used within resuduial function
nitrogen_out <- fitdata[["D401_out_nitrogen"]]

#converter
ammconverter <- converter(new("Bed", 
                                new("reaction", c(preExp = 1.0, activationEnergy = 25000, alpha = 0.55)), 
                                new("Catalyst", c(intrPorosity=0.52, porosity=0.46, pelletD=0.0047)), 
                                15.97),
                          new("Bed", 
                                new("reaction", c(preExp = 1.0, activationEnergy = 25000, alpha = 0.55)), 
                                new("Catalyst", c(intrPorosity=0.52, porosity=0.46, pelletD=0.0048)),
                                45.97),
                          new("Bed", 
                                new("reaction", c(preExp = 1.0, activationEnergy = 25000, alpha = 0.55)), 
                                new("Catalyst", c(intrPorosity=0.52, porosity=0.46, pelletD=0.0042)), 
                                74.03),
                          interchanger(d.sh.inner = 1.133, d.sh.outer = 1.692,
                                       dt.inner = 0.03175, dt.outer = 0.03175+2*0.00165, 
                                       n.tube = 740, l.tube = 8.16, baff.spacing = 0.664, 
                                       pitch = 0.03969, lambda = 43))


residualf <- function(par, observed, inletstreams, quenches2){
    require(foreach)
    require(doParallel)
    environment(ammconverter)$bed1@reaction@properties[["preExp"]] <- par$preExp1
    environment(ammconverter)$bed1@reaction@properties[["activationEnergy"]] <- par$Ea1
    environment(ammconverter)$bed2@reaction@properties[["preExp"]] <- par$preExp2
    environment(ammconverter)$bed2@reaction@properties[["activationEnergy"]] <- par$Ea2
    environment(ammconverter)$bed3@reaction@properties[["preExp"]] <- par$preExp2
    environment(ammconverter)$bed3@reaction@properties[["activationEnergy"]] <- par$Ea2
    if (length(inletstreams)!= length(quenches2))
        stop("inconsistent inlets and quenches")
    if (length(nitrogen_out) != length(inletstreams))
        stop("Different length for streams and experimantal points")
    # cl <- makeCluster(2)
    # registerDoParallel(cl)
    result <- foreach(i=1:length(inletstreams)) %do%{
       res <- ammconverter(inletstreams[[i]], list(NULL, quenches2[[i]]))
       res <- sapply(res, function(x) x[nrow(x), "1"])
       names(res) <- NULL
       bed3_out <- recalculate.stream(recalculate.stream(recalculate.stream(inletstreams[[i]], res[[1]], T),
                                                         res[[2]], T),
                                      res[[3]], T)
       as.double(fraction(bed3_out, "nitrogen"))
       
    }
    # stopCluster(cl)
    return((nitrogen_out - unlist(result))*10000)
}
residualf(parlist(preExp1=10^14.7102, preExp2 = 10^14.7102, Ea1 = 1.635e5, Ea2=1.635e5),
          nitrogen_out,
          inletstreams,
          quenches2)
library(minpack.lm)
kinparams <- nls.lm(par=list(preExp1=10^14.7102, preExp2 = 10^14.7102, Ea1 = 1.635e5, Ea2=1.635e5),
                 fn = residualf, observed = nitrogen_out,
                 inletstreams = inletstreams,
                 quenches2 = quenches2,
                 control = nls.lm.control(nprint=1))


x_real <- seq(1,20,length.out = 10)
true_pars <- list(a=1, b=2, c=2 )
foo <- function(par){
    return(par$a + par$b*x_real + x_real^2*par$c)
}
y_real <- foo(true_pars)
y_noisy <- y_real + rnorm(length(y_real), 0, 0.3)
residualfunc <- function(pars) y_noisy - foo(pars) 
kinparams <- nls.lm(par = list(a=1, b=1, c=2),fn= residualfunc,
                        control = nls.lm.control(nprint=1))


###### example 1
## values over which to simulate data
x <- seq(0,5,length=100)
## model based on a list of parameters
getPred <- function(parS, xx) parS$a * exp(xx * parS$b) + parS$c
## parameter values used to simulate data
pp <- list(a=9,b=-1, c=6)
## simulated data, with noise
simDNoisy <- getPred(pp, nls.out,x) + rnorm(length(x),sd=.1)
## plot data
plot(x,simDNoisy, main="data")
## residual function
residFun <- function(p, observed, xx) observed - getPred(p,xx)
## starting values for parameters
parStart <- list(a=3,b=-.001, c=1)
## perform fit
nls.out <- nls.lm(par=parStart, fn = residFun, observed = simDNoisy,
                  xx = x, control = nls.lm.control(nprint=1))

