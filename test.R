conc <- c(hydrogen=0.2219,
          nitrogen=0.6703,
          ammonia=0.0276,
          methane=0.0546, 
          argon=0.0256)
conc2 <- c(0.6703,0.2219,0.0276,0.0546,0.0256)
inlet1 <- new("Stream", 385, 226, 242160, conc)
inlet2 <- new("Stream", 385, 226, 242160, conc2)

reaction1 <- new("reaction", c(preExp = 1e7, activationEnergy = 48500, alpha = 0.52))

bed1 <- new("Bed", inlet1, reaction1, 4.75)

oc <- set.ocf(c(0.2219, 0.6703, 0.0276),
              inlet1@conditions["temperature"], inlet1@conditions["pressure"], total.concentration(inlet1),
              reaction1@properties[["preExp"]], reaction1@properties[["activationEnergy"]],
              reaction1@properties[["alpha"]],
              catalyst1@geometry[["intrPorosity"]],catalyst1@geometry[["porosity"]],
              catalyst1@geometry[["pelletD"]]/2)
res1 <- oc(c(0.15, 0.17, 0.19, 0.55, 0.60, 0.62,0.013, 0.019, 0.022))

debug(pr.effectiveness.factor)
ef <- pr.effectiveness.factor(c(0.2219, 0.6703, 0.0276),
                              c(rep(0.1, times = 9)),
                              inlet1@conditions["temperature"], inlet1@conditions["pressure"], total.concentration(inlet1),
                              reaction1@properties[["preExp"]], reaction1@properties[["activationEnergy"]],
                              reaction1@properties[["alpha"]],
                              catalyst1@geometry[["intrPorosity"]],catalyst1@geometry[["porosity"]],
                              catalyst1@geometry[["pelletD"]]/2
)