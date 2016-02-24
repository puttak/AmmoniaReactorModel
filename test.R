conc <- c(hydrogen=0.6703,
          nitrogen=0.2219,
          ammonia=0.0276,
          methane=0.0546, 
          argon=0.0256)
inlet1 <- new("Stream", 385, 226, 242160, conc)
reaction1 <- new("reaction", c(preExp = 10^14.7102, activationEnergy = 1.635e5, alpha = 0.55))
catalyst1 <- new("Catalyst", c(intrPorosity=0.52, porosity=0.46, pelletD=0.0057))
#inlet2 <- new("Stream", 385, 226, 242160, conc2)

#test flow methods
flow.rate(inlet1)
flow.rate(inlet1, "mole")
fraction(inlet1)
heat.capacity(inlet1, hc.un="mole")
heat.capacity(inlet1, hc.un="mass")
fugacity(inlet1)
total.concentration(inlet1)

pr.reaction.rate.NH3(pr.kinetic.constant(reaction1@properties[["preExp"]],
                                                     reaction1@properties[["activationEnergy"]],
                                                     inlet1@conditions[["temperature"]]
                                                     ),
                     pr.equilibrium.constant(inlet1@conditions[["temperature"]]),
                     fugacity(inlet1)[1:3],
                     0.55
                     )
RNH3(inlet1, reaction1)

oc <- set.ocf(c(hydrogen=0.6703,nitrogen=0.2219, ammonia=0.0276),
              inlet1@conditions["temperature"], inlet1@conditions["pressure"], total.concentration(inlet1),
              reaction1@properties[["preExp"]], reaction1@properties[["activationEnergy"]],
              reaction1@properties[["alpha"]],
              catalyst1@geometry[["intrPorosity"]],catalyst1@geometry[["porosity"]],
              catalyst1@geometry[["pelletD"]]/2)
debug(oc)
oc(rep(0.1, times = 9))

pr.effectiveness.factor(c(hydrogen=0.6703,nitrogen=0.2219, ammonia=0.0276),
                          c(rep(0.6703, times = 3), rep(0.2219, times = 3),rep(0.0276, times = 3)),
                          inlet1@conditions["temperature"], inlet1@conditions["pressure"], total.concentration(inlet1),
                          reaction1@properties[["preExp"]], reaction1@properties[["activationEnergy"]],
                          reaction1@properties[["alpha"]],
                          catalyst1@geometry[["intrPorosity"]],catalyst1@geometry[["porosity"]],
                          catalyst1@geometry[["pelletD"]]/2
                          )
init.guess <- rep(0.1, times = 9)
pr.effectiveness.factor(fraction(inlet1)[c(2,1,3)], init.guess,
                        inlet1@conditions["temperature"], inlet1@conditions["pressure"], 
                        total.concentration(inlet1),
                        reaction1@properties[["preExp"]], reaction1@properties["activationEnergy"],
                        reaction1@properties[["alpha"]],
                        catalyst1@geometry[["intrPorosity"]],catalyst1@geometry[["porosity"]],
                        catalyst1@geometry[["pelletD"]]/2
)

debug(pr.effectiveness.factor)
effectiveness.factor(inlet1, reaction1, catalyst1)

#bed
bed1 <- new("Bed", inlet1, reaction1, catalyst1, 4.75)
#bed ode function
bedode <- bed.ode.func(bed.db(bed1))
#integrate bed ode
bedresult <- bed.calculate(bedode, 0)

bed.summary(bedresult, bed1)

