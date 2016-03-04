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

#bed2
inlet2 <- recalculate.stream(inlet1, as.double(bedresult[nrow(bedresult),"1"]), TRUE)
inlet2@conditions <-  c(temperature=433, pressure=226)
bed2 <- new("Bed", inlet2, reaction1, catalyst1, 7.2)
bedode2 <- bed.ode.func(bed.db(bed2))
bedresult2 <- bed.calculate(bedode2, 0) 
                            #as.double(bedresult[nrow(bedresult),"1"]))

#bed 3
inlet3 <- recalculate.stream(inlet2, as.double(bedresult2[nrow(bedresult2),"1"]), TRUE)
inlet3@conditions <- c(temperature=415, pressure=226)
bed3 <- new("Bed", inlet3, reaction1, catalyst1, 7.8)
bedode3 <- bed.ode.func(bed.db(bed3))
bedresult3 <- bed.calculate(bedode3, 0)
                          #  as.double(bedresult2[nrow(bedresult2),"1"]))

reactorresult<-rbind(bed.summary(bedresult, bed1), bed.summary(bedresult2, bed2), bed.summary(bedresult3, bed3))
#test for CF's plant PFD
inlet_total <- new("Stream", 270, 97.0, 752901, c(hydrogen=0.6499,
                                                  nitrogen=0.2401,
                                                  ammonia=0.0286,
                                                  methane=0.0675, 
                                                  argon=0.0139))
inlet_bed1 <- inlet_total
inlet_bed1@mdot <- inlet_bed1@mdot*(1-0.31395)
quench_bed2 <- inlet_total
quench_bed2@mdot <- quench_bed2@mdot*0.34
out_bed2 <- recalculate.stream(inlet_total, 0.25, TRUE)
out_bed2@conditions[["temperature"]] <- 445
out_bed2@conditions[["pressure"]] <- inlet_total@conditions[["pressure"]]

debug(interchanger)
amminch <- interchanger(d.sh.inner = 1.133, d.sh.outer = 1.692,
                        dt.inner = 0.03175, dt.outer = 0.03175+2*0.00165, 
                        n.tube = 740, l.tube = 8.16, baff.spacing = 0.664, pitch = 0.03969, lambda = 43)
debug(amminch)
amminch(out_bed2@mdot,
        inlet_bed1@mdot,
        out_bed2@conditions[["temperature"]], inlet_bed1@conditions[["temperature"]], 95)
#test converter

ammconverter <- converter(b1 <- new("Bed", new("Stream",0, 0, 0, c(0,0,0,0,0)), reaction1, catalyst1, 4.75),
                       b2 <- new("Bed", new("Stream",0, 0, 0, c(0,0,0,0,0)), reaction1, catalyst1, 7.2),
                       b3 <- new("Bed", new("Stream",0, 0, 0, c(0,0,0,0,0)), reaction1, catalyst1, 7.8),
                       interchanger(d.sh.inner = 1.133, d.sh.outer = 1.692,
                                    dt.inner = 0.03175, dt.outer = 0.03175+2*0.00165, 
                                    n.tube = 740, l.tube = 8.16, baff.spacing = 0.664, pitch = 0.03969, lambda = 43))
debug(ammconverter)
ammconverter(inlet_bed1, list(NULL, quench_bed2))
