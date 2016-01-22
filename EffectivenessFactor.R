int_weights <- c(0.04567809, 0.12589839, 0.13397916, 1 / 36 )

pr.effectiveness.factor <- function(surface, init.guess,
                                    t, p, C_total,
                                    preExp, Ea, alpha,
                                    intrPorosity, voidage, rParticle){
    require(rootSolve)
    #set orth coll function
    orthcoll <- set.ocf(surface, 
                        t, p, C_total,
                        preExp, Ea, alpha,
                        intrPorosity, voidage, rParticle)
    #solve for internal collocation points
    conc <- multiroot(orthcoll, init.guess)
    conc <- matrix(conc$root, ncol = 3, 
                   dimnames = list(NULL, c("hydrogen", "nitrogen", "ammonia")))
    #find kinetic parameters at the current conditions
    kin.const <- pr.kinetic.constant(preExp, Ea, t)
    eq.const <- pr.equilibrium.constant(t,p)
    #concentration matrix in coll points
    Xi <- rbind(conc, surface)
    colnames(Xi) <- c("hydrogen", "nitrogen","ammonia")
    #fugacities
    fgc.cf <- sapply(fugacity.list[c("hydrogen", "nitrogen","ammonia")],
                     function(x) do.call(x, args = list(t=t+T_abs, p=p))
    )
    fgc <- Xi%*%diag(fgc.cf)*p
    colnames(fgc) <- c("hydrogen", "nitrogen","ammonia")
    #reation rates at collocation points
    rates <- apply(fgc, 1, function(fgc) pr.reaction.rate.NH3(kin.const, eq.const, fgc, alpha))
    #return effectiveness factor
    return(3*(sum(int_weights*rates))/rates[length(rates)])
}

effectiveness.factor <- function(stream, reaction, catalyst, init.guess = rep(0.1, times = 9)){
    pr.effectiveness.factor(fraction(stream)[c("hydrogen","nitrogen","ammonia")], init.guess,
                            stream@conditions["temperature"], stream@conditions["pressure"], 
                            total.concentration(stream),
                            reaction@properties[["preExp"]], reaction@properties["activationEnergy"],
                            reaction@properties[["alpha"]],
                            catalyst@geometry[["intrPorosity"]],catalyst@geometry[["porosity"]],
                            catalyst@geometry[["pelletD"]]/2
                            )
}