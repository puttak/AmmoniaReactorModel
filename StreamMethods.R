mass.to.mole <- function(mass, element){
    mole <- mass/molWeight[[element]]*1000
    return(mole)
}
#Flow rate, returns mass or molar flow rate of stream
setGeneric(
    name = "flow.rate",
    def = function(object, ...){standardGeneric("flow.rate")}
)
setMethod("flow.rate",
          "Stream",
          definition = function(object, type = "mass"){
              if (type == "mass")
                  #mass rate in kg/s
                  fr <- sum(object@mdot)
              else if (type == "mole")
                  #molar rate in mole/s
                  fr <- sum(sapply(elementNames, function(x) mass.to.mole(object@mdot[[x]],x)) )
              return(fr)
          }
)
setGeneric(
    name = "fraction",
    def = function(object, ...){standardGeneric("fraction")}
)
setMethod("fraction",
          "Stream",
          definition = function(object, element = NULL, type = "mole"){
              if (type == "mole")
                  to.call <- function(x) {mass.to.mole(object@mdot[[x]],x)/flow.rate(object, type = "mole")}
              else if (type == "mass")
                  to.call <- function(x) {object@mdot[[x]]/flow.rate(object)}
              else 
                  stop("Wrong type")
              if (is.null(element))
                  frac <-sapply(names(object@mdot), to.call)
              else 
                  frac <- to.call(element)
              return(frac)
          }
)

setMethod("heat.capacity",
          "Stream",
          definition = function(object, hc.un = "mass", 
                                t = object@conditions[["temperature"]],
                                p = object@conditions[["pressure"]]){
                  hc <- sapply(names(object@mdot), 
                               function(x) heat.capacity(x, t, p, hc.un))
                  frac <- sapply(names(object@mdot),
                                 function(x) fraction(object, x, hc.un))
                  hc <- sum(hc*frac)
                  return(hc)
          }
)
#returns fugacity coefs at stream conditions
setMethod("fugacity.cf",
          "Stream",
          definition = function(object, element = NULL,
                                t = object@conditions[["temperature"]],
                                p = object@conditions[["pressure"]], 
                                un = c("C", "atm")){
              to.call <- function(x) fugacity.cf(x, t, p, un)
              if (is.null(element))
                  fgc <- sapply(names(object@mdot), to.call)
              else
                  fgc <- to.call(element)
              return(fgc)
          }
    
)
#returns fugacities of components at stream conditions
setGeneric("fugacity",
           def = function(object, ...){standardGeneric("fugacity")}
)
setMethod("fugacity",
          "Stream",
          definition = function(object, element = NULL,
                                t = object@conditions[["temperature"]],
                                p = object@conditions[["pressure"]], 
                                un = c("C", "atm") ){
              fgc <- fugacity.cf(object, element, t, p, un)*fraction(object, element)*p
              return(fgc)
          }
)

#returns concentration of stream kmol/m3 
setGeneric("total.concentration",
           def = function(object,...){standardGeneric("total.concentration")}
)

setMethod("total.concentration",
          "Stream",
          definition = function(object){
              concentration <- sum(fugacity(object))/(R_gas*(object@conditions[["temperature"]]+273))
              return(concentration)
          }
)

#standard diffusivity of components
setGeneric("std.diffusivity",
            def = function(object, ...){standardGeneric("std.diffusivity")}
)
setMethod("std.diffusivity",
          "Stream",
          definition = function(object, element = c("hydrogen", "nitrogen", "ammonia")){
              xH2 <- fraction(object,"hydrogen")
              xN2 <- fraction(object, "nitrogen")
              xNH3 <- fraction(object, "ammonia")
              to.call <- std.diff[element]
              diff <- sapply(to.call, function(x) 
                                        do.call(x, args = list(xH2, xN2, xNH3))
                             )
              return(diff)
          }
)
#diffusivity at stream conditions
setGeneric("diffusivity",
            def = function(object, ...){ststandardGeneric("diffusivity")}
)
setMethod("diffusivity",
          "Stream",
          definition = function(object, element = c("hydrogen", "nitrogen", "ammonia")){
              std.diff <- std.diffusivity(object, element)
              diff <- std.diff*((object@conditions[["temperature"]]+273)/273)^1.75*
                                (1/object@conditions[["pressure"]])
              return(diff)
          }
)
#effective diffusivity at stream conditions
setGeneric("effective.diffusivity",
           def = function(object, ...){standardGeneric("effective.diffusivity")}
)
setMethod("effective.diffusivity",
          "Stream",
          definition = function(object, 
                                element = c("hydrogen", "nitrogen", "ammonia"), phi = 0.52){
              eff.diff <- diffusivity(object, element)*phi/2
              return(eff.diff)
          }
)
#equilibrium constant at stream conditions
setGeneric(
    "equilibrium.constant",
    def = function(object){standardGeneric("equilibrium.constant")}
)
setMethod(
    "equilibrium.constant",
    "Stream",
    definition = function(object){
        eq <- pr.equilibrium.constant(object@conditions[["temperature"]])
        return(eq)
    }
)
#kinetic constant at stream conditions
setGeneric(
    "kinetic.constant",
    def = function(object, ...){standardGeneric("kinetic.constant")}
)
setMethod(
    "kinetic.constant",
    "Stream",
    definition = function(object, reaction){
        if ( class(reaction) != "reaction" )
            stop("Wrong class of reaction object")
        kc <- pr.kinetic.constant(reaction@properties[["preExp"]],
                                  reaction@properties[["activationEnergy"]], object@conditions[["temperature"]])
        return(kc)
    }
)
#heat of reaction
setGeneric(
    "heat.of.reaction",
    def = function(object, ...){standardGeneric("heat.of.reaction")}
)
setMethod(
    "heat.of.reaction",
    "Stream",
    definition = function(object){
        heat <- pr.heat.of.reaction(object@conditions[["temperature"]], object@conditions[["pressure"]])
        return(heat)
    }
)
#RNH3
setGeneric(
    "RNH3",
    def = function(object, ...){standardGeneric("RNH3")}
)
setMethod(
    "RNH3",
    "Stream",
    definition = function(object, reaction){
        if ( class(reaction) != "reaction" )
            stop("Wrong class of reaction object")
        kin.const <- kinetic.constant(object, reaction)
        eq.const <- equilibrium.constant(object)
        fgc <- fugacity(object)
        alpha <- reaction@properties[["alpha"]]
        rate <- pr.reaction.rate.NH3(kin.const, eq.const, fgc, alpha)
        return(rate)
    }
)

pr.std.diffusivity <- function(element = c(c("hydrogen", "nitrogen", "ammonia")),
                               xH2, xN2, xNH3){
    result <- sapply(element, function(x) std.diff[[x]](xH2, xN2, xNH3) )
    return(result)
}
    

pr.diffusivity <- function(std.diff, t, p){
    result <- sapply(std.diff, function(x) x*((t+273)/273)^1.75*
                                      (1/p)
                    )
    return(result)
}
pr.effective.diffusivity <- function(diff, theta){
    result <- sapply(diff, function(x) theta*x/2)
    return(result)
}

stream.mix <- function(stream1, stream2){
    new_stream <- stream1
    new_stream@mdot <- stream1@mdot+stream2@mdot
    fr1 <- flow.rate(stream1)
    fr2 <- flow.rate(stream2)
    cp1 <- heat.capacity(stream1)
    cp2 <- heat.capacity(stream2)
    new_stream@conditions[["temperature"]] <- as.double((fr1*cp1*stream1@conditions[["temperature"]]+
                                                        fr2*cp2*stream2@conditions[["temperature"]])/
                                                        (fr1*cp1+fr2*cp2))
    if (stream1@conditions[["pressure"]]!=stream2@conditions[["pressure"]])
        warning("stream preassure are not equal. pressure of stream1 is taken")
    return(new_stream)
}