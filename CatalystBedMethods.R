#primitive functions for catalyst bed
#return mole of nitrogen if mole0 = mole0 and mole total = mole total
fractional.conversion.nitrogen <- function(mole0N2, moletotal, hi){
    if (mole0N2 <0 || moletotal <0 || hi < 0)
        stop("One or more parameters are negative")
    if (mole0N2 > moletotal)
        stop("Element mole is larger than total")
    return( mole0N2*(1 - hi)/(moletotal - 2*mole0N2*hi))
}

fractional.conversion.hydrogen <- function(mole0N2, mole0H2, moletotal, hi){
    if (mole0N2 < 0 || mole0H2 < 0 || moletotal <0 || hi < 0)
        stop("One or more parameters are negative")
    if (mole0N2 > moletotal)
        stop("Element mole is larger than total")
    return( (mole0H2 - 3*mole0N2*hi)/(moletotal - 2*mole0N2*hi) )
}

fractional.conversion.ammonia <- function(mole0N2, mole0H2, moletotal, hi){
    if (mole0N2 < 0 || mole0H2 < 0 || moletotal <0 || hi < 0)
        stop("One or more parameters are negative")
    if (mole0N2 > moletotal)
        stop("Element mole is larger than total")
    return( 2*mole0H2*hi/(moletotal - 2*mole0N2*hi) )
}

fractional.conversion.list <- list(nitrogen = fractional.conversion.nitrogen,
                                   hydrogen = fractional.conversion.hydrogen,
                                   ammonia = fractional.conversion.ammonia)

setGeneric("fractional.conversion",
           def = function(object, ...){standardGeneric("fractional.conversion")}
)
setMethod("fractional.conversion",
          "character",
          definition = function(object, ...){
              argstocall <- list(...)
              to.call <- fractional.conversion.list[[object]]
              result <- do.call(to.call, argstocall)
              return(result)
          }
    
)

setMethod("fractional.conversion",
          "Stream",
          definition = function(object, hi){
              moles <-  flow.rate(object, type = "mole")*fraction(object)
#               moles <- list(moletotal = totalmole, mol0N2 = moles[["nitrogen"]],
#                             mole0H2 = moles[["hydrogen"]])
              arguments <- list(nitrogen = list(object = "nitrogen", mole0N2 =  moles[["nitrogen"]], 
                                                moletotal = totalmole, hi = hi),
                                hydrogen = list(object = "hydrogen", mole0N2 = moles[["nitrogen"]], 
                                                mole0H2 = moles[["hydrogen"]], moletotal = totalmole, hi = hi),
                                ammonia = list(object = "ammonia", mole0N2 = moles[["nitrogen"]], 
                                               mole0H2 = moles[["hydrogen"]], moletotal = totalmole, hi = hi)                                
                                )
              fc <- sapply(c("nitrogen", "hydrogen", "ammonia"),
                           function(x) do.call(fractional.conversion, arguments[[x]]) 
                           )
              return(fc)
          }
)

#methods to recalculate masses of components in a flow based on conversion
hi.mass.nitrogen <- function(mass0N2, hi){
    res <- mass0N2*(1-hi)
    return(res)
}
hi.mass.hydrogen <- function(mole0N2, mole0H2, hi){
    res <- (mole0H2 - 3*mole0N2*hi)*molWeight[["hydrogen"]]/1000
    return(res)
}
hi.mass.ammonia <- function(mole0N2, mole0NH3, hi){
    res <- (mole0NH3 + 2*mole0N2*hi)*molWeight[["ammonia"]]/1000
    return(res)
}
hi.mass.list <- list(nitrogen = hi.mass.nitrogen,
                     hydrogen = hi.mass.hydrogen,
                     ammonia = hi.mass.ammonia
                     )
setGeneric("recalculate.stream",
           def = function(object, ...){standardGeneric("recalculate.stream")}
)
setMethod("recalculate.stream",
          "Stream",
          definition = function(object, hi, return.stream = FALSE){
              #prepares arguments to be passed for each of hi.mass functions
              arguments <-list( nitrogen = list(mass0N2 = flow.rate(object, type = "mass")*
                                                    fraction(object, element = "nitrogen",type = "mass"),
                                                 hi = hi),
                                hydrogen = list(mole0N2 = flow.rate(object, type = "mole")*
                                                    fraction(object, element = "nitrogen", type = "mole"),
                                                mole0H2 = flow.rate(object, type = "mole")*
                                                    fraction(object, element = "hydrogen", type = "mole"),
                                                hi = hi),
                                ammonia = list(mole0N2 = flow.rate(object, type = "mole")*
                                                    fraction(object, element = "nitrogen", type = "mole"),
                                               mole0NH3 = flow.rate(object, type = "mole")*
                                                    fraction(object, element = "ammonia", type = "mole"),
                                               hi = hi)
                                )
              #recalculates mass flows of N2, H2, NH3 for a given conversion hi
              result <- sapply(c("nitrogen", "hydrogen", "ammonia"),
                                 function(x) do.call(hi.mass.list[[x]], args = arguments[[x]] )
                                 )
              #add mass rates of methane and argon to the mass vector
              if (return.stream){
                  new_stream <- c(result, object@mdot[c("methane", "argon")])
                  result <- object
                  result@mdot <- new_stream
                  result@conditions[["temperature"]] <- NA
                  result@conditions[["pressure"]] <- NA
              }
              return(result)
          }
)
# returns bed as a function of its inlet
bed.db <- function(bed){
    volume <- bed@volume
    reaction <- bed@reaction
    catalyst <- bed@catalyst
    stream <- bed@inlet
    inlet_molar_fl_N2 <- flow.rate(bed@inlet,type = "mole")*fraction(bed@inlet, "nitrogen")
    deb <- function(stream){
        #effectiveness factor in stream
        ef <- effectiveness.factor(stream, reaction, catalyst)
        rate <- RNH3(stream, reaction)
        mb <- ef * rate / (2 * 3.6 * inlet_molar_fl_N2)
        eb <- - heat.of.reaction(stream) * ef * rate / (flow.rate(stream, type = "mole")* 3.6 *
                                                            heat.capacity(stream, hc.un = "mole"))
        result <- c(mb, eb)
        #patch for now to fix names from ef function
        names(result) <- NULL
        names(result) <- c("mass", "heat")
        return(result)
    }
}

bed.ode.func <- function(bed.db){
    inlet <- environment(bed.db)$stream
    func <- function(t, y, parms = NULL){
        stream <- recalculate.stream(inlet, y[[1]], return.stream = TRUE)
        stream@conditions <- c(temperature=y[[2]], pressure=inlet@conditions[["pressure"]])
        return(list(bed.db(stream)))
    }
}

bed.calculate <- function(bed.ode.func, x0, method = "ode45", n.iter = 10){
    require(deSolve)
    vol <- environment(environment(bedode)$bed.db)$bed@volume
    t_in <- environment(environment(bedode)$bed.db)$bed@inlet@conditions[["temperature"]]
    ode(c(x0, t_in), seq(0, vol, length.out = n.iter), bed.ode.func, NULL, method = method)
}