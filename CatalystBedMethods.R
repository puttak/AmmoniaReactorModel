#primitive functions for catalyst bed
#return mole of nitrogen if mole0 = mole0 and mole total = mole total
fractional.conversion.nitrogen <- function(mole0N2, moletotal, hi){
    if (mole0N2 <0 || moletotal <0 || hi < 0)
        stop("One or more parameters are negative")
    if (mole0N2 > moletotal)
        stop("Element mole is latger than total")
    return( mole0N2*(1 - hi)/(moletotal - 2*mole0N2*hi))
}

fractional.conversion.hydrogen <- function(mole0N2, mole0H2, moletotal, hi){
    if (mole0N2 < 0 || mole0H2 < 0 || moletotal <0 || hi < 0)
        stop("One or more parameters are negative")
    if (mole0N2 > moletotal)
        stop("Element mole is latger than total")
    return( (mole0H2 - 3*mole0N2*hi)/(moletotal - 2*mole0N2*hi) )
}

fractional.conversion.ammonia <- function(mole0N2, mole0H2, moletotal, hi){
    if (mole0N2 < 0 || mole0H2 < 0 || moletotal <0 || hi < 0)
        stop("One or more parameters are negative")
    if (mole0N2 > moletotal)
        stop("Element mole is latger than total")
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
              totalmole <- flow.rate(object, type = "mole") 
              moles <-  totalmole*fraction(object)
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
          definition = function(object, hi){
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
              new_mass <- sapply(c("nitrogen", "hydrogen", "ammonia"),
                                 function(x) do.call(hi.mass.list[[x]], args = arguments[[x]] )
                                 )
              return(new_mass)
          }
)
