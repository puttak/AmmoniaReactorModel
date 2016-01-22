pr.equilibrium.constant <- function(t, p = NULL){
    #p = NULL for future use of eq const as function of pressure
    t <- t + T_abs
    log10Ka <-  - 2.691122 * log10(t) - 
                5.519265E-5 * (t) +
                1.848863E-7 * (t ^ 2) +
                2001.6 / t + 
                2.6899
    Kp <- 10^log10Ka
    return(Kp)
}


pr.kinetic.constant <- function(preExp, Ea, t){
    kc <- preExp*exp(-Ea/(8.314*(t + T_abs)))
    return(kc)
}

pr.heat.of.reaction <- function(t, p){
    t <- t + 273;
    if (p < 0 )
        stop("Negative pressure")
    if (t < 0)
        stop("Negative temperature")
    heat =  -(0.54526 + 840.609/t + 459.734E6/(t ^ 3))*p -
            5.34685*t - 0.2525E-3*(t ^ 2) +
            1.69167E-6*(t ^ 3) -
            9157.09
    return (heat * 4.183)
}

pr.reaction.rate.NH3 <- function(kin.const, eq.const, fugacity, alpha){
    #//[0] - Hydrogen, [1] - Nitrogen, [2] - Ammonia
    if ( is.null( names(fugacity) ) )
        stop("No names for fugacities")
    rate <- kin.const *
            (   eq.const^2 * fugacity[["nitrogen"]]*
                (fugacity[["hydrogen"]]^3 / fugacity[["ammonia"]] ^ 2) ^ alpha -
                (fugacity[["ammonia"]] ^ 2/fugacity[["hydrogen"]] ^ 3) ^ (1 - alpha)
             )
    return(rate)
}

# pure.RNH3 <- function(molefractions, fugacitycf, kin.const, eq.const, t, p, alpha){
#     rate <- kin.const *
#         (   eq.const^2 * (fugacitycf[["nitrogen"]]*molefractions[["nitrogen"]]*p)*
#                 ((fugacitycf[["hydrogen"]]*molefractions[["hydrogen"]]*p)^3 / 
#                      (fugacitycf[["ammonia"]]*molefractions[["ammonia"]]*p) ^ 2) ^ alpha -
#                 ((fugacitycf[["ammonia"]]*molefractions[["ammonia"]]*p) ^ 2/
#                      (fugacitycf[["hydrogen"]]*molefractions[["hydrogen"]]*p) ^ 3) ^ (1 - alpha)
#         )
#     return(rate)
# }