#return fugacity coefficient of element
setGeneric("fugacity.cf",
           def = function(object, ...){standardGeneric("fugacity.cf")}
)

setMethod("fugacity.cf",
          "character",
          definition =  function(object, t = 15, p = 1, un = c("C", "atm")){
                if (un[1] == "C")
                    temperatureK <- t + 273.15
                else if ( un[1] == "K" )
                    temperatureK <- t
                if (un[2] == "atm")
                    pressure <- p
                fgc.to.call <- fugacity.list[[object]]
                fgc.to.call(temperatureK,pressure)
          }
)

fugacity.nitrogen <- function(t, p){
    fgc <- 0.93431737 + 
            0.3101804E-3 * t +
            0.295896E-3 * p - 
            0.2707279E-6 * t ^ 2+
            0.4775207E-6 * p ^ 2
    return(fgc)
}

fugacity.hydrogen <- function(t, p){
    exp1 <- exp(-3.8402 * (t ^ 0.125) + 0.541)
    exp2 <- exp(-0.1263 * (t ^ 0.5) - 15.980)
    exp3 <- exp(-0.011901*t - 5.941)
    fgc <- exp( exp1 * p
                - exp2 * (p ^ 2)
                + 300 * exp3 * (exp(-p / 300) - 1));
    
    return(fgc)
}

fugacity.ammonia <- function(t, p){
    fgc <- 0.1438996 +
            0.2028538E-2 * t -
            0.4487672E-3 * p -
            0.1142945E-5 * (t ^ 2) +
            0.2761216E-6 * (p ^ 2)
    return(fgc)
}

fugacity.methane <- function(t, p){ return(1)}
fugacity.argon <- function (t, p) {return(1)}

fugacity.list <- list(nitrogen = fugacity.nitrogen,
                      hydrogen = fugacity.hydrogen,
                      ammonia = fugacity.ammonia,
                      methane = fugacity.methane,
                      argon = fugacity.argon)