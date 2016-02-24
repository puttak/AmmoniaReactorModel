load("ammoniaHCcoefs.RData")
#generic heat capacity
setGeneric(
    name = "heat.capacity",
    def = function(object, ...){standardGeneric("heat.capacity")}
)
# ,
# def = function(...){standardGeneric("heat.capacity")}
setMethod("heat.capacity",
          "character",
          definition = function(object, t = 15, p = 1, hc.un = "mass", un = c("C", "atm")){
              #check t units and convert to C
              if (un[1] == "C")
                  temperatureC <- t
              else if (un[1] == "K")
                  temperatureC <- t - 273
              else 
                  stop("Temperature units are not supported! Only C or K!")
              #check t values
              if (temperatureC <= -273)
                  stop("Temperature below absolute zero!")
              #check p unit and convert to atm.
              if (un[2] == "atm")
                  pressureAtm <- p
              else if (un[1 == "kPa"])
                  pressureAtm <- p/101.325
              else
                  stop("Pressure units are not supported! Only atm. or kPa!")
              if (pressureAtm <= 0)
                  stop("Pressure below zero!")
              hc.to.call <- heat.capacity.list[[object]]
              #heat capacity in J/moleK
              hc <- hc.to.call(temperatureC, pressureAtm)
              #if specific heat is required
              if (hc.un == "mass")
                  hc <- hc / molWeight[[object]] * 1000
              return(hc)
              
          }
)
# heat.capacity <- function(element, t = 15, p = 1, hc.un = "mass", un = c("C", "atm")){
#     #check t units and convert to C
#     if (un[1] == "C")
#         temperatureC <- t
#     else if (un[1] == "K")
#         temperatureC <- t - 273
#     else 
#         stop("Temperature units are not supported! Only C or K!")
#     #check t values
#     if (temperatureC <= -273)
#         stop("Temperature below absolute zero!")
#     #check p unit and convert to atm.
#     if (un[2] == "atm")
#         pressureAtm <- p
#     else if (un[1 == "kPa"])
#         pressureAtm <- p/101.325
#     else
#         stop("Pressure units are not supported! Only atm. or kPa!")
#     if (pressureAtm <= 0)
#         stop("Pressure below zero!")
#     hc.to.call <- heat.capacity.list[[element]]
#     #heat capacity in J/moleK
#     hc <- hc.to.call(temperatureC, pressureAtm)
#     #if specific heat is required
#     if (hc.un == "mass")
#         hc <- hc / molWeight[[element]] * 1000
#     return(hc)
#     
# }

#heat capacities at temperature, pressure
#t = [C], p = [atm]
heat.capacity.nitrogen <- function(t, p){
    temperature <- t + 273
    #heat capacity in J/mole K
    hc <-   6.903 - 
            0.03753E-2 * temperature + 
            0.1930E-5 * (temperature ^ 2) - 
            0.6861E-9 * (temperature ^ 3)
    return(hc*4.183)
}

heat.capacity.hydrogen <- function(t, p){
    temperature <- t + 273
    hc <-   6.952 - 
            0.04576E-2 * temperature + 
            0.09563E-5 * (temperature ^ 2) -
            0.2079E-9*(temperature ^ 3)
    return(hc*4.183)
}

heat.capacity.ammonia <- function(t, p){
#     if ( p < 50 || p > 260)
#         stop("Ammonia HC: Pressure is out of range 50-250 atm.")
#     if ( t < 230 || t > 600)
#         stop("Ammonia HC: Temperature is out of range 230-600 C")
    #convert pressure to multiple of 10
    p <- round(p/10, 0) * 10
    cf <- ammoniaHCcoefs[,as.character(p)]
    hc <- cf[["t4"]]*(t ^ 4) +
          cf[["t3"]]*(t ^ 3) +
          cf[["t2"]]*(t ^ 2) +
          cf[["t1"]]*t +
          cf[["t0"]]
    return(hc*molWeight[["ammonia"]])
}

heat.capacity.methane <- function(t, p){
    temperature <- t + 273
    hc <-   4.750 +
            1.200E-2 * temperature +
            0.3030E-5 * (temperature ^ 2) -
            2.630E-9 * (temperature ^ 3)
    return(hc*4.183)
}

heat.capacity.argon <- function(t=15, p = 1){
    hc <- 4.9675
    return(hc*4.183)
}

heat.capacity.list <- list(nitrogen = heat.capacity.nitrogen,
                           hydrogen = heat.capacity.hydrogen,
                           ammonia = heat.capacity.ammonia,
                           methane = heat.capacity.methane,
                           argon = heat.capacity.argon)
