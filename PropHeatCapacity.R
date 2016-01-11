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
    if ( p < 150 || p > 300)
        stop("Ammonia HC: Pressure is out of range 150-300 atm.")
    if ( t < 300 || t > 650)
        stop("Ammonia HC: Temperature is out of range 300-650 C")
    #convert pressure to multiple of 10
    p <- round(p/10, 0) * 10
    cf <- ammoniaHCcoef[[as.character(p)]]
    hc <- cf[["a"]]*(t ^ 4) +
          cf[["b"]]*(t ^ 3) +
          cf[["c"]]*(t ^ 2) +
          cf[["d"]]*t +
          cf[["e"]]
    return(hc)
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

ammoniaHCcoef <- list("150" = c(a = 8.8508E-10,
                              b = -1.9593E-6,
                              c = 0.0016397,
                              d = -0.60062,
                              e = 134.77),
                      "160" = c(a = 1.0034E-9,
                              b = -2.2152E-6,
                              c = 0.001849,
                              d = -0.67829,
                              e = 146.19),
                      "170" = c(a = 1.1266E-9,
                              b = -2.4814E-6,
                              c = 0.0020663,
                              d = -0.75872,
                              e = 157.96),
                      "180" = c(a = 1.2533E-9,
                              b = -2.755E-6,
                              c = 0.0022893,
                              d = -0.84111,
                              e = 169.97),
                      "190" = c(a = 1.382E-9,
                              b = -3.0326E-6,
                              c = 0.0025155,
                              d = -0.92456,
                              e = 182.11),
                      "200" = c(a = 1.5109E-9,
                              b = -3.3108E-6,
                              c = 0.002742,
                              d = -1.0081,
                              e = 194.26),
                      "210" = c(a = 1.6382E-9,
                              b = -3.5856E-6,
                              c = 0.002966,
                              d = -1.0908,
                              e = 206.28),
                      "220" = c(a = 1.7621E-9,
                              b = -3.8534E-6,
                              c = 0.0031845,
                              d = -1.1715,
                              e = 218.04),
                      "230" = c(a = 1.8808E-9,
                              b = -4.1102E-6,
                              c = 0.0033945,
                              d = -1.2494,
                              e = 229.42),
                      "240" = c(a = 1.9925E-9,
                              b = -4.3525E-6,
                              c = 0.0035932,
                              d = -1.3232,
                              e = 240.28),
                      "250" = c(a = 2.0955E-9,
                              b = -4.55767E-6,
                              c = 0.0037778,
                              d = -1.3923,
                              e = 250.52),
                      "260" = c(a = 2.1885E-9,
                              b = -4.78E-6,
                              c = 0.003946,
                              d = -1.4558,
                              e = 260.02),
                      "270" = c(a = 2.27E-9,
                              b = -4.9595E-6,
                              c = 0.004096,
                              d = -1.5128,
                              e = 268.69),
                      "280" = c(a = 2.3392E-9,
                              b = -5.1131E-6,
                              c = 0.0042258,
                              d = -1.563,
                              e = +276.46),
                      "290" = c(a = 2.3951E-9,
                              b = -5.2392E-6,
                              c = 0.0043341,
                              d = -1.6057,
                              e = 283.26),
                      "300" = c(a = 2.4374E-9,
                              b = -5.3367E-6,
                              c = 0.0044201,
                              d = -1.6407,
                              e = 289.06)
)

