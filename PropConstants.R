elementNames <- c("nitrogen",
                  "hydrogen",
                  "ammonia",
                  "methane",
                  "argon")
molWeight <- c(nitrogen = 28.0,
               hydrogen = 2.0,
               ammonia = 17.0,
               methane = 16.0,
               argon = 40.0)
diffusivity <- c(nitrogen = 0.00050128,
                 hydrogen = 0.0001305,
                 ammonia = 0.00037187,
                 methane = 0.00030204,
                 argon = 0.00068376)
std.diff <- list(hydrogen = function(xH2, xN2, xNH3){
                                sum <- xN2/0.571e-4+xNH3/0.629e-4
                                sum <- (1 - xH2)/sum
                                return(sum*3600)
                            },
                 nitrogen = function(xH2, xN2, xNH3){
                                sum <- xH2/0.571e-4+xNH3/0.161e-4
                                sum <- (1 - xN2)/sum
                                return(sum*3600)    
                            },
                 ammonia = function(xH2, xN2, xNH3){
                                sum <- xH2/0.629e-4+xN2/0.161e-4
                                sum <- (1 - xNH3)/sum
                                return(sum*3600)
                            }
                )

R_gas <- 8.205e-2
T_abs <- 273.13