conc <- c(hydrogen=0.2219,
          nitrogen=0.6703,
          ammonia=0.0276,
          methane=0.0546, 
          argon=0.0256)
conc2 <- c(0.6703,0.2219,0.0276,0.0546,0.0256)
inlet1 <- new("Stream", 385, 226, 242160, conc)
inlet2 <- new("Stream", 385, 226, 242160, conc2)

reaction1 <- new("reaction", c(preExp = 1e7, activationEnergy = 48500, alpha = 0.52))

bed1 <- new("Bed", inlet1, reaction1, 4.75)
