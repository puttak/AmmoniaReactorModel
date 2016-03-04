#ensloses 3 beds within converter and heat exchanger
converter <- function(bed1, bed2, bed3, interchanger){
    #bed1 encloses inlet stream, bed2 and bed3 have NULL inlet
    bed1ode <- bed.ode.func(bed.db(bed1))
    bed2ode <- bed.ode.func(bed.db(bed2))
    bed3ode <- bed.ode.func(bed.db(bed3))
    interchanger <- interchanger
    #variables to store internal results for beds 1 and 2
    bedresult <- vector("list",3)
    names(bedresult) <- c("1", "2", "3")
    #list of quenches to first and second bed
    quench <- NULL 
    #ode function which returns concentration os function of temperature
    #environment(environment(bedode)$bed.db)$volume
    zerof <- function(temperature_in){
        #set temperature of bed1 inlet to temperature_in
        environment(environment(bed1ode)$bed.db)$stream@conditions[["temperature"]] <- temperature_in
        #if quench for 1st bed is provided, mix inlet and quench 
        if (!is.null(quench[[1]]))
            environment(environment(bed1ode)$bed.db)$stream <- 
                stream.mix(environment(environment(bed1ode)$bed.db)$inlet, quench[[1]])
        #calculate 1st bed
        bedresult[[1]] <<- bed.calculate(bed1ode,0)
        #get stream based on conversio from previous bed
        ##create an empty stream (in a lame way)
        environment(environment(bed2ode)$bed.db)$stream <- environment(environment(bed1ode)$bed.db)$stream
        environment(environment(bed2ode)$bed.db)$stream@mdot <- recalculate.stream(
            environment(environment(bed1ode)$bed.db)$stream, as.double(bedresult[[1]][nrow(bedresult[[1]]), "1"], TRUE)
            )
        environment(environment(bed2ode)$bed.db)$stream@conditions[["temperature"]] <- 
            as.double(bedresult[[1]][nrow(bedresult[[1]]), "2"])
        #mix 1st bed outlet and 2 nd quench
        if (!is.null(quench[[2]])){
            environment(environment(bed2ode)$bed.db)$stream <- 
                stream.mix(environment(environment(bed2ode)$bed.db)$inlet, quench[[2]])
        } else
            message("No quench for 2nd bed given")
        bedresult[[2]] <<- bed.calculate(bed2ode, 0)
        return(temperature_in - as.double(bedresult[[2]][nrow(bedresult[[2]]), "2"]))
    }
    
    return(function(inlet, quench){
        environment(environment(bed1ode)$bed.db)$stream <<- inlet
        quench <<- quench
        zerof(350)
    })
}