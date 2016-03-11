#ensloses 3 beds within converter and heat exchanger
converter <- function(bed1, bed2, bed3, interchanger){
    #bed1 encloses inlet stream, bed2 and bed3 have NULL inlet
    beds <- list("1"=bed1, "2"=bed2, "3"=bed3)
    bedodes <- vector("list", 3)
    names(bedodes) <- c("1", "2", "3")
    interchanger <- interchanger
    #variables to store internal results for beds 1 and 2
    bedresult <- vector("list",3)
    names(bedresult) <- c("1", "2", "3")
    #heat_exchanher temperature profile
    HE_temperature <- NULL
    #list of quenches to first and second bed
    quench <- NULL
    unit_pressure <- NULL
    feed_temperature <- NULL
    third_bed_inlet <- NULL
    #ode function which returns concentration os function of temperature
    #environment(environment(bedode)$bed.db)$volume
    zerof <- function(temperature_in){
        #set temperature of bed1 inlet to temperature_in
        environment(bedodes[[1]])$inlet@conditions[["temperature"]] <- temperature_in
        ###delete:environment(environment(bedodes[[1]])$bed.db)$stream@conditions[["temperature"]] <- temperature_in
        #if quench for 1st bed is provided, mix inlet and quench
        if (!is.null(quench[[1]]))
            environment(bedodes[[1]])$inlet <- 
                stream.mix(environment(bedodes[[1]])$inlet, quench[[1]])
            ###delete: environment(environment(bedodes[[1]])$bed.db)$stream <- 
            ###delete:    stream.mix(environment(environment(bedodes[[1]])$bed.db)$inlet, quench[[1]])
        #calculate 1st bed
        bedresult[[1]] <<- bed.calculate(bedodes[[1]])
        #get stream based on conversio from previous bed
        ##create an empty stream (in a lame way)
        #get outlet stream from 1st bed
        temp_stream <-recalculate.stream(
            environment(bedodes[[1]])$inlet, 
            as.double(bedresult[[1]][nrow(bedresult[[1]]), "1"]), TRUE
        )
        temp_stream@conditions[["temperature"]] <-  as.double(bedresult[[1]][nrow(bedresult[[1]]), "2"])
        temp_stream@conditions[["pressure"]] <- unit_pressure
        if (is.null(bedodes[[2]]))
            bedodes[[2]] <- bed.ode.func(bed.db(beds[[2]]),temp_stream)
        else
            environment(bedodes[[2]])$inlet <- temp_stream
        rm(temp_stream)
    # environment(environment(bedodes[[2]])$bed.db)$stream <- environment(environment(bed1ode)$bed.db)$stream
    # environment(environment(bedodes[[2]])$bed.db)$stream@mdot <- recalculate.stream(
    #     environment(environment(bed1ode)$bed.db)$stream, as.double(bedresult[[1]][nrow(bedresult[[1]]), "1"], TRUE)
    #     )
    # environment(environment(bed2ode)$bed.db)$stream@conditions[["temperature"]] <-
    #     as.double(bedresult[[1]][nrow(bedresult[[1]]), "2"])
        #mix 1st bed outlet and 2nd quench
        if (!is.null(quench[[2]])){
            environment(bedodes[[2]])$inlet <- 
                stream.mix(environment(bedodes[[2]])$inlet, quench[[2]])
        } else
            message("No quench for 2nd bed given")
        bedresult[[2]] <<- bed.calculate(bedodes[[2]])
        temp_stream <- recalculate.stream(
            environment(bedodes[[1]])$inlet, as.double(bedresult[[1]][nrow(bedresult[[1]]), "1"]), TRUE
        )
        temp_stream@conditions[["temperature"]] <-  as.double(bedresult[[2]][nrow(bedresult[[2]]), "2"])
        temp_stream@conditions[["pressure"]] <- unit_pressure
        HE_temperature <<- interchanger(temp_stream@mdot,
                                        environment(bedodes[[1]])$inlet@mdot,
                                        temp_stream@conditions[["temperature"]],
                                        feed_temperature,
                                        unit_pressure)
        third_bed_inlet <<- temp_stream
        rm(temp_stream)
        return(temperature_in - as.double(HE_temperature[nrow(HE_temperature), "2"]))
    }
    
    return(function(inlet, quench){
        bedresult <<- NULL
        bedodes[[1]] <<- bed.ode.func(bed.db(beds[[1]]), inlet)
        unit_pressure <<- as.double(inlet@conditions[["pressure"]])
        feed_temperature <<- as.double(inlet@conditions[["temperature"]])
        quench <<- quench
        #solve for first two beds and interchanger
        result <- uniroot(zerof, c(230, 500))$root
        #get stream entering into 3rd bed
        
        #create bed db for 3rd bed
        bedodes[[3]] <<- bed.ode.func(bed.db(beds[[3]]), third_bed_inlet)
        bedresult[[3]] <<- bed.calculate(bedodes[[3]])
        return(bedresult)
    })
}