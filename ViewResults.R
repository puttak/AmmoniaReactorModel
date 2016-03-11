bed.summary <- function(deSolveMatrix, bed, stream){
    pressure <- stream@conditions[["pressure"]]
    get.streams <- function(x){
        stream <- recalculate.stream(stream, x[[2]], return.stream = TRUE)
        stream@conditions <- c(temperature=x[[3]], pressure=pressure)
        return(stream)
    }
    fractions <- apply(deSolveMatrix,1, get.streams)
    ef <- sapply(fractions, function(x) effectiveness.factor(x, bed@reaction, bed@catalyst))
    fractions <- sapply(fractions, fraction)
    fractions <- t(rbind(fractions, ef))
    return(fractions)
}