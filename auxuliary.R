element.match <- function(here, there){
    if (length(here)==length(there))
        result <- all(sapply(there, function(x) any(here == x)))
    else
        stop("Vectors' length does not equal!")
    return(result)
}