reaction_propertynames <- c("preExp", "activationEnergy", "alpha")

element.match <- function(here, there){
    if (length(here)==length(there))
        result <- all(sapply(there, function(x) any(here == x)))
    else
        stop("Vectors' length does not equal!")
    return(result)
}

setClass(
    Class = "reaction",
    representation = representation(
        properties = "numeric"
    ),
    validity = function(object){
        if (is.null(names(object@properties)))
            stop("No names for properties")
        else if (!element.match(names(object@properties),reaction_propertynames))
            stop("Names in property vector are wrong")
        if(any(object@properties < 0 ))
            stop("Negative reaction parameters")
        return(TRUE)
    }
)

setMethod("initialize",
          "reaction",
          definition = function(.Object, properties){
              if (length(properties) != length(reaction_propertynames))
                stop("Wrong number of reaction properties")
              else if(is.null(names(properties))){
                  .Object@properties <- properties
                  names(.Object@properties) <- reaction_propertynames
              }
              else
                  .Object@properties <- properties
              validObject(.Object)
              return(.Object)
          }
)

