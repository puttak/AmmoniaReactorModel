setClass(
    Class = "Bed",
    representation = representation(
        inlet = "Stream",
        reaction = "reaction",
        catalyst = "Catalyst",
        volume = "numeric"
    ),
    validity = function(object){
        if (object@volume < 0)
            stop("Negative volume of catalyst bed")
        if (!validObject(object@inlet))
            stop("Invalid inlet stream")
        if (!validObject(object@reaction))
            stop("Invalid reaction")
        if (!validObject(object@catalyst))
            stop("Invalid catalyst")
        return(TRUE)
    }
)

setMethod("initialize",
          "Bed",
          definition = function(.Object, inlet, reaction, catalyst, volume){
              .Object@inlet <- inlet
              .Object@catalyst <- catalyst
              .Object@reaction <- reaction
              .Object@volume <- volume
              validObject(.Object)
              return(.Object)
          }
)