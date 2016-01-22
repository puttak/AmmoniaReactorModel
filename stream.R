setClass(
    #Class stream
    #Contains stream composition, flow rate and conditions (t, p)
    Class = "Stream",
    representation = representation(
        mdot = "numeric",
        conditions = "numeric"
    ),
    validity = function(object){
        if ( any(object@mdot < 0))
            stop("Negative mass rate occured in stream")
        if ( !element.match(names(object@mdot), elementNames) )
            stop("Wrong names in mdot of stream")
        if ( any(object@conditions < 0) )
            stop("Negative stream conditions")
        if ( !element.match(names(object@conditions), c("temperature", "pressure")) )
            stop("Wrong names in conditions of stream")
        return(TRUE)
    }
)
#Stream constructor
setMethod("initialize",
          "Stream",
          definition = function(.Object, t, p, flowrate, composition){
              #flowrate = [m3/h]
              .Object@conditions <- c(t,p)
              names(.Object@conditions) <- c("temperature", "pressure")
              #check if composition length is not 5
              if ( length(composition) != 5 )
                  stop("Wrong length of concentration vector")
              #check that sum of concentrations is = 1
              else if(abs(sum(composition))-1 > 1e-3 )
                  stop("Check concentrations of feed")
              #check if names are correct or no names then convert to vol rate
              else if(all(sapply(elementNames, function(x) any(names(composition) == x)))){
                  reorder <- sapply(elementNames, function(x) which( (x==names(composition))==TRUE ) )
                  mdot <- composition[reorder]
                  mdot <- flowrate*mdot/3600
              }
              else if (is.null(names(composition))){
                  mdot <- flowrate*composition/3600
                  names(composition) <- elementNames
              }
#               else
#                   stop("Something wrong with concentrations")
              
              # this.FlowStructure[i].Mass = volRate[i]*this.FlowStructure[i].MolecularWeight/(8.205e-2*273.13);
              mdot <- mdot*molWeight/(R_gas*T_abs)
              .Object@mdot <- mdot
              validObject(.Object)
              return(.Object)
          }
          )