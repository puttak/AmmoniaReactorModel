load("ammonia.RData")
reynolds <- function(v, d, rho, u){
    return(v*d*rho/u)
}

prandtl <- function(cp, u, k){
    return(cp*u/k)    
}

nusselt <- function(Re, Pr, ub, uw=0.89040E-3, type = "heating"){
    if (Re < 10000)
        message(paste("Reynolds number in nusselt function is outside the range, Re=",
                      Re)
        )
    if ((Pr<0.7)||(Pr>170))
        message(paste("Prandtl number in nusselt function is outside the range, Pr=",
                      Pr)
        )
    if (type == "heating")
        return(0.0243*Re^0.8*Pr^0.4*(ub/ub)^0.14)
    if (type == "cooling")
        return(0.0265*Re^0.8*Pr^0.3*(ub/ub)^0.14)
}
#heat transfert coefficient
htc <- function(Nu, dpipe, k){
    return(Nu*k/dpipe)
}
#return elements viscosity by name at t and p
viscosity <- function(t,p=NULL,element = c("hydrogen", "nitrogen", "ammonia","methane", "argon")){
    return(sapply(viscosity.list[element],function(x) do.call(x, args = list(t=t, p=p))))
}

pr.viscosity <- function(C,T_0, u_0){
    C <- C
    T_0 <- T_0
    u_0 <- u_0
    function(t,p){
        t <- (t+273.15)*9/5
        return(0.0001 * u_0 *
                   (0.555*T_0+C)/(0.555*t + C) *
                   (t/T_0)^1.5
        )
    }
}

# viscosity.list <- list(
#     hydrogen = pr.viscosity(72, 528.93, 0.00876),
#     nitrogen = pr.viscosity(111, 540.99, 0.01781),
#     ammonia = pr.viscosity(370, 527.67, 0.00982),
#     methane = pr.viscosity(370, 527.67, 0.00982),
#     argon = pr.viscosity(370, 527.67, 0.00982),
#     water = pr.viscosity(5000, 527.67, 0.2)
# )
viscosity.list <- list(
    hydrogen = function(t,p){return(0.001*(9.172E-3+1.093E-8*p*101.325+1.832E-5*t))},
    nitrogen = function(t,p){return(0.001*(3.079E-2+2.279E-7*p*101.325+ 4.383E-5*t))},
    ammonia = function(t,p){
        require(mgcv)
        return(0.001*as.double(predict(ammprop$am.viscosity, list(t=t, p=p*101.325))))
        },
    methane = function(t,p){return(0.001*(1.402E-2+1.155E-7*p*101.325+1.657E-5*t))},
    argon = function(t,p){return(0.001*(2.276E-2+1.839E-7*p*101.325+2.432E-5*t))}
    
)
#create viscosity.element function and list them into viscosity.list
#W/mk
k.thermal.list<- list(
    hydrogen = function(t,p){return(1.932E-1+3.219E-7*p*101.325+3.089E-4*t)},
    nitrogen = function(t,p){return(1.877E-2+1.051E-7*p*101.325+3.616E-5*t)},
    ammonia = function(t,p){
        require(mgcv)
        return(0.001*as.double(predict(ammprop$am.k, list(t=t,p=p*101.325))))
        },
    methane = function(t,p){return(2.850E-2+4.014E-7*p*101.25+1.773E-4*t)},
    argon = function(t,p){return(2.516E-2+1.339E-7*p*101.325+4.160E-5*t)}
)

#70E-3
#25E-3

k.thermal <- function(t,p=NULL,element = c("hydrogen", "nitrogen", "ammonia", "methane", "argon")){
    return(sapply(k.thermal.list[element],function(x) do.call(x, args = list(t=t, p=p))))
}