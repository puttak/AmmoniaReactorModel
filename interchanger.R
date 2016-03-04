#interchanger
#geom li a list of geometry for interchanger
interchanger <- function(d.sh.inner, d.sh.outer, dt.inner, dt.outer, n.tube,
                         l.tube, baff.spacing, pitch, lambda){
    arglist <- as.list(match.call())
    argcheck <- ArgumentCheck::newArgCheck()
    
    if (any(sapply(arglist, function(x) x < 0)))
        ArgumentCheck::addError(
            msg = "Negative values in arguments",
            argcheck = argcheck
        )
    if (d.sh.outer < d.sh.inner)
        ArgumentCheck::addError(
            msg = "Shell outer diameter less or equal to inner diameter",
            argcheck = argcheck
        ) 
    if (dt.outer < dt.inner)
        ArgumentCheck::addError(
            msg = "Tube outer diameter less or equal to inner diameter",
            argcheck = argcheck
        )
    #ensure that all geometric params are provided through geom
    # if (!element.match(names(geom), c("baffle.spacing","baffle.cut",  "nbaffle", "Nc", "pitch","nrows", "ds", "dsb")))
    #     ArgumentCheck::addError(
    #         msg = "Not all geometric parameters are provided in geom",
    #         argcheck = argcheck
    #     ) 
        
    #n tubes, pitch, number of tube rows, d shell inner, d tube outher
    #baffle spacing, buffle cut, 
    #bundle to cut clearance, shell to bundle clearance
    ArgumentCheck::finishArgCheck(argcheck)
    
    d.sh.inner <- d.sh.inner
    d.sh.outer <- d.sh.outer
    dt.inner <- dt.inner
    dt.outer <- dt.outer
    l.tube <- l.tube
    n.tube <- n.tube
    lambda <- lambda
    Acrosssec <- pi*(d.sh.inner+d.sh.outer)/2*(pitch-dt.outer)/pitch*baff.spacing
    
   
    #globals visible for all functions
    #assign within main function
    u.tube <- 0.0
    u.sh <-0.0
    pressure <- 0.0
    #ideal area of the shell
    #a.shell.ideal <- pi*((d.sh.outer^2-d.sh.iiner^2)/4-n.tube*dt.outer^2/4)
    #Returns heat capacity of stream
    stream.hc <- function(mdot,t,p){
        hc <- sapply(names(heat.capacity.list), function(x) heat.capacity(x, t, p))
        total <- sum(mdot)
        x_mass <- sapply(mdot, function(x) x/total)
        hc <- sum(hc*x_mass)
    }
    #returns mole fractions in the stream
    stream.mf <- function(mdot){
        mole <- sapply(names(mdot), function(x) mdot[[x]]/molWeight[[x]])
        total <- sum(mole)
        mole <- sapply(mole, function(x) x/total)
    }
    stream.viscosity <- function(mdot, t,p){
        sum(stream.mf(mdot)*viscosity(t,p))
    }
    stream.density <- function(mdot, t,p){
        total <- sum(mdot)
        rho <- sapply(mdot, function(x) x/total)
        rho <- sum(sapply(names(rho), function(x) rho[x]*molWeight[[x]]))*p / 
            (R_gas * (t + T_abs))
    }
    #stream thermal conductivity
    stream.k <-function(mdot, t,p){
        sum(stream.mf(mdot)*k.thermal(t,p))
    }
    #stream molecular weight kg/kmol
    stream.mw <- function(mdot){
        total <- sum(mdot)
        #    sum of     weight %                    * molWeight
        mw <- sum(sapply(mdot, function(x) x/total) * molWeight[names(mdot)])
    }
    #stream volumetric flow rate
    stream.vfr <- function(mdot, t, p){
        #Q = mdot_gas * R * T / (p * M_gas)
        sum(mdot) * R_gas * (t + T_abs) / (p * stream.mw(mdot))
    }
    #shell side heat transfer coefficient
    #n tubes, pitch, number of tube rows, d shell inner, d tube outher
    #baffle spacing, buffle cut, 
    #bundle to shell clearance, shell to baffle clearance

    de.func <- function(t.tube, t.sh, 
                        mdot.tube, mdot.sh){
        #densities
        rho.tube <- stream.density(mdot.tube, t.tube, pressure)
        rho.sh <- stream.density(mdot.sh, t.tube, pressure)
        #viscosities
        vis.tube <- stream.viscosity(mdot.tube, t.tube,pressure)
        vis.sh <- stream.viscosity(mdot.sh, t.tube,pressure)
        #thermal condictivities
        k.tube <- stream.k(mdot.tube,t.tube, pressure)
        k.sh <- stream.k(mdot.sh,t.sh,pressure)
        #heat capacities
        hc.tube <- stream.hc(mdot.tube,t.tube, pressure)
        hc.sh <- stream.hc(mdot.sh,t.sh, pressure)
        # gas vol flow rate /n tube * area of one tube
        u.tube <<- stream.vfr(mdot.tube,t.tube, pressure)/n.tube/(pi*dt.inner^2/4)
        #gas flow rate /average cross sectional area in radial direction
        u.sh <<- stream.vfr(mdot.sh,t.sh,pressure)/Acrosssec
                                                        #l.tube/(n.buffles-1))
        #HT from tube side
        alpha.tube <- htc(nusselt(reynolds(u.tube, dt.inner, rho.tube, vis.tube),
                                  prandtl(hc.tube, vis.tube,k.tube),
                                  vis.tube, 
                                  #viscosity(t.tube, element = "water"),
                                  type = "cooling"),
                          dpipe = dt.inner,
                          k.tube
        )
        #HT from shell side
        alpha.sh <- htc(nusselt(reynolds(u.sh, dt.outer, rho.sh, vis.sh),
                                prandtl(hc.sh, vis.sh,k.sh),
                                vis.sh, 
                                #viscosity(t.sh, element = "water"),
                                type = "heating"),
                        dpipe = dt.inner,
                        k.tube
        )
        #overall HT coefficient
        k.overall <- 1/(1/(alpha.tube*dt.inner/2)+1/lambda*log(dt.outer/dt.inner)+1/(alpha.sh*dt.outer/2))
        return(
            c(
              2*pi*n.tube*k.overall*(t.tube-t.sh)/(sum(mdot.tube)*hc.tube),
              2*pi*n.tube*k.overall*(t.tube-t.sh)/(sum(mdot.sh)*hc.sh)
              )
        )
    }
    function(mdot.tube, mdot.sh, t.tube, t.sh, p, n.step = 11){
        require(bvpSolve)
        if (n.step <= 0)
            stop("Negative number of n.step")
        else if (n.step < 11)
            stop("n.step of 11 is used")
        pressure <<- p
        #linear velocity in the tube
        # gas vol flow rate /n tube * area of one tube
#         u.tube <<- stream.vfr(mdot.tube,t.tube, pressure)/n.tube/(pi*dt.inner^2/4)
#         #gas flow rate /average cross sectional area in radial direction
#         u.sh <<- stream.vfr(mdot.sh,t.sh,pressure)/(pi*(d.sh.inner+d.sh.outer)/2*l.tube/(n.buffles-1))
        #ode function as a function of t tube and shell sides only
        odef <- function(l, t, params, ...) list(de.func(t[[1]], t[[2]], mdot.tube, mdot.sh))
        temps <- bvptwp(x = seq(0, l.tube, length.out = n.step), func = odef, yini = c(NA, t.sh), yend = c(t.tube, NA))
        return(temps)
        
    }
}

# shell.geom <- shell.geom(n, tube, geom[["pitch"]], geom[["Nc"]], d.sh.inner, dt.outer,
#                          geom[["baffle.spacing"]], geom[[""]])
#n tubes, pitch, number of tube rows, d shell inner, d tube outher
#baffle spacing, buffle cut, 
#bundle to shell clearance, shell to baffle clearance
shell.geom <- function(n.tubes, pitch, Nc, Ds, D0,
                      ls, lc,
                      db, dsb){
    #parallel pitch pitch * cos 30
    pp <- pitch*0.866025403
    # normal pitch pitch *sin30
    pn <- pitch*0.5
    Dotl <- Ds - db
    Fc <- 1/pi*
            (pi + 2*(Ds-2*lc)/Dotl*sin(acos((Ds-2*lc)/Dotl))-
                 2*acos((Ds-2*lc)/Dotl)
            )
    Ncw <- 0.8*lc/pp
    #cross flow area near centerline
    Sm <- ls*(Ds-Dotl +(Dotl+D0)/pitch*(pitch-D0))
    #fraction of cross flow available for bybass flow
    Fbr <- (Ds-Dotl)*ls/Sm
    #tube-to-baffle leakage area
    Stb <- 6.223E-4*D0*n.tubes*(1+Fc)
    #shell-to-baffle leakage area
    Ssb <- Ds*dsb/2*(pi-acos(1 - 2*lc/Ds))
    #area for flow through window
    Swg <- Ds^2/4*(acos(1-2*lc/Ds)-
                       (1-2*lc/Ds)*sqrt(1-(1-2*lc/Ds)^2))
    Swt <- n.tubes/8*(1-Fc)*pi*D0^2
    Sw <- Swg-Swt
    return(c(Sm=Sm, Ssb=Ssb, Stb=Stb, Sw=Sw))
}
S <- shell.geom(740, 0.03969, 41, 1.912, 0.03175+2*0.00165, 0.664, 1.692, 1.912-1.62729, 1.912-1.692)
Sx <- S["Ssb"]+S["Stb"]/S["Sm"]
Sy <- S["Ssb"]/(S["Ssb"]+S["Stb"])
#Jc = 0.7
