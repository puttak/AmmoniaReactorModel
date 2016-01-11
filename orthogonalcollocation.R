ocfunction <- function( xH2, xN2, xNH3, 
                        t, p, C_total,
                        preExp, Ea, alpha,
                        intrPorosity, voidage, rParticle){
    kin.const <- pr.kinetic.constant(preExp, Ea, t)
    xMatr <- matrix(data = c(xH2, xN2, xNH3), ncol = 3, dimnames = list(NULL,c("hydrogen",
                                                                               "nitrogen",
                                                                               "ammonia")
                                                                        )
                    )
    fgc.cf <- sapply(fugacity.list[c("hydrogen", "nitrogen","ammonia")],
                     function(x) do.call(x, args = list(t=t, p=p))
                    )
    fgc <- xMatr%*%diag(fgc.cf)
    #second derivative
    #d2ydx2 <- apply(B %*% xMatr, 1, sum)
    d2ydx2 <- B %*% c(xMatr)
    #first derivative
    #dydx <- apply(A %*% xMatr, 1, sum)
    dydx <- A %*% c(xMatr)
    #2rd term in equation
    term2 <- -2 * (gamma_i + c(xMatr)) * dydx ^2
    #3rd term in equation
    term3 <- 2 * (1 / w_1) * dydx
    #4th term
    
    
}
debug(ocfunction)
ocfunction(seq(0.7, fraction(inlet1, "hydrogen"), length.out = 4), 
           seq(0.4, fraction(inlet1, "nitrogen"), length.out = 4), 
           seq(0.1, fraction(inlet1, "ammonia"), length.out = 4),
           385, 226, total.concentration(inlet1),
           reaction1@properties[["preExp"]], reaction1@properties[["activationEnergy"]],
           reaction1@properties[["alpha"]],
           catalyst1@geometry[["intrPorosity"]],catalyst1@geometry[["porosity"]],
           catalyst1@geometry[["pelletD"]]/2)

A <- matrix(data = c(-4.1308947, 6.8819128, -4.5475385, 1.7965206,
                             -1.3388630, -2.2150478, 5.2890548, -1.7351440,
                             0.62570332, -3.7406161, -1.6671150, 4.7820278,
                             -1.0727282, 5.3255695, -20.752841, 16.5),
                    nrow = 4,
                    byrow = TRUE
)

B <- matrix(data = c(-23.853065, 30.593651, -9.74629544, 3.0057072,
                             11.099906, -43.237662, 40.818768, -8.6810122,
                             -3.3228457, 38.356814, -125.40927, 90.375305,
                             -33.675598, 152.37521, -311.19961, 192.5),
                    nrow = 4,
                    byrow = TRUE
                    
)

A <- A[-4,]
B <- B[-4,]
B <- Matrix::bdiag(B,B,B)
A <- Matrix::bdiag(A,A,A)

w_i <- c(0.36311746, 0.67718628, 0.89975799)
gamma_i <- c(rep(-1.5,3), rep(-0.5,3), rep(1,3))



ocfunction <- function(x,stream,reaction,catalyst){
    #require("matrixcalc")
    require("Matrix")
#     if (!is.square.matrix(A) || !is.square.matrix(B))
#         stop("Non-square collocatiopn matrix")
#     if (dim(A)!= dim(B))
#         stop("Dimensions of A and B does not match")
    A <- ocAmatrix[-4,]
    B <- ocBmatrix[-4,]
    number_eq <- dim(A)[1]
    number_var <- dim(A)[2]
#     if (length(x) != length(y))
#         stop("Dimensions of x and y does not match")
#     if ( floor(x/3) != 0 ||floor(number_eq/3) != 0)
#         stop("Number of equations is not a multiple of 3")
    #vector of fugacities
    fgccf_i <- fugacity.cf(stream)[1:3]
    #vector of collocation points
    w_i <- oc.points[1:3]
    #reaction rate in collocation points at bulk conditions
    #and concentrations xi in collocation points
    equilibrium_constant <- pr.equilibrium.constant(stream@conditions[["temperature"]])
    kinetic_constant <- pr.kinetic.constant(preExp = reaction@properties[["preExp"]],
                                            Ea = reaction@properties[["activationEnergy"]],
                                            t = stream@conditions[["temperature"]])
    C_total <- total.concentration(stream)
    D_ie <- effective.diffusivity(stream, phi = catalyst@geometry[["intrPorosity"]])
    r_particle <- catalyst@geometry[["pelletD"]]
    voidage <- catalyst@geometry[["porosity"]]
    R_i <- sapply(1:3, function(i){
                            pure.RNH3(c(hydrogen = x[[i]], nitrogen = x[[i+3]], ammonia = x[[i+6]]),
                                        fgccf_i,
                                        kinetic_constant,
                                        equilibrium_constant,
                                        stream@conditions[["temperature"]],
                                        stream@conditions[["pressure"]],
                                        reaction@properties[["alpha"]]
                                       )    
                    }
                  )
    #vector of fugacities
    expand_fgccf_i <- as.vector(sapply(fgccf_i, function(x) rep(x,3)))
    expand_R_i <- rep(R_i, 3)
    expand_D_ie <- as.vector(sapply(D_ie, function(x) rep(x,3)))
    #expand matrix
    B <- Matrix::bdiag(B,B,B)
    A <- Matrix::bdiag(A,A,A)
    #second derivative
    d2ydx2 <- B %*% x
    #first derivative
    dydx <- A %*% x
    result <- d2ydx2 - dydx^2/(expand_gamma_i+x)+
                2/oc.points[1:3]*dydx +
                r_particle^2/(C_total + expand_D_ie)*(expand_gamma_i+x)*expand_R_i/(1-voidage)
    return(result)
}

ocAmatrix <- matrix(data = c(-4.1308947, 6.8819128, -4.5475385, 1.7965206,
                             -1.3388630, -2.2150478, 5.2890548, -1.7351440,
                             0.62570332, -3.7406161, -1.6671150, 4.7820278,
                             -1.0727282, 5.3255695, -20.752841, 16.5),
                    nrow = 4,
                    byrow = TRUE
)

ocBmatrix <- matrix(data = c(-23.853065, 30.593651, -9.74629544, 3.0057072,
                             11.099906, -43.237662, 40.818768, -8.6810122,
                             -3.3228457, 38.356814, -125.40927, 90.375305,
                             -33.675598, 152.37521, -311.19961, 192.5),
                    nrow = 4,
                    byrow = TRUE
    
)
B <- Matrix::bdiag(B,B,B)
A <- Matrix::bdiag(A,A,A)
gamma_i <- c(rep(-1.5,3), rep(-0.5, 3), rep(1,3))
oc.points <- c(0.36311746, 0.67718628, 0.89975799, 1)
expand_w_i <- rep(oc.points[1:3],3)
    #as.vector(sapply(oc.points[1:3], function(x) rep(x,3)))
expand_gamma_i <- c(rep(-1.5,3), rep(-0.5, 3), rep(1,3))
int.weights <- c(0.04567809, 0.12589839, 0.13397916, 1/36)
