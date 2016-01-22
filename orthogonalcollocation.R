#collocation matricies
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
gamma_i <- c(rep(-3,3), rep(-1,3), rep(2,3))

#set a orthogonal collocation function for stream
set.ocf <- function(surface, 
                    t, p, C_total,
                    preExp, Ea, alpha,
                    intrPorosity, voidage, rParticle){
    if (!element.match(names(surface), elementNames[1:3]))
        stop("No names for surface concentrations provided")
    ocf <- function(intcolpoints){
        kin.const <- pr.kinetic.constant(preExp, Ea, t)
        eq.const <- pr.equilibrium.constant(t,p)
        xMatr <- matrix(intcolpoints, ncol = 3)
        xMatr <- rbind(xMatr, surface)
        dimnames(xMatr) <- list(NULL,c("hydrogen", "nitrogen","ammonia"))
        #concentration at pellet surface
        at_surface <- dim(xMatr)[1]
        #fugacities
        fgc.cf <- sapply(fugacity.list[c("hydrogen", "nitrogen","ammonia")],
                         function(x) do.call(x, args = list(t=t+273.13, p=p))
        )
        fgc <- xMatr%*%diag(fgc.cf) * p
        colnames(fgc) <- c("hydrogen", "nitrogen","ammonia")
        #second derivative
        #d2ydx2 <- apply(B %*% xMatr, 1, sum)
        d2ydx2 <- as.vector(B %*% c(xMatr))
        #first derivative
        #dydx <- apply(A %*% xMatr, 1, sum)
        dydx <- as.vector(A %*% c(xMatr))
        #2rd term in equation
        term2 <- -2 / (gamma_i + 2*c(xMatr[-at_surface,])) * dydx ^2
        #3rd term in equation
        term3 <- 2 * (1 / w_i) * dydx
        ediff <- apply(xMatr[-at_surface,], 1,
                       function(fgc)
                           pr.effective.diffusivity(
                               pr.diffusivity(
                                   pr.std.diffusivity(xH2 = fgc[["hydrogen"]],
                                                      xN2 = fgc["nitrogen"],
                                                      xNH3 = fgc["ammonia"]),
                                   t, p
                               ),
                               intrPorosity
                           ) 
                       )
        ediff <- c(t(ediff))
#         ediff <- pr.effective.diffusivity(
#             pr.diffusivity(
#                 pr.std.diffusivity(xH2 = xMatr[at_surface,"hydrogen"],
#                                    xN2 = xMatr[at_surface, "nitrogen"],
#                                    xNH3 = xMatr[at_surface, "ammonia"]),
#                 t, p
#             ),
#             intrPorosity
#         )
        #names(ediff) <- c("hydrogen", "nitrogen","ammonia")
        #ediff <- c(sapply(ediff, function(x) rep(x,times = 3)))
        #reaction rate
        rate <- apply(fgc[-at_surface,], 1, function(fgc) pr.reaction.rate.NH3(kin.const, eq.const, fgc, alpha)) 
        #4th term
        term4 <- rParticle^2 / (C_total * ediff * (1 - voidage)) * rate * (gamma_i + 2*c(xMatr[-at_surface,]))
        #result
        return(d2ydx2 + term2 + term3 + term4)
    }
    return(ocf)
}
