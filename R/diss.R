
################################################################################
#######################  AUTOCORRELATION AND PARTIAL ###########################
################################################################################


#check if series have equal length, a requisite of some functions
.check.equal.length.ts <- function(x,y) {
    if (length(x) != length(y)) {
        stop("Time series must have the same length")
    }	
    
}

#weighted distance of acf and pacf coefficients
.internal.autocorr.dist <- function( rhox, rhoy, p=NULL, omega=NULL ) {
    if (is.null(omega)) { #if there is no weighting matrix
        if (!is.null(p)) { #check if there is gemoetrical decay parameter
            omega <- diag(p*(1-p)**(1:length(rhox)))
        }
        else { #no weightinh matrix and no geomtrical decay parameter, use identity matrix
            omega <- diag(length(rhox))
        }
    }
    t(rhox - rhoy) %*% omega %*% (rhox - rhoy) #weighted euclidean distance
}


diss.ACF <- function( x, y ,  p = NULL,  omega=NULL, lag.max = 50  ) {
    rhox <- acf(x, lag.max=lag.max, plot=FALSE)$acf[-1]
    rhoy <- acf(y, lag.max=lag.max, plot=FALSE)$acf[-1]
    .internal.autocorr.dist( rhox, rhoy, p, omega)
}



diss.PACF <- function(x, y, p = NULL, omega=NULL, lag.max=50) {
    rhox <- as.vector(pacf(x, plot=FALSE)$acf)
    rhoy <- as.vector(pacf(y, plot=FALSE)$acf)
    .internal.autocorr.dist( rhox, rhoy, p, omega)
}

#######################################################
##########  distance Piccolo  #########################
#######################################################

diss.AR.PIC <- function(x, y, order=NULL) {
    order.x <- NULL
    order.y <- NULL
    if (is.null(order)) { #order for the ARIMA modeling
        order <- rbind( c(NA,NA,NA), c(NA,NA,NA))
    }
    if (!sum( is.na(order[1,]))) {
        order.x <- order[1,]
    }
    if (!sum( is.na(order[2,]))) {
        order.x <- order[2,]
    }
    
    PIx <- NULL
    if (is.null(order.x)) { #no ARIMA order, use AR
        PIx <- ar(x)$ar
    } else { 
        if ((order.x[1]) < 1) stop("The arima order must have AR coefficients, they are used for the distance")
        arim <- arima(x, order.x)
        PIx <- arim$coef[1:order.x[1]] #get the AR coeff off the arima model
    }
    PIy <- NULL
    if (is.null(order.y)) { #no ARIMA order, use AR
        PIy <- ar(y)$ar
    } else {
        if ((order.y[1]) < 1) stop("The arima order must have AR coefficients, they are used for the distance")
        arim <- arima(y, order.y)
        PIy <- arim$coef[1:order.y[1]] #get the AR coeff off the arima model
    }
    k <- max(c(length(PIx), length(PIy))) #get the maximun order
    if (k < 1) { #if no AR model found, impose an AR(1) model
        warning("Error in automatic selection of AR order, 0 selected by AIC, forcing 1")
        ARx <- ar(x, aic=FALSE, order.max = 1)
        ARy <- ar(y, aic=FALSE, order.max = 1)
        PIx <- ARx$ar
        PIy <- ARy$ar
        k <- min(c(length(PIx), length(PIy)))
        if (k < 1) stop("Could not find any AR coefficients")
    }
    PRIMAx <- rep(0,k) #fill with zeroes to the greates AR order between series x and y (k)
    if (length(PIx) > 0) {
        PRIMAx[1:length(PIx)] <- PIx
    }
    PRIMAy <- rep(0,k)
    if (length(PIy) > 0) {
        PRIMAy[1:length(PIy)] <- PIy
    }
    as.numeric( dist(rbind(PRIMAx,PRIMAy)) ) #compute the euclidean distance between the zero padded AR coefficients
}


#prototyping multiple parameter per series distance
multidiss.AR.PIC <- function( series, order=NULL) {
    distances <- matrix(0, nrow(series), nrow(series))
    rownames(distances) <- rownames(series)
    
    for ( i in 1: (nrow(series)-1) ) {
        for (j in (i+1):nrow(series) ) {
            distance <- diss.AR.PIC(series[ i,], series[j,], rbind(order[i,], order[j,]) )
            distances[ i, j] <- distance
            distances[ j, i] <- distance
        }
    }
    as.dist((distances))
}


##################################################
########### distance Maharaj  ####################
##################################################

#regression model
maharajahextended <- function( x, k ) {
    
    X <- x[-(1:k)]
    
    TT <- length(x)
    
    Wx <- matrix(ncol=k, nrow=(TT-k))
    
    for (i in 1:(TT -k) ) {
        Wx[i,] <- x[(k +i - 1):(i)]
    }
    result <- list()
    result$X <- X
    result$Wx <- Wx
    result
}

#extended distance see reference article in the documentation
distance.MAH.EXT <- function( x, y, k) {
    MX <- maharajahextended(x, k)
    MY <- maharajahextended(y, k)
    w <- dim(MX$Wx)[2]
    h <- dim( MX$Wx)[1]
    bigW <- matrix(0, nrow=2*h, ncol=2*w)
    for ( j in 1:w) {
        for (i in 1:h) {
            bigW[i,j] <- MX$Wx[i,j]  
        }
    }
    for ( j in 1:w) {
        for (i in 1:h) {
            bigW[h+i,w+j] <- MY$Wx[i,j]  
        }
    }
    Epsil <- cov(cbind(x,y))
    bigW
    Epsil
    Iden <- diag(1, length(x) - k)
    Iden
    V <- kronecker(Epsil, Iden  )
    IV <- solve(V)
    tryCatch( {
        PI <- solve(t(bigW) %*% IV %*% bigW) %*% t(bigW) %*% IV %*% c(MX$X, MY$X)
        R <- cbind( diag(1, k), diag(-1,k) )
        result <- list()
        result$statistic <- t(R %*% PI) %*% solve( R %*% (t(bigW) %*% IV %*% bigW ) %*% t(R)) %*% ( R%*%PI)
        result$p_value <- pchisq(result$statistic, k, lower.tail=F)
        result
    }, error = function(e) { if (k>1) {
        distance.MAH.EXT(x,y,k-1)
    }
                             else { stop("Could not find valid AR order")}
    })
}



distance.MAH.SIMP = function( x, y ) {
    ARx <- ar(x)
    ARy <- ar(y)
    PIx <- ARx$ar
    PIy <- ARy$ar
    k <- max(c(length(PIx), length(PIy)))
    if (k < 1) {
        warning("Error in automatic selection of AR order, 0 selected by AIC, forcing 1")
        ARx <- ar(x, aic=FALSE, order.max = 1)
        ARy <- ar(y, aic=FALSE, order.max = 1)
        PIx <- ARx$ar
        PIy <- ARy$ar
        k <- min(c(length(PIx), length(PIy)))
        if (k < 1) stop("Could not find any AR coefficient")
    }
    PRIMAx <- rep(0,k)   #fill with zeroes
    if (length(PIx) > 0) {
        PRIMAx[1:length(PIx)] <- PIx
    }
    PRIMAy <- rep(0,k)
    if (length(PIy) > 0) {
        PRIMAy[1:length(PIy)] <- PIy
    }
    
    covx <- acf(x, lag.max=k-1, type="covariance", plot=FALSE)$acf
    covy <- acf(y, lag.max=k-1,type="covariance", plot=FALSE)$acf
    Rx <- matrix(nrow=length(covx),ncol=length(covx))
    Ry <- Rx
    for (i in 1:length(covx) ) {
        indices <- (c(i:1, 2:(length(covx)-i+1)))
        Rx[i,] <- covx[ indices[1:length(covx)]]
        Ry[i,] <- covy[ indices[1:length(covx)]]
    }
    V <- (solve(solve(Rx)*(ARx$var.pred) + solve(Ry)*(ARy$var.pred)))
    dif <- (PRIMAx - PRIMAy)
    
    D <- sqrt(length(x)) * ( dif %*% V %*% dif)
    
    list(statistic=D, p_value=pchisq(D, k, lower.tail=F))
}


diss.AR.MAH = function( x, y, dependence = FALSE ) {
    .check.equal.length.ts(x,y)
    ARx <- ar(x)
    ARy <- ar(y)
    PIx <- ARx$ar
    PIy <- ARy$ar
    k <- max(c(length(PIx), length(PIy)))
    if (k < 1) { #calculate the greater AR order found, if order 0 found, impose order 1
        warning("Error in automatic selection of AR order, 0 selected by AIC, forcing 1")
        ARx <- ar(x, aic=FALSE, order.max = 1)
        ARy <- ar(y, aic=FALSE, order.max = 1)
        PIx <- ARx$ar
        PIy <- ARy$ar
        k <- min(c(length(PIx), length(PIy)))
        if (k < 1) stop("Could not find any AR coefficient")
    }
    if (dependence) {
        distance.MAH.EXT(x, y, k)
    }
    else {
        distance.MAH.SIMP(x,y)
    }
    
}

######################################################
########  PERIODOGRAM BASED DISTANCES  ###############
######################################################

diss.PER <- function(x,y, logarithm=FALSE, normalize=FALSE) {
    .check.equal.length.ts(x,y)
    Ix <- spec.pgram(x, plot=F)$spec
    Iy <- spec.pgram(y, plot=F)$spec
    if (normalize) {
        Ix <- Ix/var(x)
        Iy <- Iy/var(y)
    }
    if (logarithm) {
        Ix <- log(Ix)
        Iy <- log(Iy)
    }
    dist(rbind(Ix,Iy))/(length(Ix))
}


diss.INT.PER <- function(x,y, normalize=TRUE) {
    Ix <- spec.pgram(x, plot=F)
    Iy <- spec.pgram(y, plot=F)
    Cx <- 1
    Cy <- 1
    if (normalize) {
        Cx <- sum(Ix$spec)
        Cy <- sum(Iy$spec)
    }
    sum ( abs(cumsum(Ix$spec)/Cx - cumsum(Iy$spec)/Cy) )
}


################################################
###  SPECTRAL DENSITY APPROXMATION DISTANCES ###
################################################


### maximum likelihood functions ###

#kernel function
funcionKh <- function( value, h) {
    value <- value/h
    dnorm( -(value**2) ) / h
}

#function to be optimized for maximum likelihood, see referenced papers
Spectral.AB <- function (  ABvec, lambda, Yks, lambdas, h) {
    acum <- 0
    a <- ABvec[1]
    b <- ABvec[2]
    
    -sum ( ( -exp(Yks - a -b*(lambdas - lambda) ) + Yks - a -b*(lambdas - lambda) ) * funcionKh( lambdas - lambda, h) )
    
}

likelihood.optim <- function(  lambda, Yks, lambdas, h) {
    
    startA = lambdas[1]
    startB = 0
    
    optim(c(startA,startB), Spectral.AB, lambda=lambda, Yks = Yks, lambdas=lambdas, method="L-BFGS-B",lower = c(min(Yks) , -101), upper= c(max(Yks),101), h =h)$par[1]
}
### end of maximul likelihood


#diveregence function
divergenceW <- function( x, alpha) {
    if ((alpha > 0) & (alpha< 1)) {
        log(alpha*x + 1 - alpha) - alpha*log(x)
    }
    else {
        stop("condition 0 < alpha < 1 not fulfilled")
    }
}

simetricDivergenceW <- function(x,alpha) {
    divergenceW(x,alpha) + divergenceW(1/x,alpha)  
}

#plot the soothed spectral density
plotsmoothspec <- function ( lambdasX, YksX, lambdasY, YksY, myf, hX, hY , method="Maximum Likelihood") {
    
    baseX <- seq(min(lambdasX) + 0.001, max(lambdasX) - 0.001, length.out=min(500, 3*length(lambdasX)) )
    specX <- NULL
    if (pmatch(method , c("Maximum Likelihood", "Least Squares")) == 2) {
        specX <- (myf(YksX, lambdasX,  baseX, hX))
    } else {
        specX <- exp(myf(baseX, YksX, lambdasX, hX))
    }
    baseY <- seq(min(lambdasY)+ 0.001, max(lambdasY)- 0.001, length.out=min(500, 3*length(lambdasY)) )
    specY <- NULL
    if (pmatch(method , c("Maximum Likelihood", "Least Squares")) == 2) {
        specY <- (myf(YksY,lambdasY, baseY, hY))
    } else {
        specY <- exp(myf(baseY, YksY, lambdasY, hY))
    }
    plot.default(baseX, specX,type="l", col="red", ylim=c( min(specX,specY), max(specX, specY) ) ,
                 main=paste("Approx. spectral density by ", method), xlab="frequency",ylab="spectrum")
    lines(baseY, specY,type="l", col="blue")
    legend("topright", pch=16, col=c("red", "blue"), legend= c("x", "y") )
    
    
}


distance.W.LK <- function(x,y, alpha, plot=FALSE) {
    myf <- Vectorize(likelihood.optim,"lambda") #needed for likelihood.optim to accept a vector, required for integrate
    
    YksX <- log(spec.pgram(x,plot=FALSE)$spec)
    lambdasX <- spec.pgram(x,plot=FALSE)$freq
    YksY <- log(spec.pgram(y,plot=FALSE)$spec)
    lambdasY <- (spec.pgram(y,plot=FALSE)$freq)
    
    hX <- 0.93*dpill(lambdasX, YksX)
    hY <- 0.93*dpill(lambdasY, YksY)
    
    integrateaux <- function( lambda ) {
        xx <- exp(myf(lambda, YksX, lambdasX, hX))
        yy <- exp(myf(lambda, YksY, lambdasY, hY))
        simetricDivergenceW(  xx / yy, alpha)
    }
    
    tryCatch( {
        a <- integrate(integrateaux, min(lambdasX), max(lambdasX))
    }, error = function (e) {
        warning("Failed approximation with window from plug.in method, increasing window...")
        hX <- 2*hX
        hY <- 2*hY
        a <- integrate(integrateaux, min(lambdasX), max(lambdasX))
    })
    
    if (plot) {
        plotsmoothspec(lambdasX, YksX, lambdasY, YksY, myf, hX, hY)
    }
    
    a$value
}


distance.W.DLS <- function(x, y, alpha, plot=FALSE) {
    
    YksX <- (spec.pgram(x,plot=FALSE)$spec)
    lambdasX <- spec.pgram(x,plot=FALSE)$freq
    YksY <- (spec.pgram(y,plot=FALSE)$spec)
    lambdasY <- (spec.pgram(y,plot=FALSE)$freq)
    hX <- dpill(lambdasX, YksX)
    hY <- dpill(lambdasY, YksY)
    myf <- function(Yk, lambdas, lambda, h) {
        d <- data.frame( Yk)
        d$lambdas <- lambdas
        lp <- locpol(Yk~lambdas, d, bw=h,kernel=gaussK, xeval=lambda )
        ord <- order(lambda) #trick to get the original order of the lambas, locpol sorts the input vector
        ord2 <- order(ord)   #second part of the trick
        lp$lpFit$Yk[ord2]
    }
    
    integrateaux <- function( lambda ) {
        xx <- (myf(YksX, lambdasX, lambda, hX))
        yy <- (myf(YksY, lambdasY, lambda, hY))
        xx[xx<0.0001] <- 0.0001
        yy[yy<0.0001] <- 0.0001
        simetricDivergenceW(  xx / yy, alpha)
    }
    lambdas <- spec.pgram(x, plot=F)$freq
    a <- integrate(integrateaux, min(lambdas), max(lambdas), subdivisions=100)
    if (plot) {
        plotsmoothspec(lambdasX, YksX, lambdasY, YksY, myf, hX, hY, "Least Squares")
    }
    
    a$value
}


diss.SPEC.LLR <- function(x,y, alpha=0.5, method="DLS", plot=FALSE) {
    .check.equal.length.ts(x,y)
    typedist = 0
    type <-  (pmatch(method, c("DLS", "LK" )))
    if (is.na(type)) {
        stop(paste("Unknown method", method))
    } else if (type == 1) {
        typedist <- distance.W.DLS(x,y, alpha, plot)
    }
    else if (type == 2) {
        typedist <- distance.W.LK(x,y, alpha, plot)
    }
    typedist
}



diss.SPEC.GLK <- function(x,y, plot=FALSE) {
    .check.equal.length.ts(x,y)
    myf <- Vectorize(likelihood.optim,"lambda") #needed for likelihood.optim to accept a vector, required for integrate
    
    YksX <- log(spec.pgram(x,plot=FALSE)$spec)
    lambdasX <- spec.pgram(x,plot=FALSE)$freq
    YksY <- log(spec.pgram(y,plot=FALSE)$spec)
    lambdasY <- (spec.pgram(y,plot=FALSE)$freq)
    
    hX <- 0.93*dpill(lambdasX, YksX)
    hY <- 0.93*dpill(lambdasY, YksY)
    
    Z <- YksX - YksY
    mx <- myf(lambdasX, YksX, lambdasX, hX )
    my <- myf(lambdasY, YksY, lambdasY, hY )
    mu <- mx- my
    
    if (plot) {
        plotsmoothspec(lambdasX, YksX, lambdasY, YksY, myf, hX, hY)
    }
    #distance GLK
    sum(Z - mu - 2*log(1 + exp(Z - mu))) - sum( Z - 2*log(1 + exp(Z)))
}



#distancia ISD
diss.SPEC.ISD <- function(x,y, plot=FALSE) {
    .check.equal.length.ts(x,y)
    myf <- Vectorize(likelihood.optim,"lambda") #needed for likelihood.optim to accept a vector, required for integrate
    
    YksX <- log(spec.pgram(x,plot=FALSE)$spec)
    lambdasX <- spec.pgram(x,plot=FALSE)$freq
    YksY <- log(spec.pgram(y,plot=FALSE)$spec)
    lambdasY <- (spec.pgram(y,plot=FALSE)$freq)
    
    hX <- 0.93*dpill(lambdasX, YksX)
    hY <- 0.93*dpill(lambdasY, YksY)
    
    integraISDaux <- function(lambda) {
        (myf(lambda, YksX, lambdasX  , hX) - myf(lambda,YksY, lambdasY,  hY))**2
    }
    tryCatch( {
        a <- integrate(integraISDaux, min(lambdasX), max(lambdasX) )$value
    }, error = function(e) {
        hX <- 2*hX
        hY <- 2*hY
        a <- integrate(integraISDaux, min(lambdasX), max(lambdasX) )$value
    })
    if (plot) {
        plotsmoothspec(lambdasX, YksX, lambdasY, YksY, myf, hX, hY)
    }  
    a
}


#########################################################
##############  distance CEPSTRAL  ######################
#########################################################


.calc.cepstral.coef <- function( ARx, h ) {
    CEPSTRALx <- 1:h
    CEPSTRALx[1] <- -ARx[1]
    if (length(ARx) >= 2) {
        for (i in 2:length(ARx)) {
            acum <- 0
            for (m in 1:(i-1) ) {
                acum <- acum + ( 1- m/i)*ARx[m]*CEPSTRALx[i-m]
            }
            CEPSTRALx[i] <- -ARx[i] + -acum
        }
    }
    
    if (h > length(ARx)) {
        for (i in (length(ARx)+1):h) {
            acum <- 0
            for (m in 1:length(ARx)) {
                acum <- acum + (1 - m/i)*ARx[m]*CEPSTRALx[i-m]
            }
            CEPSTRALx[i] <- -acum
        }
    }
    CEPSTRALx
    
}


cepstral <- function(x, h, order=NULL, seasonal) {
    ARx <- NULL
    valid.order <- !sum(is.na(order))
    SAR <- NULL
    if (!valid.order) {
        ARx <- ar(x,order.max=min(length(x)-1,h))
    } else {
        if ((order[1]) < 1) stop("The arima order must have AR coefficients, they are used for the distance")
        arim <- arima(x, order, seasonal)
        ARx$ar <- arim$coef[1:order[1]]
        if (seasonal[[1]][1]>0) { #the seasonal part as in kalpakis
            SAR <- arim$coef[ (1 + order[1] + order[2]):(1 + order[1] + order[2]+ seasonal[[1]][1]) ]
        }
    }
    
    if (length(ARx$ar) < 1) {
        warning("Cepstral distance, error on the selection of the AR order, 0 by AIC, forcing 1")
        ARx  <- ar(x, aic=FALSE, order.max=1)
        if (length(ARx)<1) {
            stop("Could not find any AR coefficient")      
        }
    }
    
    cepst.SAR <- NULL #compute the sesonal part
    if (is.null(SAR)) {
        cepst.SAR <- rep.int(0, h)
    }
    else {
        cepst.SAR <- .calc.cepstral.coef(SAR, h)
    }
    
    rbind( .calc.cepstral.coef(ARx$ar, h), cepst.SAR )
    
}


diss.AR.LPC.CEPS <- function(x, y, k=50 , order=NULL, seasonal=list( list(order = c(0, 0, 0), period = NA), list(order = c(0, 0, 0), period = NA)) ) {
    if (is.null(order)) {
        order <- rbind(c(NA,NA,NA), c(NA,NA,NA))
    }
    cpx <- cepstral(x,k, order[1,], seasonal[[1]])
    cpy <- cepstral(y,k, order[2,], seasonal[[2]])
    as.numeric(dist(rbind(cpx[1,],cpy[1,])) + dist(rbind(cpx[2,],cpy[2,])))
}


#prototyping multiple parameter per series distance
multidiss.AR.LPC.CEPS <- function( series, k=50, order=NULL, seasonal=NULL) {
    distances <- matrix(0, nrow(series), nrow(series))
    rownames(distances) <- rownames(series)
    if (is.null(seasonal)) {
        seasonal <- rep( list( list(order=c(0,0,0), period=NA)), nrow(series) )
    }
    for ( i in 1:(nrow(series)-1) ) {
        for (j in (i+1):nrow(series) ) {
            distance <- diss.AR.LPC.CEPS(series[ i,], series[j,], k , rbind(order[i,], order[j,]) , list(seasonal[[i]], seasonal[[j]]) )
            distances[ i, j] <- distance
            distances[ j, i] <- distance
        }
    }
    as.dist((distances))
}



#############################################################################
#################   Temporal Correlation Distance   #########################
#############################################################################

##CHOUAKRIA-DOUZAL

corrtemporder1 <- function (x, y) {
    p <- length(x)
    sum((x[2:p] - x[1:(p-1)]) * (y[2:p] - y[1:(p-1)])) / ( sqrt( sum((x[2:p] - x[1:(p-1)])^2) ) * sqrt( sum((y[2:p] - y[1:(p-1)])^2) ))
}

diss.CORT <- function( x, y, k=2, deltamethod="Euclid") {
    .check.equal.length.ts(x,y)
    corrt <- corrtemporder1(x,y)
    type <-  (pmatch(deltamethod, c("Euclid", "Frechet", "DTW")))
    typedist <- 0
    if (is.na(type)) {
        stop(paste("Unknown method", deltamethod))
    } else if (type == 1) {
        typedist <- as.numeric( dist(rbind(x,y)) )
    }
    else if (type == 2) {
        typedist <- distFrechet(x,y)
    }
    else if (type == 3) {
        typedist <- dtw(x,y,distance.only=T)$distance
    }
    
    (2/( 1+ exp(k*corrt)))*typedist
    
}



##################################################
######  maharaj clustering algorithm #############
##################################################
#input, distance matrix
pvalues.clust <- function(pvalues,significance) {
    distancias <- pvalues
    significacion <- significance
    distancias <- as.matrix(distancias)
    tam <- dim(distancias)[1]
    distancias
    combinaciones <- combn(tam,2)
    plandist <- 1:ncol(combinaciones)
    #create a vector with the distances
    for (i in 1:length(plandist)) {
        plandist[i] <- distancias[combinaciones[1,i], combinaciones[2,i]]
    }
    
    ord <- order(plandist, decreasing=TRUE)
    plandist <- plandist[ord]
    combinaciones <- combinaciones[,ord]
    
    is_in_setlist <- function( element, setlist) {
        for (i in 1:length(setlist)) {
            if (element %in% setlist[[i]]) {
                return(TRUE)
            }
        }
        return(FALSE)
    }
    
    find_which_set <- function ( element, setlist) {
        for (i in 1:length(setlist)) {
            if (element %in% setlist[[i]]) {
                return(i)
            }
        }
        return(0)
    }
    
    
    is_pvalue_of_element_less_than_significance_with_the_set <- function( element, setid, setlist, significance,distances) {
        sum(distances[element, setlist[[setid]]] < significance) == length(setlist[[setid]])
    }
    
    add_to_set <- function ( element, setindex, setlist) {
        setlist[[setindex]] <- union(setlist[[setindex]], element)
        setlist
    }
    
    create_new_set <- function( element, setlist) {
        setlist[[length(setlist) + 1]] <- element
        setlist
    }
    
    is_all_series_already_in_a_cluster <- function(series, setlist) {
        already = TRUE
        for (i in series) {
            already <- already & (find_which_set(i,setlist) != 0)
        }
        already
    } 
    
    are_pvalues_of_all_pairs_across_clusters_greater_than_significance <- function( clusterone, clustertwo, setlist, significance, distances) {
        seriesone <- setlist[[clusterone]]
        seriestwo <- setlist[[clustertwo]]
        aregreater <- TRUE
        for (i in seriesone) {
            aregreater <- aregreater & (sum(( distances[i, seriestwo] > significance)) == length(seriestwo) )
        }
        for (i in seriestwo) {
            aregreater <- aregreater & (sum(( distances[i, seriesone] > significance)) == length(seriesone) )
        }
        aregreater
    }
    
    merge_sets <- function (setone, settwo, setlist) {
        setlist[[settwo]] <- c( setlist[[settwo]] ,setlist[[setone]])
        setlist[-setone]
    }
    
    if (plandist[1] < significacion) {
        grupos <- as.list(1:tam)
    } else {
        grupos <- list()
        grupos[[1]] <- combinaciones[,1]
        for (i in 2:length(combinaciones[1,])) {
            if (plandist[i] < significacion) { #p-value < significance = YES
                conj <- find_which_set( combinaciones[1,i], grupos )
                if (conj == 0) {
                    grupos <- create_new_set(combinaciones[1,i], grupos)
                }
                conj <- find_which_set( combinaciones[2,i], grupos )
                if (conj == 0) {
                    grupos <- create_new_set(combinaciones[2,i], grupos)
                } 
                #each remaining serie to its own cluster
                for (j in i:length(combinaciones[1,])) {
                    conj <- find_which_set( combinaciones[1,j], grupos )
                    if (conj == 0) {
                        grupos <- create_new_set(combinaciones[1,j], grupos)
                    }
                    conj <- find_which_set( combinaciones[2,j], grupos )
                    if (conj == 0) {
                        grupos <- create_new_set(combinaciones[2,j], grupos)
                    } 
                }
                break;
            } else { #p-value < significance = NO
                conj <- is_in_setlist( combinaciones[1,i], grupos )
                conj <- conj + is_in_setlist( combinaciones[2,i], grupos )
                if (conj < 2) {  #is each(ALL) series already in a cluster = NO
                    conj <- find_which_set( combinaciones[1,i], grupos )
                    if (conj > 0) { #one of the series in the pair already in a cluster = YES (x)
                        if (is_pvalue_of_element_less_than_significance_with_the_set( combinaciones[2,i], conj, grupos, significacion, distancias  )    ) {
                            grupos <- create_new_set(combinaciones[2,i], grupos)
                        } else {
                            grupos <- add_to_set( combinaciones[2,i], conj, grupos )
                        }
                    } else {
                        conj <- find_which_set( combinaciones[2,i], grupos )
                        if (conj > 0) {#one of the series in the pair already in a cluster = YES (y)
                            if (is_pvalue_of_element_less_than_significance_with_the_set( combinaciones[1,i], conj, grupos, significacion, distancias  )    ) {
                                grupos <- create_new_set(combinaciones[1,i], grupos)
                            } else {
                               grupos <- add_to_set( combinaciones[1,i], conj, grupos )
                            }
                        } else {#one of the series in the pair already in a cluster = NO
                            grupos <- create_new_set( combinaciones[1,i], grupos ) # create a new cluster with the two series
                            setid <- find_which_set( combinaciones[1,i], grupos )
                            grupos <- add_to_set( combinaciones[2,i], setid, grupos )
                        }
                    } 
                } else {  #is each series already in a cluster = YES
                    conj1 <- find_which_set( combinaciones[1,i], grupos )
                    conj2 <- find_which_set( combinaciones[2,i], grupos )
                    if (are_pvalues_of_all_pairs_across_clusters_greater_than_significance(conj1, conj2, grupos, significacion, distancias)) {
                        #merge
                        grupos <- merge_sets( conj1, conj2, grupos)
                    } else {
                        #do nothing
                    }
                }
            }
        }  
    }
    
    result <- 1:tam
    for ( i in 1:length(grupos) ) {
        for (j in grupos[[i]]) {
            result[j] <- i
        }
    }
    result
}


#check if the maharaj algorithm is equivalent to hierarchical clustering with 
# a cut level
testIgualdadMaharajHCLUST <- function( distancias, pvalor) {
    clusterp <- pvalues.clust(as.dist(distancias), pvalor)
    clusterh <- (hclust(as.dist(-distancias), method="single"))
    print(clusterp)
    plot(clusterh)
}


####################################################################
############# Feature Extraction Based on Wavelet Transform ########
####################################################################

wavelet.feature.extraction <- function(series) {
    
    calcEnergies <- function( wavdecomp) {
        level <- length(wavdecomp$data) -1
        energyD <- rep.int(0,level)
        energyD[1] <- sum(wavdecomp$data[[1]]**2)
        for ( i in 1:level) {
            energyD[i] <- sum((wavdecomp$data[[i]])**2)
        }
        return (energyD)
    }
    
    
    max.level <- as.integer(ceiling(logb(length(series[1,]),base=2)))
    true.level <- as.integer(floor(logb(length(series[1,]),base=2)))
    if (max.level != true.level) {
        npad <- 2**max.level - length(series[1,])
        for (i in 1:npad) {
            series <- cbind(series, rep(0, nrow(series)))
        }
    }
    
    energies <- matrix(0, nrow=nrow(series), ncol = max.level)
    for ( i in 1:nrow(series) ) {
        wavdecomp <- wavDWT(series[i,], n.levels=max.level, wavelet="haar")
        energies[i,] <- calcEnergies(wavdecomp)
    }
    
    sumEnergies <- colSums(energies)
    
    final_level <- max.level
    for (i in 1:(max.level-1)) {
        if (sumEnergies[i] < sumEnergies[i+1]) {
            final_level <- i
            break
        }
    }
    wavdecomp <- wavDWT(series[1,], n.levels=final_level, wavelet="haar")
    out.series <- wavdecomp$data[[final_level+1]]
    for (i in 2:nrow(series)) {
        wavdecomp <- wavDWT(series[i,], n.levels=final_level, wavelet="haar")
        out.series <- rbind(out.series, wavdecomp$data[[final_level+1]])
    }
    out.series
}

diss.DWT <- function(series) {
    if ( length(dim(series)) == 2 ) {
        if ( dim(series)[1] < 2 ) {
            stop( "diss.DWR needs a minimum of 2 series to compute the distance, incorrect amount provided" )
        }
    } else {
        stop( "diss.DWR series matrix with incorrect dimensions" ) 
    }
    wt <- wavelet.feature.extraction( series )
    dist(wt)
}

##########################################################
############  CORRELATION BASED DISTANCES ################
##########################################################



diss.COR <- function(x,y, beta = NULL) {
    correl <- cor(x,y)
    if (is.null(beta)) {
        sqrt(2*(1- correl))
    } else {
        if (beta<0) {
            stop("beta must be greater than 0")
        }
        sqrt( ((1-correl)/(1+correl ))**beta )    
    }
}






########################################################################
################# MULTIPLE SERIES DISTANCE GENERATION ##################
########################################################################





#######################################################################
############### CLUSTER SIMILARITY INDEX ##############################
#######################################################################

#gravrilov similarity ratio, used in kalpakis
Sim <- function(Gi, Sj, G, S) {
    2* (sum ( (G==Gi) & (S==Sj) ) ) / ( sum(G==Gi) + sum(S==Sj))
}

cluster.evaluation <- function(G,S) {
    acum <- 0
    gclust <- unique(G)
    sclust <- unique(S)
    for (i in gclust) {
        maxS <- 0
        for (j in sclust) {
            ms <- Sim(i, j, G, S)
            if (ms > maxS) {
                maxS <- ms
            }
        }
        acum <- acum + maxS
    }
    acum/(length(gclust))
}


############################################################################
#######################   OLD STUFF (UNUSED)   #############################
############################################################################



#distancia basada en correlaciones cruzadas
distanciaCorCruLagK = function(x,y,k) {
    muX = mean(x)
    muY = mean(y)
    TT = length(x)
    sum( (x[1:(TT-k)] - muX)*( y[(1+k):TT] - muY)) / (sqrt( sum((x - muX)**2))*sqrt( sum((y - muY)**2)))
}

#distanciaCorCruLagK(x,y,18)

distanciaCorCruTotal = function(x,y) {
    k = length(x)-1
    denom = 0
    for (i in 1:k) {
        denom = denom + distanciaCorCruLagK(x,y,i)
        print(denom)
    }
    
    sqrt( (1 - distanciaCorCruLagK(x,y,0))/denom  )
    denom
}

#distanciaCorCruTotal(x,y)
