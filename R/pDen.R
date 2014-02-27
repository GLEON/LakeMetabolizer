# Density Plot Color Scheme
# Densities will be estimated and plottted for each column of val (rows are multiple observations of same variable)
pDen <- function(vals=NULL, mu=NULL, sig=NULL){
	if(is.null(vals) & (is.null(mu)|is.null(sig))) stop("Must provide vals or (mu and sig), but not both.")
	if(any(sig<=0)) stop("Need positive sig")
	
	N <- length(mu)
	
	if(is.null(vals)){
		vals0 <- rnorm(n=100*N, mean=mu, sd=sig)
		vals <- matrix(vals0, ncol=N, byrow=TRUE)	
	}
	limX <- range(vals)
	
	dens <- apply(vals, 2, function(x)density(x, from=limX[1], to=limX[2])[c("x","y")])
	
	dX <- lapply(dens, function(z)z$x)
	dY <- lapply(dens, function(z)z$y)
	limY <- range(dY)
	
	cLine <- rainbow(n=N)
	cFill <- rgb(t(col2rgb(cLine, alpha=TRUE)), alpha=35, maxColorValue=255)
	
	# Need to be sure that the smallest and largest densities are 0 so that the bottom border of polygons are at 0 line
	xF <- 0.01 * diff(limX) # a "factor" by which to extend the range of X
	xA <- limX + c(-1,1)*xF # "add" this "adjustment" to the start and end of the dX

	dev.new(width=5, height=3.5)
	par(mar=c(1.5, 1.5, 0.5, 0.5), ps=10, cex=1, mgp=c(1.5, 0.1, 0), tcl=0.25, las=1)
	plot(c(xA[1],dX[[1]],xA[2]), c(0,dY[[1]],0), type="l", col=cLine[1], xlab="", ylab="", xlim=limX, ylim=limY, lwd=2)
	polygon(c(xA[1],dX[[1]],xA[2]), c(0,dY[[1]],0), col=cFill[1], border=NA)
	for(i in 2:N){
		polygon(c(xA[1],dX[[1]],xA[2]), c(0,dY[[i]],0), col=cFill[i], border=cLine[i], lwd=2)
	}
}

# Example (that simulates fake data):
# set.seed(3)
# pDen(mu=runif(5, -100, 100), sig=runif(5, 1, 15))


