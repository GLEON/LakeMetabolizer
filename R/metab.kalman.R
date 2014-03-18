
# install.packages("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer", repos=NULL, type="source")
library("LakeMetabolizer")
source("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer/inProgress/ryanData.R")
library("LakeMetabolizer")


metab.kalman <- function(){
	
	# ==================
	# = Filter and Fit =
	# ==================
	# ==================
	# = Fit 1: optim() =
	# ==================
	guesses <- c(1E-4,1E-4,log(5),log(5))
	# KFnllDO(Params=c(1E-4,-1E-4,5,5), do.obs=data[,"do.obs"], do.sat=data[,"do.sat"], K=data[,"K"], Zmix=data[,"Zmix"], irr=data[,"irr"], wtr=data[,"wtr"])
	fit <- optim(guesses, fn=KFnllDO, do.obs=data[,"do.obs"], do.sat=data[,"do.sat"], K=data[,"K"], Zmix=data[,"Zmix"], irr=data[,"irr"], wtr=data[,"wtr"])
	pars0 <- fit$par
	pars <- c(pars0[1], pars0[2], exp(pars0[3:4]))
	
	# ====================
	# = Fit 2: DEoptim() =
	# ====================
	# deLow <- c(0, -10, log(0.1), log(0.1))
	# deHigh <- c(10, 0, log(1E3), log(1E3))
	# fitDE <- DEoptim(fn=KFnllDO, lower=deLow, upper=deHigh, control=list(trace=FALSE), do.obs=data[,"do.obs"], do.sat=data[,"do.sat"], K=data[,"K"], Zmix=data[,"Zmix"], irr=data[,"irr"], wtr=data[,"wtr"])$optim
	# pars0DE <- fit$par
	# pars <- c(pars0[1], pars0[2], exp(pars0[3:4]))
	
	# Answer for day 137 (resp forced -)
	# $par
	# 	[1]  0.01443165 -0.66574384  6.54485304  4.39412874
	# 
	# 	$value
	# 	[1] 627.9323
	# 
	# 	$counts
	# 	function gradient 
	# 	      39       NA 
	# 
	# 	$convergence
	# 	[1] 0
	# 
	# 	$message
	# 	NULL
	
	#Answer for day 319 (H and Q forced positive, but Resp not forced)
	# $par
	# [1]  3.202261e-05 -1.324314e-02 -6.038739e+00 -1.659718e+01
	# 
	# $value
	# [1] -228.5235
	# 
	# $counts
	# function gradient 
	#      289       NA 
	# 
	# $convergence
	# [1] 0
	# 
	# $message
	# NULL
	
	# ==========
	# = Smooth =
	# ==========
	smoothDO <- KFsmoothDO(pars, do.obs=data[,"do.obs"], do.sat=data[,"do.sat"], K=data[,"K"], Zmix=data[,"Zmix"], irr=data[,"irr"], wtr=data[,"wtr"])
	
	smoothDO2 <- KFsmoothDO(pars, do.obs=data[,"do.obs"], do.sat=data[,"do.sat"], K=data[,"K"], Zmix=data[,"Zmix"], irr=data[,"irr"], wtr=data[,"wtr"], Hfac=5)
	
	dev.new(height=6, width=3.5)
	par(mfrow=c(2,1), mar=c(2,2,0.5,0.5), ps=9, mgp=c(2,0.5,0), tcl=-0.4, family="Times")
	plot(data[,"do.obs"], type="l", lwd=5, col="gray")
	lines(smoothDO, col="red")
	plot(data[,"do.obs"], type="l", lwd=5, col="gray")
	lines(smoothDO2, col="red")
	
	GPP <- sum(pars[1]*data[,"irr"])
	R <- sum(pars[2]*log(data[,"wtr"]))
	
	
}
