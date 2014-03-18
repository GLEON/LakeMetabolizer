
# install.packages("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer", repos=NULL, type="source")
library("LakeMetabolizer")
source("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer/inProgress/ryanData.R")
library("LakeMetabolizer")


metab.kalman <- function(do.obs, do.sat, k.gas, z.mix, date.times, irr, wtr){
	# ==================
	# = Filter and Fit =
	# ==================
	# ==================
	# = Fit 1: optim() =
	# ==================
	guesses <- c(1E-4,1E-4,log(5),log(5))
	# KFnllDO(Params=c(1E-4,-1E-4,5,5), do.obs=data[,"do.obs"], do.sat=data[,"do.sat"], K=data[,"K"], Zmix=data[,"Zmix"], irr=data[,"irr"], wtr=data[,"wtr"])
	fit <- optim(guesses, fn=KFnllDO, do.obs=do.obs, do.sat=do.sat, k.gas=k.gas, z.mix=z.mix, irr=irr, wtr=wtr)
	pars0 <- fit$par
	pars <- c(pars0[1], pars0[2], exp(pars0[3:4]))
	
	# ====================
	# = Fit 2: DEoptim() =
	# ====================
	# deLow <- c(0, -10, log(0.1), log(0.1))
	# deHigh <- c(10, 0, log(1E3), log(1E3))
	# fitDE <- DEoptim(fn=KFnllDO, lower=deLow, upper=deHigh, control=list(trace=FALSE), do.obs=do.obs, do.sat=do.sat, k.gas=k.gas, z.mix=z.mix, irr=irr, wtr=wtr)$optim
	# pars0DE <- fit$par
	# pars <- c(pars0[1], pars0[2], exp(pars0[3:4]))
	
	# ==========
	# = Smooth =
	# ==========
	smoothDO <- KFsmoothDO(pars, do.obs=do.obs, do.sat=do.sat, k.gas=k.gas, z.mix=z.mix, irr=irr, wtr=wtr)
	
	# ====================================
	# = Use fits to calculate metabolism =
	# ====================================
	GPP <- sum(pars[1]*data[,"irr"])
	R <- sum(pars[2]*log(data[,"wtr"]))
	
	return(list("smoothDO"=smoothDO,"params"=pars, "metab"=c("GPP"=GPP,"R"=R)))
}

# =====================
# = Example use of KF =
# =====================
KFans <- metab.kalman(do.obs=data[,"do.obs"], do.sat=data[,"do.sat"], k.gas=data[,"k.gas"], z.mix=data[,"z.mix"], irr=data[,"irr"], wtr=data[,"wtr"])

# =========================================
# = Plot smoothed and orginal time series =
# =========================================
dev.new(height=6, width=3.5)
par(mfrow=c(2,1), mar=c(2,2.5,0.5,0.5), oma=c(0.5, 0, 0, 0), ps=9, mgp=c(1.5,0.2,0), tcl=-0.2, xpd=TRUE)
plot(data[,"do.obs"], type="l", lwd=5, col="gray", xlab="", ylab=bquote(DO))
legend("topleft", legend=c("obs", "smoothed"), col=c("gray","red"), lty=1, lwd=c(5, 1))
lines(KFans$smoothDO, col="red")
plot(diff(data[,"do.obs"]), type="l", lwd=5, col="gray", xlab="", ylab=bquote(italic(d)*DO))
mtext("time",side=1, line=-0.5, outer=TRUE)
lines(diff(KFans$smoothDO), col="red")



# BayesAns <- metab.bayesian(do.obs=data[,"do.obs"], do.sat=data[,"do.sat"], k.gas=data[,"k.gas"], z.mix=data[,"z.mix"], irr=data[,"irr"], wtr=data[,"wtr"])

