### R code from vignette source 'Case_Study_4.xRnw'

###################################################
### code chunk number 1: RUNFIRST
###################################################
require(MARSS)
options(prompt=" ", continue=" ")
library(xtable)
tabledir="tables/"


###################################################
### code chunk number 2: read.in.data
###################################################
# load the data
data(lakeWAplankton)
# use only the data from 1980 onward
dat.spp.1980 = lakeWAplankton[lakeWAplankton[,"Year"]>=1980,]
# create vector of phytoplankton group names
phytoplankton = c("Cryptomonas", "Diatoms", "Greens",
                   "Bluegreens", "Unicells", "Other.algae")
# get only the phytoplankton
dat.spp.1980 = dat.spp.1980[,phytoplankton]


###################################################
### code chunk number 3: transpose.data
###################################################
# transpose data so time goes across columns
dat.spp.1980 = t(dat.spp.1980)
# get number of time series
N.ts = dim(dat.spp.1980)[1]
# get length of time series
TT = dim(dat.spp.1980)[2] 


###################################################
### code chunk number 4: z.score
###################################################
Sigma = sqrt(apply(dat.spp.1980, 1, var, na.rm=TRUE))
y.bar = apply(dat.spp.1980, 1, mean, na.rm=TRUE)
dat.z = (dat.spp.1980 - y.bar) * (1/Sigma)
rownames(dat.z) = rownames(dat.spp.1980)


###################################################
### code chunk number 5: plotdata
###################################################
spp = rownames(dat.spp.1980)
par(mfcol=c(3,2), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
for(i in spp){
  plot(dat.z[i,],xlab="",ylab="Abundance index", bty="L", xaxt="n", pch=16, col="blue", type="b")
  axis(1,12*(0:dim(dat.spp.1980)[2])+1,1980+0:dim(dat.spp.1980)[2])
  title(i)
  }


###################################################
### code chunk number 6: set.up.Z
###################################################
Z.vals = list(
"z11",  0  ,  0  ,
"z21","z22",  0  ,
"z31","z32","z33",
"z41","z42","z43",
"z51","z52","z53",
"z61","z62","z63")
Z = matrix(Z.vals, nrow=N.ts, ncol=3, byrow=TRUE)


###################################################
### code chunk number 7: print.Z
###################################################
print(Z)


###################################################
### code chunk number 8: set.up.A
###################################################
A.vals = list("a1","a2","a3","a4","a5","a6")
A = matrix(A.vals, nrow=N.ts, ncol=1)


###################################################
### code chunk number 9: set.up.QR
###################################################
Q = B = diag(1,3)


###################################################
### code chunk number 10: set.up
###################################################
R.vals = list(
"r11",0,0,0,0,0,
0,"r22",0,0,0,0,
0,0,"r33",0,0,0,
0,0,0,"r44",0,0,
0,0,0,0,"r55",0,
0,0,0,0,0,"r66")

R = matrix(R.vals, nrow=N.ts, ncol=N.ts, byrow=TRUE)


###################################################
### code chunk number 11: print
###################################################
print(R)


###################################################
### code chunk number 12: set.up.R.short
###################################################
R = "diagonal and unequal"


###################################################
### code chunk number 13: set.up.U
###################################################
x0 = U = matrix(0, nrow=3, ncol=1)
x0 = U = "zero"


###################################################
### code chunk number 14: set.up.x0
###################################################
V0 = diag(5,3)


###################################################
### code chunk number 15: set.model
###################################################
dfa.model = list(Z=Z, A="zero", R=R, B=B, U=U, Q=Q, x0=x0, V0=V0)
cntl.list = list(maxit=50)


###################################################
### code chunk number 16: load.results.if.present
###################################################
# This is being done to speed up building the user guide
if("CS4--model.fits.Rdata" %in% dir()) {
	load("CS4--model.fits.Rdata")
	saved.res = TRUE
} else { saved.res = FALSE }


###################################################
### code chunk number 17: fit.data
###################################################
if(!saved.res) {
	fit.cntl.list = list(minit=200, maxit=1200, allow.degen=FALSE)
	kemz.3 = MARSS(dat.z, model=dfa.model, control=fit.cntl.list)
	}


###################################################
### code chunk number 18: fit.data.echo (eval = FALSE)
###################################################
## kemz.3 = MARSS(dat.z, model=dfa.model, control=cntl.list)


###################################################
### code chunk number 19: plotfits
###################################################
fit = kemz.3
spp = rownames(dat.z)
par(mfcol=c(3,2), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
for(i in 1:length(spp)){
	plot(dat.z[i,],xlab="",ylab="abundance index",bty="L", xaxt="n", ylim=c(-4,3), pch=16, col="blue")
	axis(1,12*(0:dim(dat.z)[2])+1,1980+0:dim(dat.z)[2])
	par.mat=coef(fit,type="matrix")
	lines(as.vector(par.mat$Z[i,,drop=FALSE]%*%fit$states+par.mat$A[i,]), lwd=2)
	title(spp[i])
	}


###################################################
### code chunk number 20: set.up.two.trends
###################################################
if(!saved.res) {
	model.list = list(m=2, R="diagonal and unequal")
	kemz.2 = MARSS(dat.spp.1980, model=model.list,
	    z.score=TRUE, form="dfa", control=fit.cntl.list)
	}


###################################################
### code chunk number 21: set.up.two.trends.echo (eval = FALSE)
###################################################
## model.list = list(m=2, R="diagonal and unequal")
## kemz.2 = MARSS(dat.spp.1980, model=model.list,
##     z.score=TRUE, form="dfa", control=cntl.list)


###################################################
### code chunk number 22: compare.mods.2n3
###################################################
print(cbind(model=c("3 trends", "2 trends"),
      AICc=round(c(kemz.3$AICc, kemz.2$AICc))),
      quote=FALSE)


###################################################
### code chunk number 23: set.up.many.trends.no.echo
###################################################
#This is being done to speed up building the user guide
if(!saved.res) {
	# set up forms of R matrices
	levels.R = c("diagonal and equal",
	              "diagonal and unequal",
	              "unconstrained")
	model.data = data.frame()
	# fit lots of models & store results
	# NOTE: this will take a long time to run because minit is big!
	for(R in levels.R) {
	    for(m in 1:(N.ts-1)) {
	        dfa.model = list(A="zero", R=R, m=m)
	        kemz = MARSS(dat.z, model=dfa.model, control=fit.cntl.list, form="dfa", z.score=TRUE)
	        model.data = rbind(model.data,
	                            data.frame(R=R,
	                                       m=m,
	                                       logLik=kemz$logLik,
	                                       K=kemz$num.params,
	                                       AICc=kemz$AICc,
	                                       stringsAsFactors=FALSE))
	        assign(paste("kemz", m, R, sep="."), kemz)
	        } # end m loop
	    } # end R loop
	}


###################################################
### code chunk number 24: set.up.many.trends.echo (eval = FALSE)
###################################################
## # set new control params
## cntl.list = list(minit=200, maxit=1200, allow.degen=FALSE)
## # set up forms of R matrices
## levels.R = c("diagonal and equal",
##               "diagonal and unequal",
##               "unconstrained")
## model.data = data.frame()
## # fit lots of models & store results
## # NOTE: this will take a long time to run!
## for(R in levels.R) {
##     for(m in 1:(N.ts-1)) {
##         dfa.model = list(A="zero", R=R, m=m)
##         kemz = MARSS(dat.z, model=dfa.model, control=cntl.list, 
##             form="dfa", z.score=TRUE)
##         model.data = rbind(model.data,
##                             data.frame(R=R,
##                                        m=m,
##                                        logLik=kemz$logLik,
##                                        K=kemz$num.params,
##                                        AICc=kemz$AICc,
##                                        stringsAsFactors=FALSE))
##         assign(paste("kemz", m, R, sep="."), kemz)
##         } # end m loop
##     } # end R loop


###################################################
### code chunk number 25: makemodeltable
###################################################
# calculate delta-AICc
model.data$delta.AICc = model.data$AICc - min(model.data$AICc)
# calculate Akaike weights
wt = exp(-0.5*model.data$delta.AICc)
model.data$Ak.wt = wt/sum(wt)
# sort results
model.tbl = model.data[order(model.data$AICc),-4]
# drop AICc from table
# calculate cumulative wts
model.tbl$Ak.wt.cum = cumsum(model.tbl$Ak.wt)
model.tbl = model.tbl[,-4]
tmpaln = "c" #figure out the number of cols automatically
for(i in 1:ncol(model.tbl)) {tmpaln = paste(tmpaln,"c",sep="")}
thetable = xtable(model.tbl,caption='Model selection results.', label='tab:tablefits', align=tmpaln, digits=c(1,1,1,1,1,2,2))
align(thetable) = "cp{3.5cm}p{0.7cm}p{1.5cm}p{1.75cm}cc"
print(thetable, type = "latex", file = paste(tabledir,"tablefit.tex",sep=""), include.rownames=FALSE,include.colnames=TRUE, caption.placement="top",table.placement="htp", sanitize.text.function = function(x){x},hline.after = c(-1,0,nrow(model.data)))


###################################################
### code chunk number 26: getbestmodel
###################################################
# get the "best" model
best.model = model.tbl[1,]
fitname = paste("kemz",best.model$m,best.model$R,sep=".")
best.fit = get(fitname)


###################################################
### code chunk number 27: varimax
###################################################
# get the inverse of the rotation matrix
H.inv = varimax(coef(best.fit, type="matrix")$Z)$rotmat


###################################################
### code chunk number 28: rotations
###################################################
# rotate factor loadings
Z.rot = coef(best.fit, type="matrix")$Z %*% H.inv   
# rotate trends
trends.rot = solve(H.inv) %*% best.fit$states


###################################################
### code chunk number 29: plotfacloadings
###################################################
spp = rownames(dat.z)
minZ = 0.05
ylims = c(-1.1*max(abs(Z.rot)), 1.1*max(abs(Z.rot)))
par(mfrow=c(2,2), mar=c(2,4,1.5,0.5), oma=c(0.4,1,1,1))
for(i in 1:best.model$m) {
	plot(c(1:N.ts)[abs(Z.rot[,i])>minZ], as.vector(Z.rot[abs(Z.rot[,i])>minZ,i]),
		type="h", lwd=2, xlab="", ylab="", xaxt="n", ylim=ylims, xlim=c(0,N.ts+1))
	for(j in 1:N.ts) {
		if(Z.rot[j,i] > minZ) {text(j, -0.05, spp[j], srt=90, adj=1, cex=0.9)}
		if(Z.rot[j,i] < -minZ) {text(j, 0.05, spp[j], srt=90, adj=0, cex=0.9)}
	abline(h=0, lwd=1, col="gray")
    } # end j loop
    mtext(paste("Factor loadings on trend",i,sep=" "),side=3,line=.5)
} # end i loop


###################################################
### code chunk number 30: plottrends
###################################################
# get ts of trends
ts.trends = t(trends.rot)
par(mfrow=c(ceiling(dim(ts.trends)[2]/2),2), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
# loop over each trend
for(i in 1:dim(ts.trends)[2]) {
	# set up plot area
	plot(ts.trends[,i],
		ylim=c(-1.1,1.1)*max(abs(ts.trends)), 
		type="n", lwd=2, bty="L", 
		xlab="", ylab="", xaxt="n", yaxt="n")
	# draw zero-line
	abline(h=0, col="gray")
	# plot trend line
	par(new=TRUE)
	plot(ts.trends[,i],
		ylim=c(-1.1,1.1)*max(abs(ts.trends)), 
		type="l", lwd=2, bty="L", 
		xlab="", ylab="", xaxt="n")
	# add panel labels
	mtext(paste("Trend",i,sep=" "), side=3, line=0.5)
	axis(1,12*(0:dim(dat.spp.1980)[2])+1,1980+0:dim(dat.spp.1980)[2])
	} # end i loop (trends)


###################################################
### code chunk number 31: plotbestfits
###################################################
par.mat=coef(best.fit, type="matrix")
fit.b = par.mat$Z %*% best.fit$states + matrix(par.mat$A, nrow=N.ts, ncol=TT)
spp = rownames(dat.z)
par(mfcol=c(3,2), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
for(i in 1:length(spp)){
  plot(dat.z[i,],xlab="",ylab="abundance index",bty="L", xaxt="n", ylim=c(-4,3), pch=16, col="blue")
  axis(1,12*(0:dim(dat.z)[2])+1,1980+0:dim(dat.z)[2])
  lines(fit.b[i,], lwd=2)
  title(spp[i])
  }


###################################################
### code chunk number 32: set-up-covar
###################################################
temp = t(lakeWAplankton[lakeWAplankton[,"Year"]>=1980,"Temp"])
TP = t(lakeWAplankton[lakeWAplankton[,"Year"]>=1980,"TP"])


###################################################
### code chunk number 33: fit.covar
###################################################
if(!saved.res) {
	model.list=list(m=4, R="unconstrained")
	kemz.temp = MARSS(dat.spp.1980, model=model.list, z.score=TRUE, form="dfa",
	   control=fit.cntl.list, covariates=temp)
	kemz.TP = MARSS(dat.spp.1980, model=model.list, z.score=TRUE, form="dfa",
	   control=fit.cntl.list, covariates=TP)
	kemz.both = MARSS(dat.spp.1980, model=model.list, z.score=TRUE, form="dfa",
	   control=fit.cntl.list, covariates=rbind(temp,TP))
   }


###################################################
### code chunk number 34: fit.covar.echo (eval = FALSE)
###################################################
## model.list=list(m=4, R="unconstrained")
## kemz.temp = MARSS(dat.spp.1980, model=model.list, z.score=TRUE,
## 	form="dfa", control=cntl.list, covariates=temp)
## kemz.TP = MARSS(dat.spp.1980, model=model.list, z.score=TRUE,
## 	form="dfa", control=cntl.list, covariates=TP)
## kemz.both = MARSS(dat.spp.1980, model=model.list, z.score=TRUE,
## 	form="dfa", control=cntl.list, covariates=rbind(temp,TP))


###################################################
### code chunk number 35: save.results
###################################################
save(file="CS4--model.fits.Rdata",list=c("model.data", ls(pattern="^kemz.")))


###################################################
### code chunk number 36: covar-AICs
###################################################
print(cbind(model=c("no covars", "Temp", "TP", "Temp & TP"),
      AICc=round(c(best.fit$AICc, kemz.temp$AICc, kemz.TP$AICc,
      kemz.both$AICc))), quote=FALSE)


###################################################
### code chunk number 37: D-par (eval = FALSE)
###################################################
## D = matrix(c("D(Crypto,Temp)", "D(Diatoms,Temp)",
## 	"D(Grn,Temp)", "D(BG,Temp)", "D(Uni,Temp)",
## 	"D(Other,Temp)"), nrow=6, ncol=1)


