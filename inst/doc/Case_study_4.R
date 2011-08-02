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
# use only the data from 1977 onward
dat.spp.1977 <- lakeWAplankton[lakeWAplankton[,"Year"]>=1977,]
# create vector of phytoplankton group names
phytoplankton <- c("Cryptomonas","Diatom", "Green",
    "bluegreen","Unicells","Otheralgae")
# get only the phytoplankton
dat.spp.1977 <- dat.spp.1977[,phytoplankton]
# transpose data so time goes across columns
dat.spp.1977 <- t(dat.spp.1977)
# get number of time series
N.ts <- dim(dat.spp.1977)[1] 


###################################################
### code chunk number 3: z.score
###################################################
dat <- dat.spp.1977
Sigma <- sqrt(apply(dat,1,var,na.rm=TRUE))
Mean <- apply(dat,1,mean,na.rm=TRUE)
dat.z.scored <- (dat-Mean)*(1/Sigma)
rownames(dat.z.scored) <- rownames(dat)


###################################################
### code chunk number 4: plotdata
###################################################
spp <- rownames(dat.spp.1977)
par(mfcol=c(3,2), mar=c(3,4,1.5,0.5), oma=c(1.5,2,1,1))
for(i in spp){
  plot(dat.z.scored[i,],xlab="",ylab="Abundance index", bty="L", xaxt="n", pch=16, col="blue")
  axis(1,12*(0:dim(dat.spp.1977)[2])+1,1977+0:dim(dat.spp.1977)[2])
  title(i)
  }


###################################################
### code chunk number 5: set.up.Z
###################################################
Z.vals <- list(
"gam11",   0,      0,
"gam21","gam22",   0,
"gam31","gam32","gam33",
"gam41","gam42","gam43",
"gam51","gam52","gam53",
"gam61","gam62","gam63")
Z <- matrix(Z.vals,N.ts,3,byrow=TRUE)


###################################################
### code chunk number 6: print.Z
###################################################
print(Z)


###################################################
### code chunk number 7: set.up.A
###################################################
A.vals <- list("a1","a2","a3","a4","a5","a6")
A <- matrix(A.vals,N.ts,1)


###################################################
### code chunk number 8: set.up.QR
###################################################
Q <- B <- diag(1,3)


###################################################
### code chunk number 9: set.up
###################################################
R.vals <- list(
"r11",0,0,0,0,0,
0,"r22",0,0,0,0,
0,0,"r33",0,0,0,
0,0,0,"r44",0,0,
0,0,0,0,"r55",0,
0,0,0,0,0,"r66")

R <- matrix(R.vals,N.ts,N.ts,byrow=TRUE)


###################################################
### code chunk number 10: print
###################################################
print(R)


###################################################
### code chunk number 11: set.up.R.short
###################################################
R <- "diagonal and unequal"


###################################################
### code chunk number 12: set.up.U
###################################################
x0 <- U <- matrix(0,3,1)
x0 <- U <- "zero"


###################################################
### code chunk number 13: set.up.x0
###################################################
V0 <- diag(5,3)


###################################################
### code chunk number 14: set.model
###################################################
dfa.model <- list(Z=Z, A=A, R=R, B=B, U=U, Q=Q, x0=x0, V0=V0)


###################################################
### code chunk number 15: fit.data
###################################################
kemz.3 <- MARSS(dat.z.scored, model=dfa.model,
  control=list(demean.states=TRUE))


###################################################
### code chunk number 16: plotfits
###################################################
dat=dat.z.scored
fit=kemz.3
spp = rownames(dat)
par(mfcol=c(3,2), mar=c(3,4,1.5,0.5), oma=c(1.5,2,1,1))
for(i in 1:length(spp)){
  plot(dat[i,],xlab="",ylab="abundance index",bty="L", xaxt="n", ylim=c(-4,3), pch=16, col="blue")
  axis(1,12*(0:dim(dat)[2])+1,1977+0:dim(dat)[2])
  lines(as.vector(fit$par$Z[i,,drop=FALSE]%*%fit$states+fit$par$A[i,]))
  title(spp[i])
  }


###################################################
### code chunk number 17: set.up.two.trends
###################################################
m <- 2
n <- 6
Q <- B <- diag(1,m)
x0 <- U <- matrix(0,m,1)
V0 <- diag(5,m)
# set up Z as a list matrix
Z <- matrix(list(),n,m)
# replace all the Z values with "1" to "12" sequentially
Z[,] <- as.character(1:(m*n))
# set the correct i,j values in Z to numeric 0
for(i in 1:(m-1)) {Z[i,(i+1):m] <- 0}
dfa.model <- list(Z=Z, A=A, R=R, B=B, U=U, Q=Q, x0=x0, V0=V0)


###################################################
### code chunk number 18: fit.two.trends
###################################################
dat <- dat.z.scored
kemz.2 <- MARSS(dat,
                model=dfa.model,
                control=list(demean.states=TRUE))
print(c(kemz.3$AICc, kemz.2$AICc))


###################################################
### code chunk number 19: set.up.many.trends.no.echo
###################################################
#This is being done to speed up building the user guide
if("CS6--set.up.many.trends.Rdata" %in% dir()){
load("CS6--set.up.many.trends.Rdata")
}else{
# set up the model
levels.R <- c("diagonal and equal",
              "diagonal and unequal",
              "unconstrained")
A.vals <- list("a1","a2","a3","a4","a5","a6")
A <- matrix(A.vals,N.ts,1)
model.data <- data.frame()
dat <- dat.z.scored
n <- dim(dat)[1]
for(R in levels.R) {
    for(m in 1:(n-1)) {
        Z <- matrix(list(),n,m)
        Z[,] <- as.character(1:(m*n))
        if(m>1) {for(i in 1:(m-1)) {Z[i,(i+1):m] <- 0}}
        x0 <- U <- matrix(0,m,1)
        Q <- B <- diag(1,m)
        V0 <- diag(5,m)
        dfa.model <- list(Z=Z, A=A, R=R, B=B, U=U, Q=Q, x0=x0, V0=V0)
        kemz <- MARSS(dat,
                      model=dfa.model,
                      control=list(demean.states=TRUE))
        model.data <- rbind(model.data,
                            data.frame(R=R,
                                       m=m,
                                       logLik=kemz$logLik,
                                       K=kemz$num.params,
                                       AICc=kemz$AICc,
                                       stringsAsFactors=FALSE))
        assign(paste("kemz",m,R,sep="."), kemz)
        } # end m loop
    } # end R loop
save(file="CS6--set.up.many.trends.Rdata",list=c("model.data",ls(pattern="^kemz.")))
}


###################################################
### code chunk number 20: set.up.many.trends.echo (eval = FALSE)
###################################################
## # set up the model
## levels.R <- c("diagonal and equal",
##               "diagonal and unequal",
##               "unconstrained")
## A.vals <- list("a1","a2","a3","a4","a5","a6")
## A <- matrix(A.vals,N.ts,1)
## model.data <- data.frame()
## dat <- dat.z.scored
## n <- dim(dat)[1]
## for(R in levels.R) {
##     for(m in 1:(n-1)) {
##         Z <- matrix(list(),n,m)
##         Z[,] <- as.character(1:(m*n))
##         if(m>1) {for(i in 1:(m-1)) {Z[i,(i+1):m] <- 0} }
##         x0 <- U <- matrix(0,m,1)
##         Q <- B <- diag(1,m)
##         V0 <- diag(5,m)
##         dfa.model <- list(Z=Z, A=A, R=R, B=B, U=U, Q=Q, x0=x0, V0=V0)
##         kemz <- MARSS(dat,
##                       model=dfa.model,
##                       control=list(demean.states=TRUE))
##         model.data <- rbind(model.data,
##                             data.frame(R=R,
##                                        m=m,
##                                        logLik=kemz$logLik,
##                                        K=kemz$num.params,
##                                        AICc=kemz$AICc,
##                                        stringsAsFactors=FALSE))
##         assign(paste("kemz",m,R,sep="."), kemz)
##         } # end m loop
##     } # end R loop


###################################################
### code chunk number 21: makemodeltable
###################################################
# calculate delta-AICc
model.data$delta.AICc <- model.data$AICc - min(model.data$AICc)
# calculate Akaike weights
wt <- exp(-0.5*model.data$delta.AICc)
model.data$Ak.wt <- wt/sum(wt)
# sort results
model.tbl <- model.data[order(model.data$AICc),-4]
# calculate cumulative wts
model.tbl$Ak.wt.cum <- cumsum(model.tbl$Ak.wt)
tmpaln <- "c" #figure out the number of cols automatically
for(i in 1:ncol(model.tbl)) {tmpaln <- paste(tmpaln,"c",sep="")}
thetable <- xtable(model.tbl,caption='Model selection results.', label='ref:tablefits', align=tmpaln, digits=c(1,1,1,1,0,1,2,2))
print(thetable,type = "latex", file = paste(tabledir,"tablefit.tex",sep=""), include.rownames=FALSE,include.colnames=TRUE, caption.placement="top",table.placement="htp", sanitize.text.function = function(x){x},hline.after = c(-1,0,nrow(model.data)))


###################################################
### code chunk number 22: varimax
###################################################
best.model <- model.tbl[1,]
fitname <- paste("kemz",best.model$m,best.model$R,sep=".")
best.fit <- get(fitname)
Zrot <- varimax(best.fit$par$Z)$loadings   # Z%*%H^{-1}


###################################################
### code chunk number 23: plotfacloadings
###################################################
spp <- rownames(dat)
minZ <- 0.1
ylims <- c( -1.1*max(abs(Zrot)), 1.1*max(abs(Zrot)) )
par(mfrow=c(2,2), mar=c(4,4,0.5,0.5), oma=c(1.5,2,1,1))
for(i in 1:best.model$m) {
	plot(c(1:N.ts)[abs(Zrot[,i])>minZ], as.vector(Zrot[abs(Zrot[,i])>minZ,i]),
		type="h", lwd=2, xlab="", ylab="", xaxt="n", ylim=ylims, xlim=c(0,N.ts+1))
	for(j in 1:N.ts) {
		if(Zrot[j,i] > minZ) {text(j,-0.05,spp[j], srt=90, adj=1)}
		if(Zrot[j,i] < -minZ) {text(j,0.05,spp[j], srt=-90, adj=0)}
	abline(h=0, lwd=1, col="gray")
    } # end j loop
    mtext(paste("Factor loadings on trend",i,sep=" "),side=3,line=.5)
} # end i loop


###################################################
### code chunk number 24: plottrends
###################################################
# get ts of trends
ts.trends <- t(best.fit$states)
par(mfrow=c(ceiling(dim(ts.trends)[2]/2),2), mar=c(3,4,1.5,0.5), oma=c(1.5,2,1,1))
# loop over each trend
for(i in 1:dim(ts.trends)[2]) {
	# set up plot area
	plot(ts.trends[,i],
		ylim=c(-1.1,1.1)*max(abs(ts.trends)), 
		type="n", lwd=2,
		xlab="", ylab="", xaxt="n", yaxt="n")
	# draw zero-line
	abline(h=0, col="gray")
	# plot trend line
	par(new=TRUE)
	plot(ts.trends[,i],
		ylim=c(-1.1,1.1)*max(abs(ts.trends)), 
		type="l", lwd=2,
		xlab="", ylab="", xaxt="n")
	# add panel labels
	mtext(paste("Trend",i,sep=" "), side=3, line=0.5)
	axis(1,12*(0:dim(dat.spp.1977)[2])+1,1977+0:dim(dat.spp.1977)[2])
	} # end i loop (trends)


