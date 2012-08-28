### R code from vignette source 'Case_Study_7.xRnw'

###################################################
### code chunk number 1: RUNFIRST
###################################################
require(MARSS)
options(prompt=" ", continue=" ")
library(xtable)
tabledir="tables/"


###################################################
### code chunk number 2: load.wolf.data
###################################################
royale.dat = log(t(isleRoyal[,2:3]))


###################################################
### code chunk number 3: figwolf
###################################################
matplot(isleRoyal[,1],log(isleRoyal[,2:3]),
    ylab="log count",xlab="Year",type="l",
    lwd=3,bty="L",col="black")
legend("topright",c("Wolf","Moose"), lty=c(1,2), bty="n")


###################################################
### code chunk number 4: bad.wolf.model (eval = FALSE)
###################################################
## royale.model.0=list(B="unconstrained",Q="diagonal and unequal",
##       R="diagonal and unequal",U="unequal")
## kem.0=MARSS(royale.dat,model=royale.model.0)


###################################################
### code chunk number 5: z.score.wolf.data
###################################################
#if missing values are in the data, they should be NAs
z.royale.dat=(royale.dat-apply(royale.dat,1,mean,na.rm=TRUE))/
     sqrt(apply(royale.dat,1,var,na.rm=TRUE))


###################################################
### code chunk number 6: fit.model
###################################################
royale.model.1=list(Z="identity", B="unconstrained",
      Q="diagonal and unequal", R="diagonal and unequal",
      U="zero", tinitx=1)
cntl.list=list(allow.degen=FALSE,maxit=200)
kem.1=MARSS(z.royale.dat, model=royale.model.1, control=cntl.list)


###################################################
### code chunk number 7: fit.model.R0
###################################################
royale.model.2=list(Z="identity", B="unconstrained",
    Q="diagonal and unequal", R="zero", U="zero")
kem.2=MARSS(z.royale.dat, model=royale.model.2)


###################################################
### code chunk number 8: print.wolf.B
###################################################
wolf.B=parmat(kem.2)$B
rownames(wolf.B)=colnames(wolf.B)=rownames(royale.dat)
print(wolf.B, digits=2)


###################################################
### code chunk number 9: prep.cov.wolf.moose
###################################################
clim.dat= t(isleRoyal[1:52,c(4,10,6)])
z.score.clim.dat=(clim.dat-apply(clim.dat,1,mean,na.rm=TRUE))/
     sqrt(apply(clim.dat,1,var,na.rm=TRUE))


###################################################
### code chunk number 10: cov.wolf.moose.model
###################################################
royale.model.3=list(Z="identity", B="unconstrained",
    Q="diagonal and unequal", R="zero", U="zero",
    C=matrix(list(0,"Moose win temp",0,"Moose win precip",
         0,"Moose sum temp"),2,3),
    c=z.score.clim.dat)


###################################################
### code chunk number 11: cov.wolf.moose.model
###################################################
kem.3=MARSS(z.royale.dat[,2:53], model=royale.model.3)


###################################################
### code chunk number 12: figwolfcov
###################################################
cor.fun=function(x, y){text(0.5,0.5,format(cor(x,y),digits=2),cex=2)}
pairs(t(z.score.clim.dat),lower.panel=cor.fun)


###################################################
### code chunk number 13: bad.data.2 (eval = FALSE)
###################################################
## bad.data=z.royale.dat+matrix(rnorm(100,0,sqrt(.2)),2,50)
## kem.bad=MARSS(bad.data, model=model)


###################################################
### code chunk number 14: load.plankton.data
###################################################
# only use the plankton, daphnia, & non-daphnia
plank.spp = c("Large Phyto","Small Phyto","Daphnia","Non-daphnia")
plank.dat = ivesDataByWeek[,plank.spp] 
#The data are not logged
plank.dat = log(plank.dat)  
#Transpose to get time going across the columns
plank.dat = t(plank.dat)
#make a demeaned version
d.plank.dat = (plank.dat-apply(plank.dat,1,mean,na.rm=TRUE))


###################################################
### code chunk number 15: figplank
###################################################
matplot((1:(52*6))[27:295],t(d.plank.dat),type="l",lty=c(1,1,1,1),lwd=c(1,1,3,3),xlab="week of study",ylab="log biomass",xaxt="n",xlim=c(11,52*6-11),bty="L")
#axis(1,at=(1:(52*6))[seq(27,295,2)])
axis(1,at=seq(1,52*6,2))
abline(v=c(52*(1:6)))
abline(h=0)


###################################################
### code chunk number 16: plankton.model
###################################################
Z="identity"
U="zero"
B="unconstrained"
Q=matrix(list(0),4,4); diag(Q)=c("Phyto","Phyto","Zoo","Zoo")
R=matrix(list(0),4,4); diag(R)=c("Phyto","Phyto","Zoo","Zoo")
plank.model.0=list(Z=Z, U=U, Q=Q, R=R, B=B)


###################################################
### code chunk number 17: fit.plank.model.0
###################################################
plank.model.0$tinitx=1
kem.plank.0 = MARSS(d.plank.dat, model=plank.model.0 )


###################################################
### code chunk number 18: print.B.0
###################################################
#Cleaning up the B matrix for printing
B.0 = parmat(kem.plank.0)$B[1:4,1:4]
rownames(B.0) = colnames(B.0) = c("LP","SP","D","ND")
print(B.0,digits=2)


###################################################
### code chunk number 19: print.B.Ives
###################################################
#Cleaning up the B matrix for printing
B.Ives.ML = matrix(c(.5,NA,NA,NA,-.39,.076,NA,.1,NA,-.02,.77,NA,NA,-.1,NA,.55),4,4)
B.Ives.Obs = matrix(c(.48,NA,NA,NA,-.39,.25,NA,.1,NA,-.17,.74,0,NA,-.11,0,.6),4,4)
B.Ives=B.Ives.Obs
rownames(B.Ives) = colnames(B.Ives) = c("LP","SP","D","ND")
print(B.Ives,digits=2,na.print="--")


###################################################
### code chunk number 20: test.rm.NAs (eval = FALSE)
###################################################
## test.dat=d.plank.dat[,!is.na(d.plank.dat[1,])]
## test = MARSS(test.dat, model=plank.model.0 )


###################################################
### code chunk number 21: fit.plank.model.1
###################################################
plank.model.1=plank.model.0
plank.model.1$Q="unconstrained"
kem.plank.1 = MARSS(d.plank.dat, model=plank.model.1)


###################################################
### code chunk number 22: print.B.1
###################################################
#Cleaning up the B matrix for printing
B = parmat(kem.plank.1)$B[1:4,1:4]
rownames(B) = colnames(B) = c("LP","SP","D","ND")
B[B==0]=NA
B.1=B
print(B,digits=2,na.print="--")


###################################################
### code chunk number 23: B.2
###################################################
B.2=matrix(list(0),4,4) #set up the list matrix
diag(B.2)=c("B11","B22","B33","B44") #give names to diagonals
#and names to the estimated non-diagonals
B.2[1,2]="B12"; B.2[2,3]="B23"; B.2[2,4]="B24"; B.2[4,2]="B42" 
print(B.2)


###################################################
### code chunk number 24: fit.plank.model.2
###################################################
#model 2
plank.model.2=plank.model.1
plank.model.2$B = B.2
kem.plank.2= MARSS(d.plank.dat, model=plank.model.2)


###################################################
### code chunk number 25: print.B.2
###################################################
#Cleaning up the B matrix for printing
B = parmat(kem.plank.2)$B[1:4,1:4]
rownames(B) = colnames(B) = c("LP","SP","D","ND")
B[B==0]=NA
B.2=B
print(B,digits=2,na.print="--")


###################################################
### code chunk number 26: fit.plank.model.3
###################################################
#model 3
plank.model.3=plank.model.2
plank.model.3$R=diag(c(.04,.04,.16,.16))
kem.plank.3= MARSS(d.plank.dat, model=plank.model.3)


###################################################
### code chunk number 27: prep.covariates
###################################################
#transpose to make time go across columns
#drop=FALSE so that R doesn't change our matrix to a vector
phos = t(log(ivesDataByWeek[,"Phosph",drop=FALSE]))
d.phos = (phos-apply(phos,1,mean,na.rm=TRUE))


###################################################
### code chunk number 28: add.covar.model.3
###################################################
plank.model.4=plank.model.3
plank.model.4$C=matrix(list("C11","C21",0,0),4,1)
plank.model.4$c=d.phos


###################################################
### code chunk number 29: plank.model.4
###################################################
kem.plank.4= MARSS(d.plank.dat, model=plank.model.4)


###################################################
### code chunk number 30: add.fish.to.data
###################################################
#transpose to make time go across columns
#drop=FALSE so that R doesn't change our matrix to a vector
fish = t(log(ivesDataByWeek[,"Fish biomass",drop=FALSE]))
d.fish = (fish-apply(fish,1,mean,na.rm=TRUE))
#plank.dat.w.fish = rbind(plank.dat,fish)
d.plank.dat.w.fish = rbind(d.plank.dat,d.fish)


###################################################
### code chunk number 31: B.covar
###################################################
B=matrix(list(0),5,5)
diag(B)=list("B11","B22","B33","B44","Bfish")
B[1,2]="B12";B[2,3]="B23"; B[2,4]="B24"
B[4,2]="B42"; 
B[1:4,5]=list(0,0,"C32","C42")
print(B)


###################################################
### code chunk number 32: C.covar
###################################################
C=matrix(list("C11","C21",0,0,0),5,1)


###################################################
### code chunk number 33: R.covar
###################################################
R=matrix(list(0),5,5)
diag(R)=list(0.04,0.04,0.16,0.16,0.36)


###################################################
### code chunk number 34: Q.covar
###################################################
Q=matrix(list(0),5,5); 
Q[1:4,1:4]=paste(rep(1:4,times=4),rep(1:4,each=4),sep="")
Q[5,5]="fish"  
Q[lower.tri(Q)]=t(Q)[lower.tri(Q)]
print(Q)


###################################################
### code chunk number 35: fit.covar.model
###################################################
plank.model.5=plank.model.4
plank.model.5$B=B
plank.model.5$C=C
plank.model.5$Q=Q
plank.model.5$R=R
kem.plank.5=MARSS(d.plank.dat.w.fish, model=plank.model.5)


###################################################
### code chunk number 36: print.B
###################################################
#Cleaning up the B matrix for printing
B.5 = parmat(kem.plank.5)$B[1:4,1:4]
rownames(B.5) = colnames(B.5) = c("LP","SP","D","ND")
B.5[B.5==0]=NA
print(B.5,digits=2,na.print="--")


###################################################
### code chunk number 37: makemodeltable
###################################################
B.names=c("B11","B22","B33", "B44", "B12", "B23", "B24", "B42")
#rename kem.plank.0 and kem.plank.1 B to make keeping track of params easier
rownames(kem.plank.0$par$B)=paste("B",rep(1:4,times=4),rep(1:4,each=4),sep="")
rownames(kem.plank.1$par$B)=paste("B",rep(1:4,times=4),rep(1:4,each=4),sep="")
names.ests=c("B11","B22","B33", "B44", "B12", "B23", "B24", "B42", "C11", "C21", "C32","C42")
Ives.ests.Obs = c(.48,.25,.74,.6,-.39,-.17,-.11,.1, .25,.25,-.14,-.045)
Ives.ests.ML=c(0.5,.076,.77,.55,-.39,-.02,-.1, .1,.2,.32,-.13,-.048)
model.data = cbind(
Ives.ests.Obs,
c(kem.plank.0$par$B[B.names,],NA,NA,NA,NA),
c(kem.plank.1$par$B[B.names,],NA,NA,NA,NA),
c(kem.plank.2$par$B[B.names,],NA,NA,NA,NA),
c(kem.plank.3$par$B[B.names,],NA,NA,NA,NA),
c(kem.plank.4$par$B[B.names,],kem.plank.4$par$U[c("C11","C21"),],NA,NA ),
c(kem.plank.5$par$B[B.names,],kem.plank.5$par$U[c("C11","C21"),],kem.plank.5$par$B[c("C32","C42"),] )
)
rownames(model.data)=names.ests
colnames(model.data)=c("Ives et al.","Model 0","Model 1","Model 2","Model 3","Model 4", "Model 5")
tmpaln="c" #figure out the number of cols automatically
for(i in 1:ncol(model.data)) tmpaln = paste(tmpaln,"c",sep="")
thetable = xtable(model.data, caption='The parameter estimates under the different plankton models.  Models 0 to 3 do not include covariates, so the C elements are blank.  Bij is the effect of species $i$ on species $j$. 1=large phytoplankton, 2=small phytoplankton, 3=Daphnia, 4=non-Daphnia zooplankton. The Ives et al. (2003) estimates are from their table 2 for the low planktivory lake with the observation model.', label='ref:tableplank', align=tmpaln, digits=2)
print(thetable,type = "latex", file = paste(tabledir,"tableplank.tex",sep=""), include.rownames=TRUE,include.colnames=TRUE, caption.placement="top",table.placement="htp", sanitize.text.function = function(x){x},hline.after = c(-1,0,nrow(model.data)))


###################################################
### code chunk number 38: logLik.variates
###################################################
tmp=kem.plank.5
tmp$model$data[5,]=NA
LL.variates=MARSSkf(tmp)$logLik


###################################################
### code chunk number 39: Reset
###################################################
options(prompt="> ", continue="+ ")
options(width=120)


