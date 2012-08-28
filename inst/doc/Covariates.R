### R code from vignette source './tex/Covariates.xRnw'

###################################################
### code chunk number 1: RUNFIRST
###################################################
require(MARSS)
options(prompt=" ", continue=" ", width=60)


###################################################
### code chunk number 2: load.plankton.data
###################################################
# Set up some data. 1972 onward
years=(1:396)[lakeWAplankton[,"Year"]>=1972]
dat = t(lakeWAplankton[years,c("Greens", "Bluegreens")])
#z.score the data
the.mean=apply(dat,1,mean,na.rm=TRUE)
the.sigma=sqrt(apply(dat,1,var,na.rm=TRUE))
dat=(dat-the.mean)*(1/the.sigma)


###################################################
### code chunk number 3: load.covar.data
###################################################
temp.offset=tp.offset=1
month=lakeWAplankton[years,"Month"]
covariates = rbind(lakeWAplankton[years-temp.offset,"Temp"],
   lakeWAplankton[years-tp.offset,"TP"],
   month, month^2, month^3)
   
#if you put the rownames on, MARSS can come up with automatic
# naming for the D and C matrices
rownames(covariates)=c("Temp","TP","mon","mon^2","mon^3")

#z.score the covariates; they are already log-transformed
the.mean=apply(covariates,1,mean,na.rm=TRUE)
the.sigma=sqrt(apply(covariates,1,var,na.rm=TRUE))
covariates=(covariates-the.mean)*(1/the.sigma)


###################################################
### code chunk number 4: covar.model.0
###################################################
Q=U=x0="zero"; B=Z="identity"
d=covariates
A="zero"
D="unconstrained"
y=dat #to show relationship between dat and the equation
model.list=list(B=B,U=U,Q=Q,Z=Z,A=A,D=D,d=d,x0=x0)
kem = MARSS(y, model=model.list)


###################################################
### code chunk number 5: covar.model.0b
###################################################
Q="unconstrained"
B="diagonal and unequal"
A=U=x0="zero"
R="diagonal and equal"
d=covariates
D="unconstrained"
y=dat #to show the relation between dat and the model equations
model.list=list(B=B,U=U,Q=Q,Z=Z,A=A,R=R,D=D,d=d,x0=x0)
kem = MARSS(y, model=model.list)


###################################################
### code chunk number 6: covar.model.1
###################################################
R=A=U="zero"; B=Z="identity"
Q="equalvarcov"
C="unconstrained"
x=dat  #to show the relation between dat and the equations
model.list=list(B=B,U=U,Q=Q,Z=Z,A=A,R=R,C=C,c=covariates)
kem = MARSS(x, model=model.list)


###################################################
### code chunk number 7: covar.model.3 (eval = FALSE)
###################################################
## matplot(t(dat))


###################################################
### code chunk number 8: covar.model.4
###################################################
model.list$B="diagonal and unequal"
kem = MARSS(dat, model=model.list)


###################################################
### code chunk number 9: covar.model.2
###################################################
x0=dat[,1,drop=FALSE]
model.list$tinitx=1
model.list$x0=x0
kem = MARSS(dat, model=model.list)


###################################################
### code chunk number 10: covar.model.5
###################################################
model.list$R=diag(0.16,2)
kem = MARSS(dat, model=model.list)


###################################################
### code chunk number 11: reset
###################################################
options(prompt="> ", continue=" +", width=120)


