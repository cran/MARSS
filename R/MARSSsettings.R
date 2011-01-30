
alldefaults = list()
alldefaults$kem = alldefaults$kemsafe = list(
      inits=list(B=1, U=0, Q=0.05, Z=1, A=0, R=0.05, x0=-99, V0=1),
      model=list(Z="identity", A="scaling", R="diagonal and equal", B="identity", U="unconstrained", Q="diagonal and unequal", x0="unconstrained", V0="zero"),
      miss.value=NA,  
      control=list(minit=15, maxit=500, abstol=NULL, trace=0, kf.x0="x00", 
                   safe=FALSE, allow.degen=TRUE, min.degen.iter=50, degen.lim=1.0e-04, MCInit=FALSE, 
                   numInits = 500, numInitSteps = 10, min.iter.conv.test=15, conv.test.deltaT=9, conv.test.slope.tol= 0.5,  demean.states=FALSE,
                   boundsInits=list(B=c(-1,1), U=c(-1,1), Q = c(1,1),
                                    Z=c(0,1), A=c(-1,1), R = c(1,1) ) 
                   )
)
alldefaults$BFGS = list(
      inits=list(B=1, U=0, Q=0.05, Z=1, A=0, R=0.05, x0=-99, V0=10),
      model=list(Z="identity", A="scaling", R="diagonal and equal", B="identity", U="unconstrained", Q="diagonal and unequal", x0="unconstrained", V0="zero"),
      miss.value=NA,  
      control=list(maxit=5000, kf.x0="x00", diffuse=FALSE, MCInit=FALSE, numInits = 500, numInitSteps = 10, boundsInits=list(B=c(0,1), U=c(-1,1), Q = c(1,1), Z=c(0,1), A=c(-1,1), R = c(1,1) ),
      REPORT=NULL, reltol=NULL, fnscale=NULL, parscale=NULL, ndeps=NULL, alpha=NULL, beta=NULL, gamma=NULL, type=NULL, lmm=NULL, factr=NULL,
      pgtol=NULL, tmax=NULL, temp=NULL, lower=NULL, upper=NULL ))
alldefaults$BFGSkf = alldefaults$BFGS

model.elem =      c("A","U","B","R","Q","Z","x0","V0")
            
#method=="kemsafe"   This is the safe list
allowed=list()
allowed$kemsafe= list(
        A=c("scaling", "zero"),
        B=c("identity"),
        Q=c("unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov"),
        R=c("unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov"),
        U=c("unconstrained", "unequal", "equal", "zero"),
        x0=c("unconstrained", "unequal", "equal", "zero"),
        Z=c("identity"),
        V0=c("zero"),
        factors = c("Z"),
        matrices = c("A","B","Q","R","U","x0","Z","V0")
      )
#a fairly unconstrained list for MARSS 2.0
# unequal == unconstrained
allowed$kem = list(
        A=c("scaling", "unconstrained", "unequal", "equal", "zero"),
        B=c("identity", "zero", "unconstrained", "unequal", "diagonal and unequal", "diagonal and equal", "equalvarcov"),
        Q=c("identity", "zero", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov"),
        R=c("identity", "zero", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov"),
        U=c("unconstrained", "equal", "unequal", "zero"),
        x0=c("unconstrained", "equal", "unequal", "zero"),
        V0=c("identity", "zero", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov"),
        Z=c("identity","unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov", "onestate"),
        factors = c("Z"),
        matrices = c("A","B","Q","R","U","x0","Z","V0") 
      )

allowed$BFGS= list(
        A=c("scaling", "unconstrained", "unequal", "equal", "zero"),
        B=c("identity", "zero", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov"),
        Q=c("identity", "zero", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov"),
        R=c("identity", "zero", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov"),
        U=c("unconstrained", "equal", "unequal", "zero"),
        x0=c("unconstrained", "equal", "unequal", "zero"),
        V0=c("identity", "zero", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov"),
        Z=c("identity","unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov", "onestate"),
        factors = c("Z"),
        matrices = c("A","B","Q","R","U","x0","Z","V0") 
      )
allowed$BFGSkf=allowed$BFGS

kem.methods=c("kem","kemsafe")
optim.methods=c("BFGS","BFGSkf")      
allowed.methods = c(kem.methods, optim.methods)
