
alldefaults = list()
alldefaults$kem = alldefaults$kemdev = list(
      inits=list(B=1, U=0, Q=0.05, Z=1, A=0, R=0.05, x0=-99, V0=10),
      constraint=list(Z="identity", A="scaling", R="diagonal and equal", B="identity", U="unconstrained", Q="diagonal and unequal", x0="unconstrained", V0="zero"),
      miss.value=-99,  
      control=list(minit=15, maxit=500, abstol=NULL, iter.V0=10, trace=0,
safe=FALSE, MCInit=FALSE, numInits = 500, numInitSteps = 10, min.iter.conv.test=15, conv.test.deltaT=9, conv.test.slope.tol= 0.5,
boundsInits=list(B=c(0,1), U=c(-1,1), logQ = c(log(1.0e-05),log(1.0)),
Z=c(0,1), A=c(-1,1), logR = c(log(1.0e-05),log(1.0)) ) )
)
alldefaults$BFGS = list(
      inits=list(B=1, U=0, Q=0.05, Z=1, A=0, R=0.05, x0=-99, V0=10),
      constraint=list(Z="identity", A="scaling", R="diagonal and equal", B="identity", U="unconstrained", Q="diagonal and unequal", x0="unconstrained", V0="zero"),
      miss.value=-99,  
      control=list(maxit=5000, iter.V0=10, MCInit=FALSE, numInits = 500, numInitSteps = 10,
boundsInits=list(B=c(0,1), U=c(-1,1), logQ = c(log(1.0e-05),log(1.0)),
Z=c(0,1), A=c(-1,1), logR = c(log(1.0e-05),log(1.0)) ) )
)
      
model.elem = c("A","U","B","R","Q","Z","x0")
model.elem.w.V0 = c("A","U","B","R","Q","Z","x0","V0")
            
#method=="kem"   This is the safe list
allowed=list()
allowed$kem= list(
        A=c("scaling", "zero", "use fixed/free"),
        B=c("identity"),
        Q=c("unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov", "use fixed/free"),
        R=c("unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov", "use fixed/free"),
        U=c("unconstrained", "unequal", "equal", "zero", "use fixed/free"),
        x0=c("unconstrained", "unequal", "equal", "zero", "use fixed/free"),
        Z=c("identity"),
        V0=c("zero"),
        factors = c("Q","R","U","Z","x0"),
        matrices = c("A","B","Q","R","U","x0","Z","V0")
      )
#a fairly unconstrained list for testing purposes
# unequal == unconstrained
allowed$kemdev = list(
        A=c("scaling", "unconstrained", "unequal", "equal", "zero", "ones", "use fixed/free"),
        B=c("identity", "ones", "zero", "unconstrained", "unequal", "diagonal and unequal", "diagonal and equal", "use fixed/free"),
        Q=c("unconstrained", "unequal", "diagonal and unequal", "diagonal and equal", "equalvarcov", "ones", "use fixed/free"),
        R=c("unconstrained", "unequal", "diagonal and unequal", "diagonal and equal", "equalvarcov", "ones", "use fixed/free"),
        U=c("unconstrained", "equal", "unequal", "zero", "ones", "use fixed/free"),
        x0=c("unconstrained", "equal", "unequal", "zero", "ones", "use fixed/free"),
        Z=c("identity", "ones", "use fixed/free"),
        factors = c("A","B","Q","R","U","x0","Z"),
        matrices = c("A","B","Q","R","U","x0","Z") 
      )

allowed$BFGS= list(
        A=c("scaling", "zero", "use fixed/free"),
        B=c("identity"),
        Q=c("diagonal and unequal", "diagonal and equal", "use fixed/free"),
        R=c("diagonal and unequal", "diagonal and equal", "use fixed/free"),
        U=c("unconstrained", "equal", "unequal", "zero", "use fixed/free"),
        x0=c("unconstrained", "equal", "unequal", "zero", "use fixed/free"),
        Z=c("identity"),
        factors = c("Q","R","U","Z","x0"),
        matrices = c("A","B","Q","R","U","x0","Z")
      )

kem.methods=c("kem","kemdev")
optim.methods=c("BFGS")      
allowed.methods = c(kem.methods, optim.methods)
