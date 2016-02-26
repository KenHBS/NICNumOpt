###############################################################################
##################    Numerical Optimization    ###############################
###############################################################################


# ------------------------------------------------------------------------------
# Three self-implemented algorithms for numerical optimization. 
# ------------------------------------------------------------------------------
# Usage: Optimize a function with global optima
# ------------------------------------------------------------------------------
# Keywords: Nelder-Mead, Newton-Raphson, Gradient-Descent
# ------------------------------------------------------------------------------ 
# Author: Michael Lebacher, Ken Schröder, Johannes Stoiber, 2015/11/18
# ------------------------------------------------------------------------------



######################################################################
#################         Structure          #########################
######################################################################

## Loading required packages
## Data import and preparation:
##      - Set working directory
##      - Importing data
##      - Extract variables & scaling
##      - Data matrix X and output vector y
## Estimation with R-function GLM
## Basic helper functions:
##      - logit operator
##      - (negative) log likelihood function
##      - Score
##      - Hessians
## Estimation with R-package optim

##############    SELF WRITTEN ALGORITHMS     ###########
## Newton-Raphson algorithm:
##      - Algorithm code
##      - Time taken
##      - Visualizing parameter estimate convergence
## Nelder-Mead algorithm:
##      - Algorithm code
##      - Time taken
##      - Visualizing parameter estimate convergence
## Gradient-Descent (fixed learning parameter):
##      - Algorithm code
##      - Time taken
##      - Visualizing parameter estimate convergence
## Gradient-Descent (flexible learning parameter):
##      - Algorithm code
##      - Time taken
##      - Visualizing parameter estimate convergence

##############     ALGORITHM COMPARISON      ############
## Number of iterations until convergence:
## Sensitivity of accuracy w.r.t. computation time

##############   FUNCTION WITH MULTIPLE      ############
##############       LOCAL OPTIMA            ############
## Helper functions:
##      - f(x,y) = 2x^2 - 4xy + y^4
##      - Score
##      - Hessian
## Adapted algorithms for f(x,y):
##      - Adapted Newton-Raphson
##      - Adapted Gradient-Descent, fix
##      - Adapted Gradient-Descent, flex
##      - Adapted Nelder-Mead

########################################################
###########     Loading required packages     ##########
########################################################
library(foreign)
library(stats)

########################################################
#########     Data import and preparations    ##########
########################################################

## Set working directory:
  # Johannes 
  setwd ( "C:/Users/Johannes/Dropbox/Numerical Introductory" )
  # Ken
  setwd("/Users/Ken/Dropbox/Numerical Introductory")
  # Michael
  setwd("/Users/apple/Dropbox/Numerical Introductory")

## Importing data (mroz.dta):
  # Johannes
    data <- read.dta("mroz.dta")
  # Ken
    data<-read.dta("/Users/Ken/Dropbox/Numerical Introductory/mroz.dta")
  # Michael
    data<-read.dta("/Users/apple/Dropbox/Numerical Introductory/mroz.dta")

## Extract variables & scaling:
    head(data)
  
  # define variables from the data set
    inlf <- data$inlf
    educ <- data$educ
    huswage <- data$huswage

  # Scaling (improves the numerical tractability):
    inlf <- (inlf - min(inlf))/(max(inlf)-min(inlf))
    educ <- (educ - min(educ))/(max(educ)-min(educ))
    huswage <- (huswage - min(huswage))/(max(huswage)-min(huswage))

  # Define the x vector (Nxk)
    X <- cbind(matrix(seq(1, 1, length.out = 753), ncol =1 ),
               matrix(educ, ncol = 1), matrix(huswage, ncol =1 ))

  # Define the y vector (Nx1)
    y <- matrix(inlf,ncol=1)
    
########################################################
###########   Estimation with R-function GLM  ##########
########################################################
    
solution_mle <- glm(formula= inlf~educ + huswage,
                    family = binomial(link = "logit"))
summary(solution_mle)

########################################################
###########     Basic helper functions      ############
########################################################

# Logit operator:
  logit <- function(X, beta) {
   return(exp(X %*% beta)/(1+ exp(X %*% beta)) )
  }

# (Negative) Log likelihood function:
  # stable parametrisation of the log-likelihood function
  # Note: The negative of the log-likelihood is being returned, 
  # since we will minimize the function
  ll <- function(beta, X, y) {
    return(-sum(y*(X %*% beta - log(1+exp(X %*% beta)))
      + (1-y)*(-log(1 + exp(X %*% beta)))
      ) 
    ) 
  }

# Score of the log likelihood
  score <- function(beta, X, y){
    s <- t(X) %*% (y - logit(X,beta)) 
    return(s)
  }

# Hessian of the log likelihood
  # Hessian1:
  hessian <- function(beta, X){
    h<-0
    for (i in dim(X)[2]){
      h<-h+ logit(X,beta)[i]*(1-logit(X,beta)[i])*X[i,]%*%t(X[i,])
    }
    return(h)
  }

  # Hessian2:
  hessian2 <- function(beta, X, y){
    p <- logit(X, beta)
    diag <- diag(diag(p %*% t(1-p)))    # the first diag extracts the diagonal elements into vector. the second puts the vector on the diagonal of a zero matrix.
    h <- -t(X) %*% diag %*% X           # Hessian is equal to -X' W X , where W is diag(p(1-p))
    return(h)
  }

  
########################################################
###########     Solution via optim package      ########
########################################################

# Initial set of parameters:
  beta_start <- runif(dim(X)[2],0,1)

# Solution via optim:
  solution <-optim(beta_start, ll, beta, X,y )
  solution

  solution_nm <- optim(beta_start, ll, beta, X, y, method="Nelder-Mead")
  solution_nm

# Time taken for computation:
  start_time <- Sys.time()
  optim(beta_start, ll, beta, X,y )
  end_time <- Sys.time()
  time_taken_optim <- end_time - start_time
  time_taken_optim

# Compare parameter to the solution of the glm package
  solution_mle
  
  
################################   SELF WRITTEN ALGORITHMS   ###########################  

########################################################
############         Newton-Raphson         ############
########################################################

# Newton-Raphson algorithm code:
  newton_raphson <- function(beta,           # k x 1 matrix of starting values for beta
                            X,               # n x k matrix of regressors
                            y,               # n x 1 vector of dependent variables
                            tol = 0.0001     # desired tolerance level
    ){
    update_size <- 1                                                  # some initial value for the update size
    res <- c()                                                        # create a container for the results of the iteration steps
    while (update_size > tol){                                        # iteration step, until the improvement is smaller than tolerance.
      updater <- solve(hessian2(beta, X, y)) %*% score(beta, X, y)    # The inverse hessian multiplied by the score for some X, y and beta.     
      update_size <- sqrt(sum((updater)^2))                           # size of the update in absolute value
      beta <- beta - updater                                          # updating the old parameter value
      res <- matrix(c(res, beta), nrow=dim(X)[2])                     # saving the iteration values at every iteration.
    }
    return(res)
  }

# Time taken for computation:
  start <- Sys.time()
  newton <- newton_raphson(beta_start, X, y, tol = 0.0001)             # run an optimization and store the results
  time.nr <- Sys.time() - start
  time.nr

# visualize the convergence of the parameters:
  pdf("Convergence_newton_raphson_const.pdf")
  plot(newton[1,], xlab = "Iteration", ylab = "Parameter value", main = "Convergence of the constant (Newton-Raphson), tolerance = 1e^-7)", type = "l")
  dev.off()

  pdf("Convergence_newton_raphson_educ.pdf")
  plot(newton[2,], xlab = "Iteration", ylab = "Parameter value", main = "Convergence of the 'educ' parameter (Newton-Raphson), tolerance = 1e^-7)", type = "l")
  dev.off()

  pdf("Convergence_newton_raphson_huswage.pdf")
  plot(newton[3,], xlab = "Iteration", ylab = "Parameter value", main = "Convergence of the 'huswage' parameter (Newton-Raphson), tolerance = 1e^-7)", type = "l")
  dev.off()


  
########################################################
############         Nelder-Mead        ################
########################################################
  
# Nelder-Mead algorithm code:
  nelder_mead <- function(ll,             # function to optimize        
                          X,              # nxk matrix of regressors
                          y,              # nx1 vector of dependend variable
                          tol = 0.0001,   # desired tolerance level
                          alpha = 1,      # standard alpha for Nelder-Mead
                          beta = 1/2,     # standard beta for Nelder-Mead
                          gamma = 2,      # standard gamma for Nelder-Mead
                          xmin = -3,      # min of interval where random number is drawn for initial value of parametervector
                          xmax = 3,       # max of interval where random number is drawn for initial value of parametervector
                          tau1 = 0.05,    # usual tau to compute initial simplex
                          tau0 = 0.00025  # tau to compute initial simplex if if value in the paramtervector on i-th place is zero
  ) {
                      
  # Define the loglikelihood as a function that solely depends on beta                      
    f <- function(beta){
      return(ll(beta, X, y))
    }

  # Define a variable n that gives informations about the number of parameters to estimate
    n <- dim(X)[2]

  # Draw random numbers for the initail value
    x_0 <- runif(seq(1:n), min = xmin, max = xmax)

  # construct matrix NM, where starting values for the initial simplex are stored, first columns is x_0,
  # all the other values are set to zero and later replace
    NM = matrix(c(x_0, rep(0,n*n)), nrow = n)  

  # construct identity matrix with dim n, where each column is used as unity vector
    I = diag(n)

  # Compute initial values for the Simplex
    for (i in 1:n){
      if(NM[i,1] != 0) {NM[,1+i] = NM[,1] + tau1*I[,i]}  else {NM[,1+i] = NM[,1] + tau0*I[,i]}
    }

    check <- tol+1            # define a variable check which is used as termination criterion and higher than tolerance level in the initial step and adjusted with every iteration
    res <- c()                # construct a vector res, where results from optimization procedure are stored at every iteration 
    while(check>tol^2) {      # start Nelder-Mead procedure and iterate as long as tolerance level is bigger than termination criterion
      results <- apply(NM,2,f)    # define a variable results to identify best, worst and second worst evaluation of the function
      a = order(results)          # command order gives the information about the ascending order of the elements in results

      # find points x_b, x_w, x_sw
      b = NM[,a[1]]               # best: lowest value of function
      sw = NM[,a[n]]              # second worst: second highest value of function
      w = NM[,a[n+1]]             # worst: highest value of function
  
      # find centroid of the n best points
      z = 1/n*(apply(NM,1,sum)-w) # compute centroid
      r = z + alpha*(z-w)         # compute reflection point r
  
      # Reflection if if object function evaluated at r is between the evaluation of the function at b and sw
      if (f(b) <= f(r) & f(r)<= f(sw)) {
        NM[,a[n+1]]<-r            # If condition fulfilled, replace w with r
      }
  
      # Expansion if object function evaluated r is smaller than b
      if (f(r) < f(b)) {  
        e = z + beta*(z-w)      # compute expansion point e
        if (f(e)<f(r)) {        
          NM[,a[n+1]]<-e        # if f(e)<f(r) replace w with e
        } 
        else {
          NM[,a[n+1]]<-r        # replace w with r
        } 
      }
  
      # Contraction if object function evaluated at r is bigger than evaluation at sw
      if (f(r) > f(sw)) {
        if (f(r) > f(w)) {
          k = z + gamma*(w-z)   # compute contraction point k
        }
        else {
          k = z + gamma*(z-w)   # use different method to compute contraction point k
        }
        if (f(k) < f(w)) {      # compare object function evaluated at the contraction point with object function evaluated at w
          NM[,a[n+1]]<-k        # replace w with r
        }
        else {
          NM <- (NM+b)/2        # contract simplex by b
        }
      }
  
      fdash <- 1/(n+1)*sum(apply(NM,2,f))
      check<- 1/(n+1)*sum(  ( (apply(NM,2,f)-fdash)^2)  )          # compute a termination criteria
      res<-matrix(c(res,NM[,order(apply(NM,2,f))[1]]),nrow=n)      # store in res the parameters, where object function is lowest 
    } # end of iteration step

    return(res) # return matrix that shows development of the parameter vector
  } # end of function nelder_mead()

# Time taken for computation:
  start <- Sys.time()
  sol_nm <- nelder_mead(ll, X, y, tol = 0.0001)
  time.nm <- Sys.time() - start
  time.nm

# Visualizing parameter estimate convergence:
  # First extract the relevant estimates:
    sol_nm<-nelder_mead(ll,X,y, tol=0.0001)   # run an optimization problem and store results
    dim(sol_nm)                   # use dimensions of sol_nm
    paste( "number of iterations needed is", dim(sol_nm)[2])
    # access constant from last iteration of nm_procedure
    beta0_nm <- sol_nm[1,dim(sol_nm)[2]]
    beta1_nm <- sol_nm[2,dim(sol_nm)[2]]
    beta2_nm <- sol_nm[3,dim(sol_nm)[2]]
    beta0_nm
  # Make plot to show convergence of the parameter
    pdf("Convergence_constant_nelder_mead.pdf")
    plot(sol_nm[1,],xlab="Iteration",ylab="Parameter value",main="Convergence of the constant (Nelder-Mead)",type="l")
    dev.off()

    pdf("Convergence_educ_nelder_mead.pdf")
    plot(sol_nm[2,],xlab="Iteration",ylab="Parameter value",main="Convergence of the 'educ' parameter (Nelder-Mead)",type="l")
    dev.off()

    pdf("Convergence_huswage_nelder_mead.pdf")
    plot(sol_nm[3,],xlab="Iteration",ylab="Parameter value",main="Convergence of the 'huswage' parameter (Nelder-Mead)",type="l")
    dev.off()

    
########################################################
############   Gradient-Descent, fix    ################
########################################################

# Gradient-Descent, fix algorithm code:
  gd_alpha_fix<- function(score,                         # score function that should be evaluated
                          X,                             # Nxk matrix of regressors
                          y,                             # Nx1 matrix of depend variable
                          beta_start=runif(dim(X)[2],-1,1),      # Starting value for beta, by default uniformly distributed on U(-1,1)
                          tol=0.0001){                   # Termination tolerance
                  
    diff_beta <- 1                        # starting value for the differences between beta_old and beta_new
    alpha <- 1/dim(X)[1]                  # fixed value for the tuning parameter alpha according to the dimension
    res <- c()                          # create a container for results
    while(diff_beta>tol) {                               # the while loop stops if the squared difference between the values of b_old and b_new falls below the threshold defined by tol
      diff_beta <- sqrt(sum((beta_start - (beta_start + alpha*score(beta_start,X,y)) )^2))
      beta_start<- beta_start + alpha*score(beta_start,X,y) # Iteration step
      res <- matrix(c(res, beta_start),nrow=dim(X)[2])      # fills the updated results in the results matrix
    }
    return(res)                                         # returns a matrix with the solution parameters for each (!) iteration
  }


# Time taken for computation:
  start_time <- Sys.time()                     
  sol_fix<-gd_alpha_fix(score,X,y)                 # Save the result of the optimation
  end_time <- Sys.time()
  time_taken_gd_fix <- end_time - start_time
  time_taken_gd_fix                                # Evaluation of the time needed
  iterations_gd_fix <- dim(sol_fix)[2]             # Dimension of the result gives the number of iterations

# Visualizing parameter estimate convergence:
  pdf("Convergence_constant_gd_alpha_fix.pdf")
  plot(sol_fix[1,],xlab="Iteration",ylab="Parameter value",main="Convergence of the constant (Gradient-Descent)",type="l")
  dev.off()

  pdf("Convergence_educ_gd_alpha_fix.pdf")
  plot(sol_fix[2,],xlab="Iteration",ylab="Parameter value",main="Convergence of educ (Gradient-Descent)",type="l")
  dev.off()

  pdf("Convergence_huswage_gd_alpha_fix.pdf")
  plot(sol_fix[3,],xlab="Iteration",ylab="Parameter value",main="Convergence of huswage  (Gradient-Descent)",type="l")
  dev.off()

    
########################################################
#############  Gradient-Descent, flex   ################
########################################################

# Gradient-Descent, flex algorithm code:
  # Gradient-Descent use alpha as a flexible learning parameter 
  # depending on the changes in the error rate
  gd_alpha_flex<- function(score,                    # score function that sould be evaluated
                           X,                        # Nxk matrix of regressors
                           y,                        # Nx1 matrix of depend variable
                           beta_start=runif(dim(X)[2],-1,1),      # Starting value for beta, by default uniformly distributed on U(-1,1)
                           alpha=1,                  # basic value for the learning rate alpha
                           tol=0.0001){             # Termination tolerance
  
    diff_beta<-1                                     # starting value for the differences between beta_old and beta_new
    res <- c()                                       # creates a container for the results
    while(diff_beta>tol){
      diff_beta <- sqrt(sum((beta_start - (beta_start + alpha*score(beta_start,X,y)) )^2))  # evaluates the difference between beta_old and beta_new
      Error_1<-  sum( (y - logit(X,beta_start))^2 ) # Estimates the Error_Rate_old as squared residual sum
      beta_start<- beta_start + alpha*score(beta_start,X,y) # iteration step beta_new=f(beta_old)
      Error_2<-  sum( (y - logit(X,beta_start))^2 ) # Estimates the Error_Rate_new as squared residual sum
      
      # Update of the learning rate
        if ( (Error_1-Error_2)<0 ) {alpha<-alpha*0.95}
        else {alpha<-alpha}
      res <- matrix(c(res, beta_start),nrow=dim(X)[2])
    }
    return(res) 
  }

    
# Time taken for computation:
  start_time <- Sys.time()                     
  sol_flex<-gd_alpha_flex(score,X,y)                # Save the result of the optimization
  end_time <- Sys.time()
  time_taken_gd_flew <- end_time - start_time
  time_taken_gd_flew                                # Evaluation of the time needed

  iterations_gd_flex<- dim(sol_flex)[2]             # Dimension of the result gives the number of iterations

  
# Visualizing parameter estimate convergence:
  pdf("Convergence_constant_gd_alpha_flex_.pdf")
  plot(sol_flex[1,],xlab="Iteration",ylab="Parameter value",
       main="Convergence of the constant (Gradient-Descent, flexible learning)",
       type="l")
  dev.off()

  pdf("Convergence_educ_gd_alpha_flex_.pdf")
  plot(sol_flex[2,],xlab="Iteration",ylab="Parameter value",
       main="Convergence of the slope (Gradient-Descent, flexible learning",
       type="l")
  dev.off()

################################     ALGORITHM  COMPARISON     ##################################  

  
########################################################
##########   Comparing nr. of iterations   #############
########################################################
    
## Comparing the number of iterations needed for convergence:
  # save results about parameters and number of iterations in specific variables
    nm <-       c(((sol_nm)[,dim(sol_nm)[2]]), dim(sol_nm)[2], ll((sol_nm)[,dim(sol_nm)[2]],X,y))
    gd_flex <-  c(((sol_flex)[,dim(sol_flex)[2]]), dim(sol_flex)[2], ll((sol_flex)[,dim(sol_flex)[2]],X,y))     
    gd_fix <-   c(((sol_fix)[,dim(sol_fix)[2]]), dim(sol_fix)[2], ll((sol_fix)[,dim(sol_fix)[2]],X,y))
    nr <-       c(((newton)[,dim(newton)[2]]), dim(newton)[2], ll((newton)[,dim(newton)[2]],X,y))
    opt_nm <-   c(solution_nm$par, solution_nm$counts[1], ll(solution_nm$par,X,y))

  # Merge that information in a matrix
    Res_Mat <- as.matrix(cbind(opt_nm, nm,  nr, gd_fix, gd_flex))
    rownames(Res_Mat) <- c("beta_0", "beta_1","beta_2", "iterations", "likelihood")
    Res_Mat <- t(Res_Mat)
    Res_Mat

  # Provide a table for LaTeX, with the needed information
    library("stargazer") 
    stargazer(Res_Mat, title="Model Comparison", covariate.labels=c("","\beta_{0}" , "$\beta_{1}$", "Iterations", "log-likelihood"))
    
    
########################################################
#########   Comparing accuracy sensitivity   ###########
########################################################
    
## Comparing the calculation time at different levels of accuracy:
  # a function that returns the average time of 15 loops for different tollerance levels.
    accuracy.sensitivity <- function(acc_range=5,
                                     method = c("nr",         # "nr" for Newton-Raphson, 
                                                "nm",         # "nm" for Nelder-Mead, 
                                                "gd.fix",     # "gd.fix" for fix Gradient-Descent, 
                                                "gd.flex",    # "gd.flex" for flex Gradient-Descent   
                                                "optim"))
    { avg_time <- matrix(rep(1), 
                         ncol = length(method), 
                         nrow = acc_range)                    # container for the average time every tolerance takes to compute
      rowname <- 1/(10^(1:acc_range))
      for (k in 1:length(method)){
        methods <- method[k]
        for (j in 1:acc_range){                                # determines up to how many decimal points
          toll <- (1 / (10^j))
          record <- c()                                        # container for the time every single loop takes to compute
          for (i in 1:10){                                     # the number of times the function is clocked. The speeds will be averaged.
            start_time <- Sys.time()
            if (methods == "nr"){
              newton_raphson(beta_start, X, y, tol = toll)     # runs the newton_raphson function for different tolerance levels.
            }
            if (methods == "nm"){
              nelder_mead(ll, X, y, tol = toll)
            }
            if (methods == "gd.fix"){
              gd_alpha_fix(score, X, y)
            }
            if (methods == "gd.flex"){
              gd_alpha_flex(score, X, y)
            }
            if (methods == "optim"){
              optim(beta_start, ll, beta, X, y, method="Nelder-Mead", control = list(reltol = toll))
            }
            time_taken <- Sys.time() - start_time
            record <- c(record, time_taken) 
          }
          avg_time[j, k] <- mean(record)
          colnames(avg_time) <- method
          rownames(avg_time) <- rowname
        }
      }
      return(avg_time)
    } # end of sensitivity function
    
  # Saving the accuracy sensitivity    
    tol.sensitivity <- accuracy.sensitivity()
    
  # Exporting to LaTeX:
    library("stargazer") 
    stargazer(tol.sensitivity, title="Accuracy Sensitivity", 
              covariate.labels=c("Accuracy","Newton-Raphson", "Nelder-Mead", 
                                 "Gradient-Descent, fix", 
                                 "Gradient-Descent, flex", "\texttt{optim}"))
    
    
#########################    FUNCTION WITH MULTIPL LOCAL MAXIMA   ###################
    
########################################################
##########      Basic Helper Functions     #############
########################################################

# We have to adapt the optimization algorithms to optimize other problems than the logit-model.

# Target function f(x,y) = 2x^2 - 4xy + y^4 + 2
  funct_mult <- function(X){              # X contains x_1 and x_2
      x_1 <- X[1]
      x_2 <- X[2]
      return(2*x_1^2 - 4*x_1*x_2 + x_2^4 + 2)
      }

# Score function of f(x,y):
  score_mult <- function(X){              # X contains x_1 and x_2
    # The original function is 2x^2 - 4xy + y^4 + 2
    x_1 <- X[1]
    x_2 <- X[2]
    derivative1 <- 4*x_1 - 4*x_2
    derivative2 <- -4*x_1 + 4*x_2^3
    return(c(derivative1, derivative2))
  }

# Hessian of f(x,y):
  hessian_mult <- function(X){
    x_1 <- X[1]
    x_2 <- X[2]
    derivative11 <- 4
    derivative12 <- -4
    derivative21 <- -4
    derivative22 <- 12*x_2^2
    hessian <- rbind(c(derivative11, derivative12), c(derivative21, derivative22))
    return(hessian)
  }

  
########################################################
##########       Adapted Algorithms         ############
########################################################
  
# Adapted Newton-Raphson:
    # Output is a table containing values for every iteration. The values
    # included are: 
      # x1 and x2
      # Objective function evaluated right BEFORE x1, x2
      # Objective function evaluated at x1, x2
      # Objective function evaluated right AFTER x1, x2  (to see whether we have a local min, max or saddle point)
    NR_mult <- function(X,                   # X is c(x1, x2) 
                        tol = 0.00001
                        ){
      res <- c(X, funct_mult(X-c(0.1,0.1)), funct_mult(X), funct_mult(X+c(0.1, 0.1)))
      update_size <- 1
      while(update_size > tol){
        updater <- solve(hessian_mult(X)) %*% score_mult(X)
        update_size <- sqrt(sum((updater)^2))
        X <- X - updater
        subres <- c(X, funct_mult(X - c(0.1, 0.1)), funct_mult(X), funct_mult(X + c(0.1, 0.1) ))
        res <- matrix(rbind(res, subres), 
                      ncol=5, 
                      dimnames = list(c(), c("x_1", "x_2", "f(X-0.1)", "f(X)", "f(X+0.1)")))
      }
      return(res)
    }

    NR_mult(c(19,-4))
    NR_mult(c(5,5))
    NR_mult(c(0,0))        # stays at the local maximum of 2
    NR_mult(c(0,0.1))      # also fails to get away from the local maximum..

  
  
# Adapted Gradient-Descent, fix
  # WARNING: This function keeps on making overshoots and eventually shoots to -inf and +inf, except for value in the (-2,2) range
  # NOTE: The convergence is especially affected by high/low values for x_2. The function is much more robust for values for x_1
  # Output is a table containing values for every iteration. The values
  # included are: 
  # x1 and x2
  # Objective function evaluated right BEFORE x1, x2
  # Objective function evaluated at x1, x2
  # Objective function evaluated right AFTER x1, x2  (to see whether we have a local min, max or saddle point)
  gd_alpha_fix_mult <- function(X,              # X is c(x1, x2)
                                alpha = 0.1, 
                                tol = 0.00001){
    update_size <- 1
    res <- c(X, funct_mult(X-c(0.1,0.1)), funct_mult(X), funct_mult(X+c(0.1, 0.1)))
  
    while(update_size > tol){
      updater <- alpha * score_mult(X)
      update_size <- sqrt(sum((updater)^2))
      X <- X - updater 
      subres <- c(X, funct_mult(X - c(0.1, 0.1)), funct_mult(X), funct_mult(X + c(0.1, 0.1) ))
      res <- matrix(rbind(res, subres), 
                    ncol=5, 
                    dimnames = list(c(), c("x_1", "x_2", "f(X-0.1)", "f(X)", "f(X+0.1)")))
    }
    return(res)                              # returns a matrix with the solution parameters for each (!) iteration
  }
  
  gd_alpha_fix_mult(c(8, 2), tol=0.0001)     # It will converge to something, as long as the initial values are in the (-2, 2) range.
  
  
# Adapted Gradient-Descent with flexible learning parameter:
  # Output is a table containing values for every iteration. The values
  # included are: 
  # x1 and x2
  # Objective function evaluated right BEFORE x1, x2
  # Objective function evaluated at x1, x2
  # Objective function evaluated right AFTER x1, x2  (to see whether we have a local min, max or saddle point)
  gd_alpha_flex_mult <- function(X,                 # X is c(x1, x2)
                                 alpha=0.1, 
                                 tol=0.00001){
    update_size <- 1
    res <- c()
    while(update_size > tol){
      updater <- alpha * score_mult(X)
      update_size <- sqrt(sum((updater)^2))
      X_new <- X - updater
      
      if ((funct_mult(X) - funct_mult(X_new)) < 0){        # f(x) - f(x_new) < 0 means we are moving in the wrong direction
        alpha <- alpha * 0.5
      }
      else{
        alpha <- alpha * 1.05
        X <- X_new
      }
      
      subres <- c(X, funct_mult(X - c(0.1, 0.1)), funct_mult(X), funct_mult(X + c(0.1, 0.1) ))
      res <- matrix(rbind(res, subres), 
                    ncol=5, 
                    dimnames = list(c(), c("x_1", "x_2", "f(X-0.1)", "f(X)", "f(X+0.1)")))
    }
    return(res)
  }
  
  gd_alpha_flex_mult(c(8,4), tol = 0.00001)
  
  
# Adapted Nelder-Mead:
  nelder_mead_mult <- function(f,              # function to optimize
                               n = 2,          # number of variables in fx
                               tol = 0.0001,   # desired tolerance level
                               alpha = 1,      # standard alpha for Nelder-Mead
                               beta = 1/2,     # standard beta for Nelder-Mead
                               gamma = 2,      # standard gamma for Nelder-Mead
                               xmin = -3,      # min of interval where random number is drawn for initial value of parametervector
                               xmax = 3,       # max of interval where random number is drawn for initial value of parametervector
                               tau1 = 0.05,    # usual tau to compute initial simplex
                               tau0 = 0.00025  # tau to compute initial simplex if if value in the paramtervector on i-th place is zero
  ){
    # initial values for the first simplex:
    x_0 <- runif(seq(1:n), min = xmin, max = xmax)
    # construct matrix NM, where starting values for the initial simplex are stored, first column is x_0,
    # all the other values are set to zero and will later be replaced
    NM <- matrix(c(x_0, rep(0,n*n)), nrow = n)
    I <- diag(n)
    
    for (i in 1:n){
      if(NM[i,1] != 0){
        NM[,1+i] = NM[,1] + tau1*I[,i]
      }  
      else{
        NM[,1+i] = NM[,1] + tau0*I[,i]
      }
    }
    check <- tol+1            # define a variable check which is used as termination criterion and higher than tolerance level in the initial step and adjusted with every iteration
    res <- c()                # construct a vector res, where results from optimization procedure are stored at every iteration
    
    
    # Start iteration procedure:
    while(check > tol^2){         # start Nelder-Mead procedure and iterate as long as tolerance level is bigger than termination criterion
      results <- apply(NM,2,f)    # define a variable results to identify best, worst and second worst evaluation of the function
      a <-  order(results)        # command order gives the information about the ascending order of the elements in results
      
      # find points x_b, x_w, x_sw
      b <-  NM[,a[1]]               # best: lowest value of function
      sw <-  NM[,a[n]]              # second worst: second highest value of function
      w <-  NM[,a[n+1]]             # worst: highest value of function
      
      # find centroid of the n best points
      z <- 1/n * (apply(NM,1,sum)-w) # compute centroid
      r <- z + alpha*(z-w)           # compute reflection point r
      
      # Reflection if if object function evaluated at r is between the evaluation of the function at b and sw
      if (f(b) <= f(r) & f(r) <= f(sw)){
        NM[,a[n+1]] <- r            # If condition fulfilled, replace w with r
      }
      
      # Expansion if object function evaluated r is smaller than b
      if (f(r) < f(b)) {
        e <-  z + beta*(z-w)      # compute expansion point e
        if (f(e)<f(r)) {
          NM[,a[n+1]]<-e        # if f(e)<f(r) replace w with e
        }
        else {
          NM[,a[n+1]]<-r        # replace w with r
        }
      }
      # Contraction if object function evaluated at r is bigger than evaluation at sw
      if (f(r) > f(sw)) {
        if (f(r) > f(w)) {
          k <-  z + gamma*(w-z)   # compute contraction point k
        }
        else {
          k <-  z + gamma*(z-w)   # use different method to compute contraction point k
        }
        if (f(k) < f(w)) {      # compare object function evaluated at the contraction point with object function evaluated at w
          NM[,a[n+1]]<-k        # replace w with r
        }
        else {
          NM <- (NM+b)/2        # contract simplex by b
        }
      }
      fdash <- 1/(n+1)*sum(apply(NM,2,f))
      check <- 1/(n+1)*sum(  ( (apply(NM,2,f)-fdash)^2)  )          # compute a termination criteria
      res <- matrix(c(res,NM[,order(apply(NM,2,f))[1]]), nrow=n)      # store in res the parameters, where object function is lowest
    } # end of iteration step
    return(res) # return matrix that shows development of the parameter vector
  }
  
  
  # End of code