
# This file contains the necessary code to run projetc_ML.R



# Function to calculate the heteroscedastic standard error: 
hetero_se <- function(X, e, XpXinv){
    p <- ncol(X)
    n <- nrow(X)
    
    V = t(X*((e^2)%*%matrix(1, nrow = 1, ncol = p)))%*%X
    
    vhetero = (n/(n-p))*(XpXinv%*%V%*%t(XpXinv))
    vhetero
    return(sqrt(diag(vhetero)))
}


# Function to perform cross-validation and return the variables selected by LASSO using the function glmnet(): 
betterlasso <- function(X, Y, folds = 10){
    cv.lasso <- cv.glmnet(x = X, y = Y, alpha = 1, nfold = folds)
    lambda.min <- cv.lasso$lambda.min
    lasso.model <- glmnet(x = X, y = Y, alpha = 1, lambda = lambda.min)
    c <- coef(lasso.model)
    variables <- rownames(c)
    c <- data.table(variable = variables, coefficients = as.matrix(c))
    c <- rename(c, coefficients = coefficients.s0)
    var <- c[ coefficients != 0, ]
    
    l = list(lambda = lambda.min, variables = var)
    return(l)
}


# Function to estimate LASSO coefficients (different than the glmnet function from glmnet package): 
LassoShooting2 <- function(X, y, lambda, Ups, maxIter = 10000, b0 = NULL){
    n <- nrow(X)
    
    p <- ncol(X)
    
    XX = t(X)%*%X
    
    Xy = t(X)%*%y
   
    if(is.null(b0) == TRUE){
        # beta inicial (start form the OLS solution):
        beta = qr.solve(XX + diag(as.matrix(lambda))*diag(p), Xy)
    } else{
        beta = b0
    }
    
    XX2 <- XX*2
    Xy2 <- Xy*2    
    i <- 0
    
    while(i < maxIter){
        
        beta_old <- beta
        
        for(j in 1:p){
            
            # Compute the Shoot and Update the variable
            S0 = sum(XX2[j, ]%*%beta) - XX2[j,j]*beta[j] - Xy2[j];
            
            if (S0 > lambda*Ups[j]){
                beta[j] = (lambda*Ups[j] - S0)/XX2[j,j]
            } else if(S0 < -lambda*Ups[j]){
                beta[j] = (-lambda*Ups[j] - S0)/XX2[j,j]
            } else if(abs(S0) <= lambda*Ups[j]){
                beta[j] = 0
            }
        }
        
        i <- i + 1
        
        # Check termination
        if(sum(abs(beta - beta_old)) < 1e-5){
            break
        }
        
    } # End of the while()
    
    return(round(beta, 5))   
}





feasibleLasso <- function(y, x, lambda, MaxIter = 100){
    n <- nrow(x)
    p <- ncol(x)
    
    # Score 0 inicial:
    Syx = x*(y%*%matrix(1, ncol = p, nrow = 1))
    
    # Ups0 inicial:
    Ups0 <- sqrt(apply(Syx^2, 2, mean))
    
    # beta inicial:
    b <- LassoShooting2(x, y, 0.5*lambda, Ups0)
    
    use <- abs(b) > 0
    
    # residual
    e <- y - x[, use]%*%qr.solve(x[, use], y)
    
    # inicio da iteracao
    i <- 1
    
    # score 1:
    Syx <- x*(e%*%matrix(1, ncol = p, nrow = 1))
    
    # ups 1:
    Ups1 <- sqrt(apply(Syx^2, 2, mean))
    
    #MaxIter = 100
    while( norm(as.matrix(Ups0 - Ups1), type = "2") > 1e-6 & i < MaxIter){
        d0 <- norm(as.matrix(Ups0 - Ups1), type = "2")
        
        # novo beta:
        b = LassoShooting2(x, y, lambda, Ups1, b0 = b)
        
        use <- abs(b) > 0
        
        # novo residuo
        e = y - x[, use]%*%qr.solve(x[, use], y)
        
        Ups0 <- Ups1
        
        Syx <- x*(e%*%matrix(1, ncol = p, nrow = 1))
        
        Ups1 <- sqrt(apply(Syx^2, 2, mean))
        
        i <- i + 1
        
        # Verifica se a nova rodada do beta melhorou:
        d1 = norm(as.matrix(Ups0 - Ups1), type = "2")
        
        if(d1 == d0){
            Ups0 = Ups1
        }
        
    } # end of the while()
    
    return(round(b, 5))
    
}


