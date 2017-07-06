###########################################################################
# PROGRAM: project_ML.R
#
# LANGUAGE*: R
#
# DESCRIPTION: This script replicated the results of the paper:
#              "High Dimensional Methods and Inference on Structural
#              and Treatment Effects".
#
# JOURNAL: Jornal of Economic Perspectives
#
# YEAR: 2014
#
# DETAILS: This code makes use of the LASSO method to select control variables
#          in the context of an Intrumental Variable method. To illustrate
#          how this can be achieved, the method is applied to problem 
#          raised in the paper "The Colonial Origins of Comparative 
#          Development: An Empirical Investigation", which investigate
#          the effects of institutions on economic performance. 
#
# WRITTEN BY: Cleyton Farias
#
# DATE: Jun 08, 2017
#
# LAST MODIFIED: Jun 23, 2017 
#
###########################################################################

source('functions_ML.R')
library(dplyr)
library(data.table)

foo <- fread("~/data.txt")

GDP = foo[, (GDP)];
Exprop = foo[, (Exprop)];	
Mort = foo[, (log(Mort))];
Latitude = foo[ ,(Latitude)];
Neo = foo[,(Neo)];
Africa = foo[,(Africa)];
Asia = foo[,(Asia)];
Namer = foo[, (Namer)];
Samer = foo[, (Samer)]



x <- data.table(Latitude, 
                Latitude2 = Latitude^2,  
                Latitude3 = Latitude^3, 
                (Latitude-.08)*((Latitude-.08) > 0),
                (Latitude-.16)*((Latitude-.16) > 0),
                (Latitude-.24)*((Latitude-.24) > 0),
                ((Latitude-.08)*((Latitude-.08) > 0))^2,
                ((Latitude-.16)*((Latitude-.16) > 0))^2, 
                ((Latitude-.24)*((Latitude-.24) > 0))^2,
                ((Latitude-.08)*((Latitude-.08) > 0))^3,
                ((Latitude-.16)*((Latitude-.16) > 0))^3,
                ((Latitude-.24)*((Latitude-.24) > 0))^3)


x = data.table(Africa, Asia, Namer, Samer, x)



#========================= BASELINE MODEL ==================================
z_B <- data.table(Mort, Latitude, rep(1, length(Mort)))
colnames(z_B) <- c('Mort', 'Latitude', 'Intercept')
z_B <- as.matrix(z_B)


#--------------------------- FIRST STAGE ----------------------------------
# First Stage: Exprop ~ Mort + CONTROLS  
FS_B <- qr.solve(z_B, Exprop) # Solving a non-square system: 
FS_B <- as.matrix(FS_B)
round(FS_B, 4)

Exprop <- as.matrix(Exprop)

# First stage estimator standard errors:
SEF_B <- hetero_se(z_B, Exprop-z_B%*%FS_B, solve(t(z_B)%*%z_B))
SEF_B <- as.matrix(SEF_B)
round(SEF_B, 4)


#---------------------------- SECOND STAGE --------------------------------- 
x_B <- data.table(Exprop, Latitude,
                     matrix(1, nrow = length(Mort), ncol = 1))

colnames(x_B) <- c('Exprop', 'Latitude', 'intercept')
x_B <- as.matrix(x_B)


# 
SS_B <- solve((t(z_B)%*%x_B), (t(z_B)%*%GDP))
round(SS_B, 4)

# 
SES_B <- hetero_se(z_B, GDP - x_B%*%SS_B, solve(t(z_B)%*%x_B))
SES_B <- as.matrix(SES_B)
round(SES_B, 4)
#==========================================================================







#=========================== ALL CONTROLS =================================
z_all <- data.table(Mort, x, 
                    intercept = rep(1, nrow = length(Mort)))

z_all <- as.matrix(z_all)


#-------------------------- FIRST STAGE -----------------------------------
FS_all <- qr.solve(z_all, Exprop)
FS_all <- as.matrix(FS_all)
round(FS_all, 4)

# Standard Errors: 
SEF_all <- hetero_se(z_all, Exprop - z_all%*%FS_all, 
                     solve(t(z_all)%*%z_all) )
SEF_all <- as.matrix(SEF_all)
round(SEF_all, 4)


#------------------------- SECOND STAGE -----------------------------------
x_all = data.table(Exprop = Exprop, x, rep(1, length(Mort)))
x_all <- as.matrix(x_all)

SS_all = solve(t(z_all)%*%x_all, t(z_all)%*%GDP)
SS_all <- as.matrix(SS_all)
round(SS_all, 4)

# Standard Errors:
SES_all = hetero_se(z_all, GDP - x_all%*%SS_all, solve(t(z_all)%*%x_all))
SES_all <- as.matrix(SES_all)
round(SES_all, 4)
#==========================================================================







#========================= VARIABLE SELECTION ==============================
n <- nrow(x)
p <- ncol(x)

# Demeaning the variables to perform variable selection: 
My <- as.matrix(GDP) - mean(GDP)
My <- as.matrix(My)

Mz <- as.matrix(Mort) - mean(Mort)
Mz <- as.matrix(Mz)

Md <- Exprop - mean(Exprop)
Md <- as.matrix(Md)

mean_x <- sapply(x, function(x) mean(x))
mean_x <- rep(1, n)%*%t(mean_x)
Mx <- x - mean_x
Mx <- as.matrix(Mx)


#=============================== NAIVE LASSO ===============================
library(glmnet) 
set.seed(2017)

lasso.gdp <- betterlasso(X = Mx, Y = My)
var.gdp <- lasso.gdp$variables[, 1] 
var.gdp

lasso.mort <- betterlasso(X = Mx, Y = Mz)
var.mort <- lasso.mort$variables[,1]
var.mort 

lasso.exprop <- betterlasso(X = Mx, Y = Md)
var.exprop <- lasso.exprop$variables[,1]
var.exprop

x <- data.table(x)
xSel <- x[, .(Africa, Asia, Namer, Samer, V4, V5, V6, V12)]


#--------------------------- FIRST STAGE --------------------------------
z_sel = data.table(Mort, xSel, rep(1, length(Mort)))
z_sel = as.matrix(z_sel)
FS_sel = qr.solve(z_sel, Exprop)
round(FS_sel, 4)

# FIRST STAGE Standard Errors: 
SEF_sel <- hetero_se(z_sel, Exprop - z_sel%*%FS_sel, 
                     solve(t(z_sel)%*%z_sel))

SEF_sel <- as.matrix(SEF_sel)
round(SEF_sel, 4)


#--------------------------- SECOND STAGE -------------------------------
x_sel = data.table(Exprop = Exprop, xSel, rep(1, length(Mort)))
x_sel <- as.matrix(x_sel)

SS_sel = solve(t(z_sel)%*%x_sel, t(z_sel)%*%GDP)
SS_sel <- as.matrix(SS_sel)
round(SS_sel, 4)

SES_sel = hetero_se(z_sel, GDP - x_sel%*%SS_sel, solve(t(z_sel)%*%x_sel))
SES_sel <- as.matrix(SES_sel)
round(SES_sel, 4)
#========================================================================




#======================= USING Cross-validation ==========================
n <- nrow(x)
p <- ncol(x)

# Demeaning the variables:
My = GDP - mean(GDP)
Mz = Mort - mean(Mort)
Md = Exprop - mean(Exprop)

mean_x <- sapply(x, function(x) mean(x))
mean_x <- rep(1, n)%*%t(mean_x)
mean_x
Mx <- x - mean_x

Mx <- as.matrix(Mx)
My <- as.matrix(My)
Mz <- as.matrix(Mz)
Md <- as.matrix(Md)


# Number of folders: 
k = 10

# Setting a seed: 
set.seed(666333)

# Folder Index:
folds <- sample(x =  1:k,size =  n, replace = T)


# Grid of penalty terms:
lambdaGrid = seq(20, 70, 1)
lambdaGrid <- as.matrix(lambdaGrid)
nG = length(lambdaGrid)


# Variables: 
y = GDP
d = as.numeric(Exprop)
z = Mort


# Cross-Validation Error Matrixes:               
CVy = matrix(NA, nrow = nG, ncol = k)
CVd = matrix(NA, nrow = nG, ncol = k)
CVz = matrix(NA, nrow = nG, ncol = k)
    

x <- data.table(x)
z <- data.table(Mort = z)
d <- data.table(Exprop = d)
y <- data.table(GDP = y)

tempo.inicio = Sys.time()

for(i in 1:nG){
    for(j in 1:k){
        
        # training set:      
        train_set <-  folds != j
        
        x_train <- x[train_set, ]
        z_train <- z[train_set, ]
        d_train <- d[train_set, ]
        y_train <- y[train_set, ]
        
        # Test set:
        test_set <- folds == j
        
        x_test <- x[test_set, ]
        z_test <- z[test_set, ]
        d_test <- d[test_set, ]
        y_test <- y[test_set, ]
                
        
        # demeaning training set:
        mean_x_train <- sapply(x_train, function(x) mean(x))
        mean_x_train <- rep(1, length(mean_x_train))%*%t(mean_x_train)
        Mx_train <- x_train - mean_x_train 
        
        Mz_train <- z_train - mean(z_train$Mort)
        Md_train <- d_train - mean(d_train$Exprop)
        My_train <- y_train - mean(y_train$GDP)
                
        Mx_train <- as.matrix(Mx_train)
        Mz_train <- as.matrix(Mz_train)
        Md_train <- as.matrix(Md_train)
        My_train <- as.matrix(My_train)
    
        # Performing Post-Lasso:
        PL.Z_train <- feasibleLasso(y = Mz_train, x = Mx_train,
                                        lambda = lambdaGrid[i], MaxIter = 100)

        PL.D_train <- feasibleLasso(y = Md_train, x = Mx_train,
                                        lambda = lambdaGrid[i], MaxIter = 100)

        PL.Y_train <- feasibleLasso(y = My_train, x = Mx_train,
                                        lambda = lambdaGrid[i], MaxIter = 100)
        
        # Variable Selection: 
        index.Z <- abs(PL.Z_train) > 0
        index.D <- abs(PL.D_train) > 0
        index.Y <- abs(PL.Y_train) > 0
       
        # intercept vector for the training set
        vone_train <- rep(1, nrow(x_train))
        
        # intercept vector for the test set
        vone_test <- rep(1, nrow(x_test))
        
        
        
        #------------ TEST SET
        x_test <- as.matrix(x_test)
        
        # Test set for Z:
        x_test.Z <- x_test[, index.Z, drop = F]
        if(length(x_test.Z) != 0){
            x_test.Z <- as.matrix(data.table(x_test.Z, intercept = vone_test))
        } else{
            x_test.Z <- as.matrix(data.table(intercept = vone_test))
        }
                
        # Test set for D:
        x_test.D <- x_test[, index.D, drop = F]
        if(length(x_test.D) != 0){
            x_test.D <- as.matrix(data.table(x_test.D, intercept = vone_test))
        } else{
            x_test.D <- as.matrix(data.table(intercept = vone_test))
        }
        
        # Test set for Y:
        x_test.Y <- x_test[, index.Y, drop = F]
        if(length(x_test.Y) != 0){
            x_test.Y <- as.matrix(data.table(x_test.Y, intercept = vone_test))
        } else{
            x_test.Y <- as.matrix(data.table(intercept = vone_test))
        }
        
        
        #----------------TRAINING SET
        x_train <- as.matrix(x_train)
        
        # Training set for Z:
        x_train.Z <- x_train[, index.Z, drop = F]
        if(length(x_train.Z) != 0){
            x_train.Z <- as.matrix(data.table(x_train.Z, intercept = vone_train))
        } else{
            x_train.Z <- as.matrix(data.table(intercept = vone_train))
        }
                
        # Training set for D:
        x_train.D <- x_train[, index.D, drop = F]
        if(length(x_train.D) != 0){
               x_train.D <- as.matrix(data.table(x_train.D, intercept = vone_train))
           } else{
               x_train.D <- as.matrix(data.table(intercept = vone_train))
           }
        
        # Training set for Y:
        x_train.Y <- x_train[, index.Y, drop = F]
        if(length(x_train.Y) != 0){
            x_train.Y <- as.matrix(data.table(x_train.Y, intercept = vone_train))
        } else{
            x_train.Y <- as.matrix(data.table(intercept = vone_train))
        }
        
           
         
        # Using least square estimator to predict (Post-LASSO): 
        # Gven E(e) = 0, then Y = XB -> B = (X_train^-1)*Y_train
        # Y^ = XB^ -> Y^ = X_test*[(X_train^-1)*Y_train]
        z_hat <- x_test.Z%*%qr.solve(x_train.Z, z_train)
        d_hat <- x_test.D%*%qr.solve(x_train.D, d_train)
        y_hat <- x_test.Y%*%qr.solve(x_train.Y, y_train)
        
        # Cross-validation MSE
        y_test <- as.matrix(y_test)
        d_test <- as.matrix(d_test)
        z_test <- as.matrix(z_test)

        CVy[i,j] = mean((y_test - y_hat)^2)
        CVd[i,j] = mean((d_test - d_hat)^2)
        CVz[i,j] = mean((z_test - z_hat)^2)
        }
}

tempo.final <- Sys.time()

tempo.execucao = tempo.final - tempo.inicio
tempo.execucao
                

# Média dos CV-predicted- values for each lambda:
meanCVy <- apply(CVy, 1, mean)
meanCVd <- apply(CVd, 1, mean)
meanCVz <- apply(CVz, 1, mean)


# Standard Deviation of CV-predicted-value by each lambda: 
stdCVy <- apply(CVy, 1, sd)
stdCVd <- apply(CVd, 1, sd)
stdCVz <- apply(CVz, 1, sd)


# t-inverse (Quantil function of t-student)
c68 <- qt(1-.16,9) # one standard error


# take the minimum CV for each regression:
min.CVy <- meanCVy == min(meanCVy)
min.CVd <- meanCVd == min(meanCVd)
min.CVz <- meanCVz == min(meanCVz)


# the lambda with the smallest mse +  1 standard deviation
by = max(meanCVy[min.CVy] + c68*stdCVy[min.CVy])
bd = max(meanCVd[min.CVd] + c68*stdCVd[min.CVd])
bz = max(meanCVz[min.CVz] + c68*stdCVz[min.CVz])


# Optimal choice for lAMBDA:
lambdaCVy = sqrt(10/9)*max(lambdaGrid[meanCVy < by])
lambdaCVd = sqrt(10/9)*max(lambdaGrid[meanCVd < bd])
lambdaCVz = sqrt(10/9)*max(lambdaGrid[meanCVz < bz])


# VARIANT OF LASSO WITH OPTIMAL LAMBDA
PL.My_CV <- feasibleLasso(y = My, x = Mx, lambdaCVy, MaxIter = 100)
PL.Md_CV <- feasibleLasso(y = Md, x = Mx, lambdaCVd, MaxIter = 100)
PL.Mz_CV <- feasibleLasso(y = Mz, x = Mx, lambdaCVz, MaxIter = 100)


index.My_CV <- abs(PL.My_CV) > 0
index.Md_CV <- abs(PL.Md_CV) > 0
index.Mz_CV <- abs(PL.Mz_CV) > 0

index.XCV <- data.table(as.numeric(index.My_CV), as.numeric(index.Md_CV),
                        as.numeric(index.Mz_CV))

index.XCV <- apply(t(index.XCV), 2, max)
index.XCV

# Matrix with all variables: 
x <- as.matrix(x)

# Matrix with the selected variables:
xCV = x[ ,index.XCV, drop = F]
xCV



#--------------------------- FIRST STAGE ----------------------------------
z_CV = data.table(Mort, xCV, intercept = rep(1, length(Mort)))
z_CV = as.matrix(z_CV)
FS_CV = qr.solve(z_CV, Exprop)
round(FS_CV, 4)

# FIRST STAGE Standard Errors: 
SEF_CV <- hetero_se(z_CV, Exprop - z_CV%*%FS_CV, 
                     solve(t(z_CV)%*%z_CV))

SEF_CV <- as.matrix(SEF_CV)
round(SEF_CV, 4)


#------------------------- SECOND STAGE -----------------------------------
x_CV = data.table(Exprop = Exprop, xCV, intercept = rep(1, length(Mort)))
x_CV <- as.matrix(x_CV)
x_CV

SS_CV = solve(t(z_CV)%*%x_CV, t(z_CV)%*%GDP)
SS_CV <- as.matrix(SS_CV)
round(SS_CV, 4)

SES_CV = hetero_se(z_CV, GDP - x_CV%*%SS_CV, solve(t(z_CV)%*%x_CV))
SES_CV <- as.matrix(SES_CV)
round(SES_CV, 4)
#========================================================================




#============== USING THE RULE AS IN BELLONI ET AL (2012) ================
n <- nrow(x)
p <- ncol(x)

# Parameters 
gamma <- .1/log(n)

# Penalty level: 
lambda <- 1.1*2*sqrt(2*n*(log(2*(p/gamma))))


# Demeaning the variables:
My = GDP - mean(GDP)
Mz = Mort - mean(Mort)
Md = Exprop - mean(Exprop)

mean_x <- sapply(x, function(x) mean(x))
mean_x <- rep(1, n)%*%t(mean_x)
Mx <- x - mean_x
Mx <- as.matrix(Mx)

My <- as.matrix(My)
Mz <- as.matrix(Mz)


PiMy = feasibleLasso(y = My, x = Mx, lambda = lambda, MaxIter = 100)

PiMd = feasibleLasso(y = Md, x = Mx, lambda = lambda, MaxIter = 100)

PiMz = feasibleLasso(y = Mz, x = Mx, lambda = lambda, MaxIter = 100)


IndMy = abs(PiMy) > 0
IndMy <- as.numeric(IndMy)

IndMd = abs(PiMd) > 0
IndMd <- as.numeric(IndMd)

IndMz = abs(PiMz) > 0
IndMz <- as.numeric(IndMz)

IndX <- data.table(IndMy,IndMd,IndMz)

IndX <- as.matrix(IndX)

# Getting the union of the selected variables in each equation:
index <- apply(t(IndX), 2, max)

# Matrix with all 16 variables: 
x <- as.matrix(x)

# Matrix with the selected variables:
xSel <- x[, index]


#--------------------------- FIRST STAGE --------------------------------
z_Sel = data.table(Mort, xSel, intercept = rep(1, length(Mort)))
z_Sel = as.matrix(z_Sel)

# Estimating the first stage:
FS_Sel = qr.solve(z_Sel, Exprop)
round(FS_Sel, 4)

# FIRST STAGE Standard Errors: 
SEF_Sel <- hetero_se(z_Sel, Exprop - z_Sel%*%FS_Sel, 
                     solve(t(z_Sel)%*%z_Sel))

SEF_Sel <- as.matrix(SEF_Sel)
round(SEF_Sel, 4)


#--------------------------- SECOND STAGE -------------------------------
x_Sel = data.table(Exprop = Exprop[, 1], xSel, intercept = rep(1, length(Exprop)))
x_Sel <- as.matrix(x_Sel)

# Estimating the coefficients:
SS_Sel = solve(t(z_Sel)%*%x_Sel, t(z_Sel)%*%GDP)
SS_Sel <- as.matrix(SS_Sel)
round(SS_Sel, 4)

SES_Sel = hetero_se(z_Sel, GDP - x_Sel%*%SS_Sel, solve(t(z_Sel)%*%x_Sel))
SES_Sel <- as.matrix(SES_Sel)
round(SES_Sel, 4)
#========================================================================
