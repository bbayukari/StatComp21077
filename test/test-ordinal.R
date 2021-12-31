library(StatComp21077)
library(ordinalNet)
library(abess)
test.ordinal <- function(n,p,seed,support.size = 10,class.num = 4){
    v <- round(n/2)
    repeat{
        dataset <- abess::generate.data(n+v, p, support.size,
                                    family = "ordinal",
                                    class.num = class.num,
                                    seed = seed)
        if(unique(dataset$y) == max(dataset$y)+1){
            break
        }
    }

    train <- sample(1:(n+v),n)
    y <- dataset$y[train]
    x <- dataset$x[train,]
    y.test <- dataset$y[-train]
    x.test <- dataset$x[-train,]

    time.gic <- as.numeric(system.time(fit.gic <- abess(x,y,family = "ordinal",tune.type = "gic"))[3])
    #time.cv <- as.numeric(system.time(fit.cv <- abess(x,y,family = "ordinal",tune.type = "cv"))[3])
    time.net <- as.numeric(system.time(fit.net <- ordinalNet(x,as.factor(y),family = "cumulative",
                           link = "logit",
                           parallelTerms = TRUE,
                           nonparallelTerms = FALSE))[3])

    beta.idx.true <- which(dataset[["beta"]]!=0)

    coef.gic <- coef(fit.gic, support.size = fit.gic[["best.size"]])[[1]]
    beta.idx.gic <- unique(coef.gic@i)[-1]
    intercept.gic <- as.vector(coef.gic[1,])
    beta.gic <- as.vector(coef.gic[-1,1])

    # coef.cv <- coef(fit.cv, support.size = fit.cv[["best.size"]])[[1]]
    # beta.idx.cv <- unique(coef.cv@i)[-1]
    # intercept.cv <- as.vector(coef.cv[1,])
    # beta.cv <- as.vector(coef.cv[-1,1])

    coef.net <- coef(fit.net,matrix = TRUE)
    beta.net <- coef.net[-1,1]
    beta.idx.net <- which(beta.net!=0)
    intercept.net <- coef.net[1,]


    TPR.gic <- TPR(beta.idx.true,beta.idx.gic)
    #TPR.cv <- TPR(beta.idx.true,beta.idx.cv)
    TPR.net <- TPR(beta.idx.true,beta.idx.net)

    TNR.gic <- TNR(beta.idx.true,beta.idx.gic,p)
    #TNR.cv <- TNR(beta.idx.true,beta.idx.cv,p)
    TNR.net <- TNR(beta.idx.true,beta.idx.net,p)

    ReErr.gic <- ReErr(dataset[["beta"]],beta.gic)
    #ReErr.cv <- ReErr(dataset[["beta"]],beta.cv)
    ReErr.net <- ReErr(dataset[["beta"]],beta.net)

    SLE.gic <- SLE(beta.idx.true,beta.idx.gic)
    #SLE.cv <- SLE(beta.idx.true,beta.idx.cv)
    SLE.net <- SLE(beta.idx.true,beta.idx.net)

    BS.gic <- Brier.score(x.test,y.test,beta.gic,intercept.gic)
    #BS.cv <- Brier.score(x.test,y.test,beta.cv,intercept.cv)
    BS.net <- Brier.score(x.test,y.test,beta.net,intercept.net)

    MR.gic <- MR(x.test,y.test,beta.gic,intercept.gic)
    #MR.cv <- MR(x.test,y.test,beta.cv,intercept.cv)
    MR.net <- MR(x.test,y.test,beta.net,intercept.net)

    #MAE.gic <- MAE(x.test,y.test,beta.gic,intercept.gic)
    #MAE.cv <- MAE(x.test,y.test,beta.cv,intercept.cv)
    #MAE.net <- MAE(x.test,y.test,beta.net,intercept.net)


    # result <-
    # matrix(data = c(TPR.gic,TNR.gic,ReErr.gic,SLE.gic,BS.gic,MR.gic,MAE.gic,time.gic,
    #                 TPR.cv,TNR.cv,ReErr.cv,SLE.cv,BS.cv,MR.cv,MAE.cv,time.cv,
    #                 TPR.net,TNR.net,ReErr.net,SLE.net,BS.net,MR.net,MAE.net,time.net),
    #        nrow = 3,byrow = TRUE)
    result <-
    matrix(data = c(TPR.gic,TNR.gic,ReErr.gic,SLE.gic,BS.gic,MR.gic,time.gic,
                    TPR.net,TNR.net,ReErr.net,SLE.net,BS.net,MR.net,time.net),
           nrow = 2,byrow = TRUE)

    return (result)
}




## start test
set.seed(1)
R <- 1
n.total <- c(rep(500,8),100,200,400,800,1600)
p.total <- c(seq(20,40,5),500,1500,2500,rep(500,5))
num <- length(n.total)
seed <- sample(1:(num*R),num*R)
seed.idx <- 1

TPR     <- list(length=num)
TNR <- list(length=num)
ReErr <- list(length=num)
SLE <- list(length=num)
BS <- list(length=num)
MR <- list(length=num)
#MAE <- list(length=num)
time <- list(length=num)



for(test.idx in 1:num){
    n <- n.total[test.idx]
    p <- p.total[test.idx]
    TPR[[test.idx]] <- matrix(0,nrow = 2,ncol = R)
    TNR[[test.idx]] <- matrix(0,nrow = 2,ncol = R)
    ReErr[[test.idx]] <- matrix(0,nrow = 2,ncol = R)
    SLE[[test.idx]] <- matrix(0,nrow = 2,ncol = R)
    BS[[test.idx]] <- matrix(0,nrow = 2,ncol = R)
    MR[[test.idx]] <- matrix(0,nrow = 2,ncol = R)
    #MAE[[test.idx]] <- matrix(0,nrow = 2,ncol = R)
    time[[test.idx]] <- matrix(0,nrow = 2,ncol = R)

    for(i in 1:R){
        result <- test.ordinal(n,p,seed[seed.idx])
        seed.idx = seed.idx+1

        TPR[[test.idx]][,i] <- result[,1]
        TNR[[test.idx]][,i] <- result[,2]
        ReErr[[test.idx]][,i] <- result[,3]
        SLE[[test.idx]][,i] <- result[,4]
        BS[[test.idx]][,i] <- result[,5]
        MR[[test.idx]][,i] <- result[,6]
        #MAE[[test.idx]][,i] <- result[,7]
        time[[test.idx]][,i] <- result[,7]
        cat("\nidx =",test.idx,", i =",i)
    }
}








