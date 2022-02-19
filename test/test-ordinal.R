library(StatComp21077)
library(ordinalNet)
library(abess)
test.ordinal.gic <- function(n,p,seed,support.size = 10,class.num = 4){
    v <- round(n/2)
    repeat{
        dataset <- abess::generate.data(n+v, p, support.size,
                                        family = "ordinal",
                                        class.num = class.num,
                                        seed = seed)
        train <- sample(1:(n+v),n)
        y <- dataset$y[train]
        x <- dataset$x[train,]
        y.test <- dataset$y[-train]
        x.test <- dataset$x[-train,]
        if(length(unique(y)) == max(y)+1 && length(unique(y.test)) == max(y.test)+1){
            break
        }
        seed = seed + 1000
    }

    time.gic <- as.numeric(system.time(fit.gic <- abess(x,y,family = "ordinal",tune.type = "gic"))[3])

    beta.idx.true <- which(dataset[["beta"]]!=0)

    coef.gic <- coef(fit.gic, support.size = fit.gic[["best.size"]])[[1]]
    beta.idx.gic <- unique(coef.gic@i)[-1]
    intercept.gic <- as.vector(coef.gic[1,])
    beta.gic <- as.vector(coef.gic[-1,1])

    TPR.gic <- TPR(beta.idx.true,beta.idx.gic)
    TNR.gic <- TNR(beta.idx.true,beta.idx.gic,p)
    ReErr.gic <- ReErr(dataset[["beta"]],beta.gic)
    SLE.gic <- SLE(beta.idx.true,beta.idx.gic)
    BS.gic <- Brier.score(x.test,y.test,beta.gic,intercept.gic)


    MR.gic <- MR(x.test,y.test,beta.gic,intercept.gic)

    return (c(TPR.gic,TNR.gic,ReErr.gic,SLE.gic,BS.gic,MR.gic,time.gic))
}

test.ordinal.cv <- function(n,p,seed,support.size = 10,class.num = 4){
    v <- round(n/2)
    repeat{
        dataset <- abess::generate.data(n+v, p, support.size,
                                        family = "ordinal",
                                        class.num = class.num,
                                        seed = seed)
        train <- sample(1:(n+v),n)
        y <- dataset$y[train]
        x <- dataset$x[train,]
        y.test <- dataset$y[-train]
        x.test <- dataset$x[-train,]
        if(length(unique(y)) == max(y)+1 && length(unique(y.test)) == max(y.test)+1){
            break
        }
        seed = seed + 1000
    }

    time.cv <- as.numeric(system.time(fit.cv <- abess(x,y,family = "ordinal",tune.type = "cv"))[3])

    beta.idx.true <- which(dataset[["beta"]]!=0)

    coef.cv <- coef(fit.cv, support.size = fit.cv[["best.size"]])[[1]]
    beta.idx.cv <- unique(coef.cv@i)[-1]
    intercept.cv <- as.vector(coef.cv[1,])
    beta.cv <- as.vector(coef.cv[-1,1])

        TPR.cv <- TPR(beta.idx.true,beta.idx.cv)
    TNR.cv <- TNR(beta.idx.true,beta.idx.cv,p)
    ReErr.cv <- ReErr(dataset[["beta"]],beta.cv)
    SLE.cv <- SLE(beta.idx.true,beta.idx.cv)
    BS.cv <- Brier.score(x.test,y.test,beta.cv,intercept.cv)
    MR.cv <- MR(x.test,y.test,beta.cv,intercept.cv)

    return (c(TPR.cv,TNR.cv,ReErr.cv,SLE.cv,BS.cv,MR.cv,time.cv))
}

test.ordinal.net <- function(n,p,seed,support.size = 10,class.num = 4){
    v <- round(n/2)
    repeat{
        dataset <- abess::generate.data(n+v, p, support.size,
                                        family = "ordinal",
                                        class.num = class.num,
                                        seed = seed)
        train <- sample(1:(n+v),n)
        y <- dataset$y[train]
        x <- dataset$x[train,]
        y.test <- dataset$y[-train]
        x.test <- dataset$x[-train,]
        if(length(unique(y)) == max(y)+1 && length(unique(y.test)) == max(y.test)+1){
            break
        }
        seed = seed + 1000
    }

    time.net <- as.numeric(system.time(fit.net <- ordinalNet(x,as.factor(y),family = "cumulative",
                                                             link = "logit",
                                                             parallelTerms = TRUE,
                                                             nonparallelTerms = FALSE))[3])

    beta.idx.true <- which(dataset[["beta"]]!=0)

    coef.net <- coef(fit.net,matrix = TRUE)
    beta.net <- coef.net[-1,1]
    beta.idx.net <- which(beta.net!=0)
    intercept.net <- coef.net[1,]

    TPR.net <- TPR(beta.idx.true,beta.idx.net)
    TNR.net <- TNR(beta.idx.true,beta.idx.net,p)
    ReErr.net <- ReErr(dataset[["beta"]],beta.net)
    SLE.net <- SLE(beta.idx.true,beta.idx.net)
    BS.net <- Brier.score(x.test,y.test,beta.net,intercept.net)
    MR.net <- MR(x.test,y.test,beta.net,intercept.net)

    return (c(TPR.net,TNR.net,ReErr.net,SLE.net,BS.net,MR.net,time.net))
}




## start test
set.seed(123)
R <- 3
n.total <- c(rep(500,8),100,200,400,800,1600)
p.total <- c(seq(20,40,5),500,1500,2500,rep(500,5))
num <- length(n.total)
seed.total <- sample(1:(num*R*10),num*R)

gic.data <- matrix(0,nrow = num*R,ncol = 9)
for(i in 1:R){
    for(j in 1:num){
        id <- (i-1)*num+j
        gic.data[id,] <- c(test.ordinal.gic(n.total[j],p.total[j],seed.total[id]),n.total[j],p.total[j],j,1)
    }
}

cv.data <- matrix(0,nrow = num*R,ncol = 9)
for(i in 1:R){
    for(j in 1:num){
        id <- (i-1)*num+j
        cv.data[id,] <- c(test.ordinal.cv(n.total[j],p.total[j],seed.total[id]),n.total[j],p.total[j],j,2)
    }
}

net.data <- matrix(0,nrow = num*R,ncol = 9)
for(i in 1:R){
    for(j in 1:num){
        id <- (i-1)*num+j
        net.data[id,] <- c(test.ordinal.net(n.total[j],p.total[j],seed.total[id]),n.total[j],p.total[j],j,3)
    }
}

low <- data.frame(rbind(gic.data[gic.data[,10]%in%1:5,c(-10,-8)],gic.data[gic.data[,10]%in%1:5,c(-10,-8)],gic.data[gic.data[,10]%in%1:5,c(-10,-8)]))
high <- data.frame(rbind(gic.data[gic.data[,10]%in%6:8,c(-10,-8)],gic.data[gic.data[,10]%in%6:8,c(-10,-8)],gic.data[gic.data[,10]%in%6:8,c(-10,-8)]))
limit <- data.frame(rbind(gic.data[gic.data[,10]%in%9:13,c(-10,-9)],gic.data[gic.data[,10]%in%9:13,c(-10,-9)],gic.data[gic.data[,10]%in%9:13,c(-10,-9)]))

colnames(low) <- c("TPR","TNR","ReErr","SLE","BS","MR","time","p","method")
colnames(high) <- c("TPR","TNR","ReErr","SLE","BS","MR","time","p","method")
colnames(limit) <- c("TPR","TNR","ReErr","SLE","BS","MR","time","n","method")

# function(x){
#     x[x==1] <- "gic"
#     x[x==2] <- "cv"
#     x[x==3] <- "net"
#     x <- factor(x)
# }
#
#
# ordinal.data <- list(low = dat)






