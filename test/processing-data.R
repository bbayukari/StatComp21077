r <- 20
low.dim.data <- data.frame(
    TPR = c(TPR[[1]][1,1:R],TPR[[2]][1,1:R],TPR[[3]][1,1:R],TPR[[4]][1,1:R],TPR[[5]][1,1:R],TPR[[1]][2,1:r],TPR[[2]][2,1:r],TPR[[3]][2,1:r],TPR[[4]][2,1:r],TPR[[5]][2,1:r]),
    TNR = c(TNR[[1]][1,1:R],TNR[[2]][1,1:R],TNR[[3]][1,1:R],TNR[[4]][1,1:R],TNR[[5]][1,1:R],TNR[[1]][2,1:r],TNR[[2]][2,1:r],TNR[[3]][2,1:r],TNR[[4]][2,1:r],TNR[[5]][2,1:r]),
    ReErr = c(ReErr[[1]][1,1:R],ReErr[[2]][1,1:R],ReErr[[3]][1,1:R],ReErr[[4]][1,1:R],ReErr[[5]][1,1:R],ReErr[[1]][2,1:r],ReErr[[2]][2,1:r],ReErr[[3]][2,1:r],ReErr[[4]][2,1:r],ReErr[[5]][2,1:r]),
    SLE = c(SLE[[1]][1,1:R],SLE[[2]][1,1:R],SLE[[3]][1,1:R],SLE[[4]][1,1:R],SLE[[5]][1,1:R],SLE[[1]][2,1:r],SLE[[2]][2,1:r],SLE[[3]][2,1:r],SLE[[4]][2,1:r],SLE[[5]][2,1:r]),
    BS = c(BS[[1]][1,1:R],BS[[2]][1,1:R],BS[[3]][1,1:R],BS[[4]][1,1:R],BS[[5]][1,1:R],BS[[1]][2,1:r],BS[[2]][2,1:r],BS[[3]][2,1:r],BS[[4]][2,1:r],BS[[5]][2,1:r]),
    MR = c(MR[[1]][1,1:R],MR[[2]][1,1:R],MR[[3]][1,1:R],MR[[4]][1,1:R],MR[[5]][1,1:R],MR[[1]][2,1:r],MR[[2]][2,1:r],MR[[3]][2,1:r],MR[[4]][2,1:r],MR[[5]][2,1:r]),
    time = c(time[[1]][1,1:R],time[[2]][1,1:R],time[[3]][1,1:R],time[[4]][1,1:R],time[[5]][1,1:R],time[[1]][2,1:r],time[[2]][2,1:r],time[[3]][2,1:r],time[[4]][2,1:r],time[[5]][2,1:r]),
    p = c(rep(20,R),rep(25,R),rep(30,R),rep(35,R),rep(40,R),rep(20,r),rep(25,r),rep(30,r),rep(35,r),rep(40,r)),
    method = c(rep("gic",R*5),rep("net",r*5))
)
low.dim.data$method = factor(low.dim.data$method)
low.dim.data$p = factor(low.dim.data$p)




high.dim.data <- data.frame(
    TPR = c(TPR[[6]][1,1:R],TPR[[7]][1,1:R],TPR[[8]][1,1:R],TPR[[6]][2,1:r],TPR[[7]][2,1:r],TPR[[8]][2,1:r]),
    TNR = c(TNR[[6]][1,1:R],TNR[[7]][1,1:R],TNR[[8]][1,1:R],TNR[[6]][2,1:r],TNR[[7]][2,1:r],TNR[[8]][2,1:r]),
    ReErr = c(ReErr[[6]][1,1:R],ReErr[[7]][1,1:R],ReErr[[8]][1,1:R],ReErr[[6]][2,1:r],ReErr[[7]][2,1:r],ReErr[[8]][2,1:r]),
    SLE = c(SLE[[6]][1,1:R],SLE[[7]][1,1:R],SLE[[8]][1,1:R],SLE[[6]][2,1:r],SLE[[7]][2,1:r],SLE[[8]][2,1:r]),
    BS = c(BS[[6]][1,1:R],BS[[7]][1,1:R],BS[[8]][1,1:R],BS[[6]][2,1:r],BS[[7]][2,1:r],BS[[8]][2,1:r]),
    MR = c(MR[[6]][1,1:R],MR[[7]][1,1:R],MR[[8]][1,1:R],MR[[6]][2,1:r],MR[[7]][2,1:r],MR[[8]][2,1:r]),
    time = c(time[[6]][1,1:R],time[[7]][1,1:R],time[[8]][1,1:R],time[[6]][2,1:r],time[[7]][2,1:r],time[[8]][2,1:r]),
    p = c(rep(500,R),rep(1500,R),rep(2500,R),rep(500,r),rep(1500,r),rep(2500,r)),
    method = c(rep("gic",R*3),rep("net",r*3))
)
high.dim.data$method = factor(high.dim.data$method)
high.dim.data$p = factor(high.dim.data$p)


limit.data <- data.frame(
    TPR = c(TPR[[9]][1,1:R],TPR[[10]][1,1:R],TPR[[11]][1,1:R],TPR[[12]][1,1:R],TPR[[13]][1,1:R],TPR[[9]][2,1:r],TPR[[10]][2,1:r],TPR[[11]][2,1:r],TPR[[12]][2,1:r],TPR[[13]][2,1:r]),
    TNR = c(TNR[[9]][1,1:R],TNR[[10]][1,1:R],TNR[[11]][1,1:R],TNR[[12]][1,1:R],TNR[[13]][1,1:R],TNR[[9]][2,1:r],TNR[[10]][2,1:r],TNR[[11]][2,1:r],TNR[[12]][2,1:r],TNR[[13]][2,1:r]),
    ReErr = c(ReErr[[9]][1,1:R],ReErr[[10]][1,1:R],ReErr[[11]][1,1:R],ReErr[[12]][1,1:R],ReErr[[13]][1,1:R],ReErr[[9]][2,1:r],ReErr[[10]][2,1:r],ReErr[[11]][2,1:r],ReErr[[12]][2,1:r],ReErr[[13]][2,1:r]),
    SLE = c(SLE[[9]][1,1:R],SLE[[10]][1,1:R],SLE[[11]][1,1:R],SLE[[12]][1,1:R],SLE[[13]][1,1:R],SLE[[9]][2,1:r],SLE[[10]][2,1:r],SLE[[11]][2,1:r],SLE[[12]][2,1:r],SLE[[13]][2,1:r]),
    BS = c(BS[[9]][1,1:R],BS[[10]][1,1:R],BS[[11]][1,1:R],BS[[12]][1,1:R],BS[[13]][1,1:R],BS[[9]][2,1:r],BS[[10]][2,1:r],BS[[11]][2,1:r],BS[[12]][2,1:r],BS[[13]][2,1:r]),
    MR = c(MR[[9]][1,1:R],MR[[10]][1,1:R],MR[[11]][1,1:R],MR[[12]][1,1:R],MR[[13]][1,1:R],MR[[9]][2,1:r],MR[[10]][2,1:r],MR[[11]][2,1:r],MR[[12]][2,1:r],MR[[13]][2,1:r]),
    time = c(time[[9]][1,1:R],time[[10]][1,1:R],time[[11]][1,1:R],time[[12]][1,1:R],time[[13]][1,1:R],time[[9]][2,1:r],time[[10]][2,1:r],time[[11]][2,1:r],time[[12]][2,1:r],time[[13]][2,1:r]),
    n = c(rep(100,R),rep(200,R),rep(400,R),rep(800,R),rep(1600,R),rep(100,r),rep(200,r),rep(400,r),rep(800,r),rep(1600,r)),
    method = c(rep("gic",R*5),rep("net",r*5))
)
limit.data$method = factor(limit.data$method)
limit.data$n = factor(limit.data$n)


ordinal.data <- list(low = low.dim.data, high = high.dim.data, limit = limit.data)
save(ordinal.data,file = "ordinal-data.Rdata")

