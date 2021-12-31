library(StatComp21077)
library(ordinalNet)
library(abess)
n <- 1000
p <- 10
support.size <- 10
class.num <- 4

dataset <- abess::generate.data(n, p, support.size,
                         family = "ordinal",
                         class.num = class.num,
                         seed = 4)
y <- dataset$y
x <- dataset$x
coef.abess <- fit.ordinal(x,y,family = "ordinal",tune.type = "gic")
coef.abess <- fit.ordinal(x,y,family = "ordinal",tune.type = "cv")
coef.net <- ordinalNet(x,y,family = "cumulative",
                              link = "logit",
                              parallelTerms = TRUE,
                              nonparallelTerms = FALSE
                              )

cat("\n true:",dataset$intercept,dataset$beta)
cat("\n abess:",coef.abess)
cat("\n net:", coef.net$coefs)
#cat("\n ratio:",coef/c(dataset$intercept,dataset$beta))


cat("\n loss_true:",loss(x,y,dataset$beta,dataset$intercept))
cat("\n loss_abess:",loss(x,y,coef.abess[-1:-class.num+1],coef.abess[1:class.num-1]))
cat("\n loss_net:",-coef.net$loglik)

# dif_true <- dif(x,y,dataset$beta,dataset$intercept)
# dif_est <- dif(x,y,coef[-1:-class.num+1],coef[1:class.num-1])

# print("\n")
# print(dif_true)
# print(dif_est)
# #print(sum((dif_true-dif_est)*(dif_true-dif_est)))
#
# grad_true <- grad(x,y,dataset$beta,dataset$intercept)
# grad_est <- grad(x,y,coef[-1:-class.num+1],coef[1:class.num-1])

n <- seq(0.1,10,0.1)
lo <- numeric(length(n))
for(i in 1:length(n)){
    lo[i] = loss(x,y,n[i]*dataset$beta,n[i]*dataset$intercept)
}
plot(n,lo)
#
n2 <- seq(1,5,0.1)
lo2 <- numeric(length(n2))
for(i in 1:length(n2)){
    lo2[i] = loss(x,y,n2[i]*coef.abess[-1:-class.num+1],n2[i]*coef.abess[1:class.num-1])
}
plot(n2,lo2)
