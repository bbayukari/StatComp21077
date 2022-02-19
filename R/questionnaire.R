library(abess)
library(aricode)

n <- 400
m.vec <- c(20,50,100)
b.vec <- c(0.2,0.5)

for(m in m.vec){
    for(b in b.vec){
        nmi <- numeric(100)
        for(i in 1:100){
            a <- rep(0,m)
            a[sample(1:m,5)] <- 1
            e <- rnorm(n)
            q <- matrix(sample(0:4,n*m,replace = TRUE),nrow = n,ncol = m)
            z <- rnorm(n,ifelse(q[,1]>2.5,-1,1))
            y <- b * q %*% a + z + e
            fit <- abess(cbind(q,z),y,always.include = c(m+1), support.size = 6)
            coef <- coef(fit, support.size = fit[["best.size"]])@i
            indix <- coef[2:(length(coef)-1)]
            a.fit <- rep(0,m)
            a.fit[indix] <- 1
            nmi[i] <- NMI(a,a.fit)
            #nmi[i] <- fit[["best.size"]]
        }
        cat("m=",m,"b=",b,"NMI=",mean(nmi)," \n")
    }
}
