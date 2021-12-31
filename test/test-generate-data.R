library(StatComp21077)

n <- 1000
dataset <- generate.data(n, 2, 2,
                         family = "ordinal",
                         class.num = 5,
                         seed = 132)


print(table(dataset$y))
plot(dataset$x,col = dataset$y+1)
