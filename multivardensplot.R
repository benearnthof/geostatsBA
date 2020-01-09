library(MASS)
s <- matrix(data = c(1, 0.8, 0.8, 1), nrow = 2)
x <- mvrnorm(100000, mu = c(0, 0), Sigma = s)

plot(x)

den3d <- kde2d(x[,1], x[,2])

library(magrittr)
library(plotly)
plot_ly(x=den3d$x, y=den3d$y, z=den3d$z) %>% add_surface()

library(mvtnorm)

dmvnorm(c(0,0), mean = c(0, 0), sigma = s)
n <- 50
x = seq(-3, 3, length= n)
y = x

z <- matrix(nrow = n, ncol = n)
for (i in 1:n) {
  for (j in 1:n) {
    z[i, j] <- dmvnorm(c(x[i], y[j]), mean = c(0, 0), sigma = s)
  }
}

f = function(x, y) {
  dmvnorm(c(x, y), mean = c(0, 0), sigma = s)
}

df <- expand.grid(x, y)
library(purrr)
df$Var3 <- map2_dbl(df$Var1, df$Var2, f)

library("fields")
library(scatterplot3d)

#Method 1
#Perspective Plot
# z <- df$Var3
# x <- df$Var1
# y <- df$Var2
persp(x,y,z,col="lightblue",main="Perspective Plot")

#Method 2
#Contour Plot
contour(x,y,z,main="Contour Plot")
filled.contour(x,y,z,color=terrain.colors,main="Contour Plot",)

#Method 3
#Heatmap
image(x,y,z,main="Heat Map")
image.plot(x,y,z,main="Heat Map")

#Method 4
#3-D Scatter Plot
X = expand.grid(x,y)
x = X[,1]
y = X[,2]
z = c(z)
scatterplot3d(x,y,z,color="lightblue",pch=21,main="3-D Scatter Plot")
