# vector field visualization
vector_field <- function(
  f,  # Function describing the vector field
  xmin=0, xmax=1, ymin=0, ymax=1,
  width=600, height=600,
  iterations=50,
  epsilon=.01,
  trace=TRUE
) {
  z <- matrix(runif(width*height),nr=height)
  i_to_x <- function(i) xmin + i / width  * (xmax - xmin)
  j_to_y <- function(j) ymin + j / height * (ymax - ymin)
  x_to_i <- function(x) pmin( width,  pmax( 1, floor( (x-xmin)/(xmax-xmin) * width  ) ) )
  y_to_j <- function(y) pmin( height, pmax( 1, floor( (y-ymin)/(ymax-ymin) * height ) ) )
  i <- col(z)
  j <- row(z)
  x <- i_to_x(i)
  y <- j_to_y(j)
  res <- z
  for(k in 1:iterations) {
    v <- matrix( f(x, y), nc=2 )
    x <- x+.01*v[,1]
    y <- y+.01*v[,2]
    i <- x_to_i(x)
    j <- y_to_j(y)
    res <- res + z[cbind(i,j)]
    if(trace) {
      cat(k, "/", iterations, "\n", sep="")
      dev.hold()
      image(res)
      dev.flush()
    }
  }
  if(trace) {
    dev.hold()
    image(res>quantile(res,.6), col=0:1)
    dev.flush()
  }
  res
}

# Sample data
van_der_Pol <- function(x,y, mu=1) c(cos(y-pi/4), sin(x-0.5))
res <- vector_field(
  van_der_Pol,
  xmin=-3, xmax=3, ymin=-3, ymax=3,
  width=800, height=800,
  iterations=500,
  epsilon=.01
)
image(-res)
