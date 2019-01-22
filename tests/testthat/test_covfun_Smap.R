library(sparseGP)
context("test covariance function with latent space")

gp <- list(
  x = data.frame(chan = sample(c("xs1","xs2","xs3","xs4"), 10, replace=TRUE),
                 E = runif(10,0,5)),
  y = data.frame(chan = sample(c("xs1","xs2","xs3","xs4"), 10, replace=TRUE),
                 E = runif(10,0,5)),
  hyperpars = list(

    gpFuns = list(
      base1 = list(covfun=getMatrix, deriveCovfunHyp=deriveMatrixHyp),
      base2 = list(covfun=getMatrix, deriveCovfunHyp=deriveMatrixHyp)
    ),
    hyperpars = list(
      base1 = list(
        sigma = 0.5,
        len = 2,
        nugget = 0.01
      ),
      base2 = list(
        sigma = 1,
        len = 5,
        nugget = 0.01
      )
    ),
    baseChannels = c("base1", "base2"),
    allChannels = c("xs1","xs2","xs3","xs4"),
    Smat = matrix(c(
      1.0, 0.0,
      0.0, 1.0,
      1.0, 0.0,
      1.0, 1.0),
      nrow=4, ncol=2, byrow=TRUE)
  )
)


test_that("element-wise and block-wise calculation of covariance agree",{

  covmat1 <- with(gp, getMatrixSmap(x, y, hyperpars, FALSE))
  mygrid <- expand.grid(i=seq(nrow(gp$x)), j=seq(nrow(gp$y)))
  covmat2 <- with(gp, getMatrixSmap(x[mygrid$i,], y[mygrid$j,], hyperpars, TRUE))
  expect_equivalent(covmat1,covmat2)
})

# TODO: more tests...

