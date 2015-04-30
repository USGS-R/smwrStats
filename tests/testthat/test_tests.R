context("Hypothesis Test functions")

test_that("Simple Hypothesis Tests working", {
  testthat::skip_on_cran()
  X9 <- c(4.56, 8.47, 9.01, 9.26, 9.55,13.02, 14.06, 14.34, 15.37)
  Y16 <- c(3, 2, 13, 6, 1, 5, 15, 7, 10, 4, 9, 11, 8, 12, 14, 16)
  T16 <- seq(1, 16)
  # Grubbs test
  expect_that(grubbs.test(X9)$statistic, 
              equals(c(G=1.77472737)))
  expect_that(grubbs.test(X9)$p.value, 
  						equals(0.453646459))
  # Kendall-Sen test
  expect_that(kensen.test(Y16, T16)$statistic, 
              equals(c(tau=0.483333333)))
  expect_that(as.vector(kensen.test(Y16, T16)$p.value), 
  						equals(0.00595023991))
#   expect_that(as.vector(kensen.test(Y15, T15)$estimate[5L]), 
#   						equals(19.8096079))
  # PPCC test
  expect_that(ppcc.test(X9)$statistic, 
  						equals(c(r=0.962735488)))
  expect_that(as.vector(ppcc.test(X9)$p.value), 
  						equals(0.403664504))
  # Seasonal Kendall test
  expect_that(seaken(Y16, 2)$statistic, 
              equals(c(tau=0.428571433)))
  expect_that(seaken(Y16, 2)$p.value, 
  						equals(0.04421174526214599609375))
  # Serial test
  expect_that(serial.test(Y16)$statistic,
  						equals(c(S=44.5)))
  expect_that(serial.test(Y16)$p.value,
  						equals(0.0193140193))
  expect_that(serial.test(Y16, method="runs")$statistic,
  						equals(c(S=5)))
  expect_that(serial.test(Y16, method="runs")$p.value,
  						equals(0.0768127441))
})

test_that("Correlation/Regression Tests working", {
	testthat::skip_on_cran()
	XM <- cbind(X1=seq(1, 4, length.out=30),
             X2=c(1.45, 2.72, 1.38, 2.18, 3.67, 2.71, 3.36, 2.20, 3.32, 2.86, 2.57,
             	    3.66, 3.40, 3.36, 2.74, 3.14, 1.54, 1.97, 2.40, 1.01, 1.20, 2.57,
                  1.82, 2.43, 2.55, 1.36, 3.87, 1.02, 3.89, 3.78))
	YT <- c(4.25, 2.05, 0.75, 3.69, 1.11, 3.14, 3.76, -3.80, -3.96,  1.85, -4.83,
          0.50, 3.18, 3.88, -5.28, 5.39, -0.33, 3.47, 3.04, 7.75, 2.52, 1.94, 2.70,
          -3.74, -1.57, 1.21, 0.12, 2.37, 3.20, 2.03)
	YM <- XM[,1L] + XM[,2L] + YT
	LM <- lm(YM ~ XM) # linear regression
	GM <- glm(YM > 6.4 ~ X1 + X2, data=as.data.frame(XM), family=binomial)
	# Output from cor.all
	expect_that(cor.all(XM)$estimates[1L, 2L],
							equals(-0.015497903))
	expect_that(cor.all(XM)$p.values[1L, 2L],
							equals(0.935216658))
	# Hosmer-Lemshow test
	expect_that(hosmerLemeshow.test(GM)$statistic,
							equals(c("Chi-square"=6.40319888)))
	expect_that(as.vector(hosmerLemeshow.test(GM)$p.value),
							equals(0.602163681))
})