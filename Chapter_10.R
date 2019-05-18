m <- 1; c <- 1; k <- 16.25; Delta <- 0.01
a0 <- m / Delta^2 + c / Delta + k
a1 <- -2 * m / Delta^2 - c / Delta; a2 <- m / Delta^2
n <- 100000
y <- c(0, 0); x <- c(0, 0)
set.seed(1)
for (i in 3:n) {
  x[i] <- x[i-1] - 0.5 * x[i-2] + rnorm(1)
  y[i] <- (-a1 * y[i-1] - a2 * y[i-2]) / a0 + x[i] / a0
}
nax <- nay <- 1
for (i in 1:n) {
  x[i] <- x[i] + nax * rnorm(1)
  y[i] <- y[i] + nay * rnorm(1)
}
Sxx <- spectrum(x, span = 31)
Syy <- spectrum(y, span = 31)
Gemp <- sqrt( Syy$spec[1:5000] / Sxx$spec[1:5000] )
Freq <- Syy$freq[1:5000]
FreH <- Freq / Delta
Omeg <- 2 * pi * Freq
OmegH <- 2 * pi * FreH
Gth <- sqrt( 1/( (k-m*OmegH^2)^2 + c^2*OmegH^2 ))
Gar <- 1 / abs( 1 + a1/a0 * exp(-Omeg*1i) + a2/a0 * exp(-Omeg*2i) )
plot(FreH, Gth, xlab = "Frequency (Hz)", ylab = "Gain", type="l")
lines(FreH, Gemp, lty = "dashed")
lines(FreH, Gar, lty = "dotted")