set.seed(2021)
n <- 1000; n0 <- 1000; m <- 1000
data = sim_data_generator_1D_example2(sim_name = "1D_setA", n = n, n0 = n0, m = m)
xtrain = data$xtrain
ytrain = data$ytrain
xcalibration = data$xcalibration
ycalibration = data$ycalibration
xtest = data$xtest
ytest = data$ytest
PItruth = data$truePI

D = matrix(0, n0, n0)
Dnew = matrix(0, m, n0)
DnewT = matrix(0, n0, m)
for(i in 1:n0){
  D[i,] = abs(xcalibration[i] - xcalibration)
}
for(i in 1:m){
  Dnew[i,] = abs(xtest[i] - xcalibration)
  DnewT[,i] = abs(xtest[i] - xcalibration)
}
eps = abs(ycalibration)
order1 = order(eps)
alpha = .05

Vcv = abs(ytrain)
Dcv = matrix(0, n, n)
for(i in 1:n){
  Dcv[i,] = abs(xtrain[i] - xtrain)
}
max0 = max(Dcv)*2; min0 =quantile(Dcv, 0.01)
hs = exp(seq(log(min0), log(max0), length.out = 20))

myLCR = LCPmodule$new(H =D[order1, order1], V = eps[order1], h = .2, alpha = alpha, type = "distance")
####auto-tuning with independent training set and regression score
auto_ret = myLCR$LCP_auto_tune(V0 =Vcv, H0 =Dcv, hs = hs, B = 2, delta =alpha/2, lambda = 1, trace = TRUE)
myLCR$h = auto_ret$h
myLCR$lower_idx()
myLCR$cumsum_unnormalized()
# Dnew is m by n (training) ordered distance matrix for m test samples and n training samples
myLCR$LCP_construction(Hnew =Dnew[,order1], HnewT = DnewT[order1,])
deltaLCP = myLCR$band_V
qLL = -deltaLCP
qLU = +deltaLCP

ll = order(xtest)
plot(xtest[ll], ytest[ll], type = "p", pch = ".", ylim = c(-3, 4), xlab = 'X', ylab = 'Y')
points(xtest[ll], PItruth[ll,1], type = 'l')
points(xtest[ll], PItruth[ll,2], type = 'l')

points(xtest[ll], qLL[ll], type = 'l', col = "red")
points(xtest[ll], qLU[ll], type = 'l', col = "red")

legend("topleft", legend = c("Truth", "CP", "LCP"), lty = 1, col = c("black","blue","red"), bty = "n")


