#This code estimates persister frequency, death rate for non-persisters in a similar way as in Windels et al. (2019).

library(minpack.lm)

# Define the model function
x=c(0,1,2,5,8,24)
# Plug in raw survival measurements from experiments
y=c(1.000E+00,
    4.459E-04,
    2.420E-05,
    2.293E-05,
    8.790E-06,
    2.548E-07)
db=data.frame(x,y)

# Fit the model to the data
fit_result = nlsLM(log10(y) ~ log10((1-P0) * exp(-k*x) + P0 * exp(-p*x)),
                     start=list(P0=1e-4, k=5.097, p=0.3),
                     data = db)
summary(fit_result)

# Plot the data and the fitted curve
plot(x, log10(y), pch = 16, main = "Nonlinear Least Squares Fit")
curve(predict(fit_result, newdata = data.frame(x = x)), add = TRUE, col = "red")
legend("topright", legend = "Fitted Curve", col = "red", lty = 1)
