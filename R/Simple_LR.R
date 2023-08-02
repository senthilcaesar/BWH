library(dplyr)
model_equation <- function(model, ...) {
  format_args <- list(...)
  
  model_coeff <- model$coefficients
  format_args$x <- abs(model$coefficients)
  model_coeff_sign <- sign(model_coeff)
  model_coeff_prefix <- case_when(model_coeff_sign == -1 ~ " - ",
                                  model_coeff_sign == 1 ~ " + ",
                                  model_coeff_sign == 0 ~ " + ")
  model_eqn <- paste(strsplit(as.character(model$call$formula), "~")[[2]], # 'y'
                     "=",
                     paste(if_else(model_coeff[1]<0, "- ", ""),
                           do.call(format, format_args)[1],
                           paste(model_coeff_prefix[-1],
                                 do.call(format, format_args)[-1],
                                 " * ",
                                 names(model_coeff[-1]),
                                 sep = "", collapse = ""),
                           sep = ""))
  return(model_eqn)
}

# Research Question
# Does the amount of sleep hours affect how much money people spend on food shopping ?

df <-  data.frame(sleepHours = c(5, 5.5, 6, 6, 7, 7, 7.5, 8, 8.5, 9),
                  dollarsSpent = c(47, 53, 52, 44, 39, 49, 50, 38, 43, 40))

plot(df$sleepHours, df$dollarsSpent, pch = 16, cex = 1.3, col = "blue", xlim=c(4,10), ylim=c(36,54))
relation <- lm(dollarsSpent ~ sleepHours, data=df)
abline(lm(dollarsSpent ~ sleepHours, data=df))

print(coef(summary(relation)))

model_equation(relation, digits = 4, trim = TRUE)

coefficient <- coef(summary(relation))
pValue <- coefficient[2,c(4)]
print(paste0("P-values: ", pValue))

SSE <-  summary(relation)$r.squared
print(paste0("R-squared: ", SSE))

cat("The model is not a statistically significant fit to the data.
    This may be due to the small sample size.")