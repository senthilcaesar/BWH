data <- mtcars

# Function to apply
mpg_category <- function(mpg) {
  if(mpg > 30) {
    return ("High")
  } else if (mpg > 20) {
    return ("Medium")
  }
  return ("Low")
}

# Apply to each element
lapply(X = data$mpg, FUN = mpg_category)

# Use sapply to simplify the result to a vector or matrix instead of a list
sapply(X = data$mpg, FUN = mpg_category)