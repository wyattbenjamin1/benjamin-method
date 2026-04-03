# Benjamin Method - Demo
# Run this file to see the method in action

source("benjamin_method.R")

# --- Oil Filter Example ---
means <- c(14.5, 13.8, 13.3, 14.3, 13.1)
J <- 9
MSE <- 0.088

cat("Running oil filter comparison...\n")
result <- my_comparison(means, J, MSE, alpha = 0.05)
print(result)

# --- Run Full Analysis (FWER + Power) ---
run_all()
