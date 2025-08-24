# Load necessary library
library(ggplot2)

## -------Part a: Simulate Response Data ----------##

# Define sample sizes for two cancer types based on n = 1e3
n1 <- 480 # Sample size for lung cancer 
n2 <- 520  # Sample size for colorectal cancer

# Define true response rates 
pi1_true <- 0.7  
pi2_true <- 0.3  

# Simulate number of responders using Binomial distribution
x1 <- rbinom(1, n1, pi1_true)  
x2 <- rbinom(1, n2, pi2_true)  

##--------- Part b: Assume Prior------------------##

# Define weakly informative Beta prior (Beta(1,1)  == U[0,1])
a_prior <- 1
b_prior <- 1

##------ Part c: Posterior Probability------------##
#(Assuming independent beta-binomial frameworks)

# Compute posterior parameters (based on assumptions)
a_post_1 <- a_prior + x1
b_post_1 <- b_prior + (n1 - x1)

a_post_2 <- a_prior + x2
b_post_2 <- b_prior + (n2 - x2)

# Generate posterior samples using Beta distribution
samples_1 <- rbeta(10000, a_post_1, b_post_1)
samples_2 <- rbeta(10000, a_post_2, b_post_2)

# Compute probabilities P(π > 0.4)
p_greater_0.4_1 <- mean(samples_1 > 0.4)
p_greater_0.4_2 <- mean(samples_2 > 0.4)

# Print results
cat("Probability that π_1 > 0.4:", p_greater_0.4_1, "\n")
cat("Probability that π_2 > 0.4:", p_greater_0.4_2, "\n")

## ------- Decision Rule ---------- ##

decision <- ifelse(p_greater_0.4_1 > 0.8 | p_greater_0.4_2 > 0.8, "Continue Trial", "Consider Stopping")
cat("Decision:", decision, "\n")

## Step 5: Visualization ---------- ##

# Create a data frame for visualization
posterior_df <- data.frame(
  Samples = c(samples_1, samples_2),
  Cancer_Type = rep(c("Lung Cancer", "Colorectal Cancer"), each = 10000)
)

# Plot posterior distributions
ggplot(posterior_df, aes(x = Samples, fill = Cancer_Type)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0.4, linetype = "dashed", color = "red") +
  labs(title = "Posterior Distributions of Response Rates",
       x = "Response Rate (π)",
       y = "Density") +
  theme_minimal()
 

## ------ Creating Functions ------ ##
run_simu = function(pi1, pi2, n1 = 480, n2 = 520)
{
  # Simulate number of responders using Binomial distribution
  x1 <- rbinom(1, n1, pi1)  
  x2 <- rbinom(1, n2, pi2)  
  
  # Compute posterior parameters (based on assumptions)
  a_post_1 <- a_prior + x1
  b_post_1 <- b_prior + (n1 - x1)
  
  a_post_2 <- a_prior + x2
  b_post_2 <- b_prior + (n2 - x2)
  
  # Generate posterior samples using Beta distribution
  samples_1 <- rbeta(10000, a_post_1, b_post_1)
  samples_2 <- rbeta(10000, a_post_2, b_post_2)
  
  # Compute probabilities P(π > 0.4)
  p_greater_0.4_1 <- mean(samples_1 > 0.4)
  p_greater_0.4_2 <- mean(samples_2 > 0.4)
  return(c(p_greater_0.4_1,p_greater_0.4_2))
}
run_simu(0.7, 0.3)

## Running for all 9 values
p1_s = c(0.3, 0.5, 0.7)
p2_s = c(0.3, 0.5, 0.7)

outs1 = numeric(9)
outs2 = numeric(9)
count = 1
for(a in p1_s)
{
 for(b in p2_s)
 {
  out = run_simu(a,b)
  outs1[count] = out[1]
  outs2[count] = out[2]
  count = count + 1
 }
}
outs1 = data.frame(matrix(outs1, 3, 3))
colnames(outs1) = c(0.3, 0.5, 0.7)
rownames(outs1) = c(0.3, 0.5, 0.7)
outs2 = data.frame(matrix(outs2, 3, 3))
colnames(outs2) = c(0.3, 0.5, 0.7)
rownames(outs2) = c(0.3, 0.5, 0.7)
