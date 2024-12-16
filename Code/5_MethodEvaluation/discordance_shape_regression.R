library(rstan)
library(ape)
library(TreeDist)
library(treebalance)
library(dplyr)
library(bayesplot)
library(plotly)
library(boot)
library(ggplot2)


shape_calc <- function(tree) {
  shape <- collessI(tree, method="corrected")
  return(shape)
}

data_table <- read.csv("./data_table.csv")

model_data <- data.frame(
  shape = numeric(nrow(data_table)),
  gene_discordance = numeric(nrow(data_table))
)
for (i in seq_len(nrow(data_table))) {
  true_tree <- read.tree(data_table$true_tree[i])
  model_data$shape[i] <- shape_calc(true_tree)
  model_data$gene_discordance[i] <- data_table$gene_tree_error[i]
}

model <- "
  data {
    int<lower=0> N;
    vector<lower=0, upper=1>[N] shape;
    vector<lower=0, upper=1>[N] gene_discordance;
  }
  parameters {
    real intercept;
    real slope;
    real<lower=0> phi;
  }
  transformed parameters {
    vector[N] mu;
    vector[N] a;
    vector[N] b;
  
    mu = inv_logit(
        intercept + slope .* shape
    );
  
    a = mu .* phi;
    b = (1 - mu) .* phi;
  }
  model {
    intercept ~ normal(0, 1);
    slope ~ normal(0, 1);
    phi ~ gamma(5,1);
  
    gene_discordance ~ beta(a, b);
  }
  "

stan_data = list(N=nrow(model_data), shape=model_data$shape, 
                 gene_discordance=model_data$gene_discordance)

stan_fit <- stan(model_code = model, data = stan_data, iter = 2000, chains = 4)

posterior <- as.data.frame(stan_fit)
mcmc_areas(posterior, 
           pars = c("intercept", "slope"),
           prob = 0.95) + labs(title = "Posterior Distribution")


x_grid <- seq(0, 1, length.out = 100)
intercept_samples <- posterior$intercept
slope_samples <- posterior$slope
predictions <- expand.grid(
  x = x_grid,
  draw = 1:nrow(posterior)
) %>%
  mutate(
    intercept = rep(intercept_samples, each = length(x_grid)),
    slope = rep(slope_samples, each = length(x_grid)),
    predicted_y = inv.logit(intercept + slope * x)
  )

summary_predictions <- predictions %>%
  group_by(x) %>%
  summarize(
    mean_y = mean(predicted_y),
    lower = quantile(predicted_y, 0.025),
    upper = quantile(predicted_y, 0.975),
    .groups = "drop"
  )


ggplot(summary_predictions, aes(x = x)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(aes(y = mean_y), size = 1) +
  geom_point(data = model_data, aes(x = shape, y = gene_discordance), color = "black", alpha = 0.6) +
  labs(
    x = "Tree Shape", y = "Gene Tree Discordance"
  ) +
  theme_minimal()
