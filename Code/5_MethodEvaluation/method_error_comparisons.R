library(rstan)
library(ape)
library(TreeDist)
library(treebalance)
library(dplyr)
library(phangorn)
library(bayesplot)
library(plotly)

shape_calc <- function(tree) {
  shape <- collessI(tree, method="corrected")
  return(shape)
}

normalized_rf <- function(tree1, tree2) {
  rf <- RF.dist(tree1, tree2, normalize = TRUE)
  return(rf)
}

run_analysis <- function(model_data, method_string){
  correlation_model <- "
  data {
    int<lower=0> N;
    vector<lower=0, upper=1>[N] rf_error;
    vector<lower=0, upper=1>[N] shape;
    vector<lower=0, upper=1>[N] gene_discordance;
  }
  parameters {
    real intercept;
    real shape_slope;
    real discordance_slope;
    real<lower=0> phi;
  }
  transformed parameters {
    vector[N] mu;
    vector[N] a;
    vector[N] b;
  
    mu = inv_logit(
        intercept + shape_slope .* shape +
        discordance_slope .* gene_discordance
    );
  
    a = mu .* phi;
    b = (1 - mu) .* phi;
  }
  model {
    intercept ~ normal(0, 1);
    discordance_slope ~ normal(0, 1);
    shape_slope ~ normal(0, 1);
    phi ~ gamma(5,1);
  
    rf_error ~ beta(a, b);
  }
  "
  
  stan_data = list(N=nrow(model_data), rf_error=model_data$rf_error, 
                     shape=model_data$shape, gene_discordance=model_data$gene_discordance)
  
  stan_fit <- stan(model_code = correlation_model, data = stan_data, iter = 2000, chains = 4)
  
  posterior <- as.matrix(stan_fit)
  mcmc_areas(posterior, 
             pars = c("intercept", "shape_slope", "discordance_slope"),
             prob = 0.95) + labs(title = paste(method_string, "Posterior Distribution"))
  
  x <- seq(0, 1, length.out = 100)
  gene_error <- seq(0, 1, length.out = 100)
  grid <- expand.grid(x = x, gene_error = gene_error)
  extracted_fit <- extract(stan_fit, pars=c("intercept", "shape_slope", "discordance_slope"))
  mean_intercept <- mean(extracted_fit$intercept)
  mean_shape_slope <- mean(extracted_fit$shape_slope)
  mean_discordance_slope <- mean(extracted_fit$discordance_slope)
  print(mean_intercept)
  print(mean_shape_slope)
  print(mean_discordance_slope)
  grid$mu <- plogis(mean_intercept + mean_shape_slope * grid$x + mean_discordance_slope * grid$gene_error)
  z_matrix <- matrix(grid$mu, nrow = length(x), ncol = length(gene_error))
  plot_ly(
    x = x, 
    y = gene_error, 
    z = z_matrix, 
    type = "surface",
    colorscale = "Viridis",
    opacity=0.5
  ) %>% add_markers(
    x = model_data$shape,
    y = model_data$gene_discordance,
    z = model_data$rf_error,
    marker = list(color = "red", size = 10),
    name = "Observed Points"
  ) %>%
    layout(
      scene = list(
        xaxis = list(title = "Normalized Colless' I"),
        yaxis = list(title = "Gene Tree Discordance"),
        zaxis = list(title = "Predicted Mean Error")
      ),
      title = paste("Mean Posterior RF Prediction (", method_string, ")", sep="")
    )
}



data_table <- read.csv("./data_table.csv")

# ASTRAL-III Stuff
astral_3_data <- data.frame(
  shape = numeric(nrow(data_table)),
  rf_error = numeric(nrow(data_table)),
  gene_discordance = numeric(nrow(data_table))
)
for (i in seq_len(nrow(data_table))) {
  true_tree <- read.tree(data_table$true_tree[i])
  
  astral_tree <- read.tree(data_table$astral.iii[i])
  astral_tree$tip.label <- sub("_0_0$", "", astral_tree$tip.label)
  astral_3_data$shape[i] <- shape_calc(true_tree)
  astral_3_data$gene_discordance[i] <- data_table$gene_tree_error[i]
  astral_3_data$rf_error[i] <- normalized_rf(true_tree, astral_tree) + 1e-6
}
run_analysis(astral_3_data, "ASTRAL-III")


# Tree-QMC Stuff
treeqmc_data <- data.frame(
  shape = numeric(nrow(data_table)),
  rf_error = numeric(nrow(data_table)),
  gene_discordance = numeric(nrow(data_table))
)
for (i in seq_len(nrow(data_table))) {
  true_tree <- read.tree(data_table$true_tree[i])
  
  treeqmc_tree <- read.tree(data_table$tree.qmc[i])
  treeqmc_tree$tip.label <- sub("_0_0$", "", treeqmc_tree$tip.label)
  treeqmc_data$shape[i] <- shape_calc(true_tree)
  treeqmc_data$gene_discordance[i] <- data_table$gene_tree_error[i]
  treeqmc_data$rf_error[i] <- normalized_rf(true_tree, treeqmc_tree) + 1e-6
}
run_analysis(treeqmc_data, "Tree-QMC")

# CASTER Stuff
caster_data <- data.frame(
  shape = numeric(nrow(data_table)),
  rf_error = numeric(nrow(data_table)),
  gene_discordance = numeric(nrow(data_table))
)
for (i in seq_len(nrow(data_table))) {
  true_tree <- read.tree(data_table$true_tree[i])
  
  caster_tree <- read.tree(data_table$caster[i])
  caster_tree$tip.label <- sub("_0_0$", "", caster_tree$tip.label)
  caster_data$shape[i] <- shape_calc(true_tree)
  caster_data$gene_discordance[i] <- data_table$gene_tree_error[i]
  caster_data$rf_error[i] <- normalized_rf(true_tree, caster_tree) + 1e-6
}
run_analysis(caster_data, "CASTER")


# ASTRID Stuff
astrid_data <- data.frame(
  shape = numeric(nrow(data_table)),
  rf_error = numeric(nrow(data_table)),
  gene_discordance = numeric(nrow(data_table))
)
for (i in seq_len(nrow(data_table))) {
  true_tree <- read.tree(data_table$true_tree[i])
  
  astrid_tree <- read.tree(data_table$astrid[i])
  astrid_tree$tip.label <- sub("_0_0$", "", astrid_tree$tip.label)
  astrid_data$shape[i] <- shape_calc(true_tree)
  astrid_data$gene_discordance[i] <- data_table$gene_tree_error[i]
  astrid_data$rf_error[i] <- normalized_rf(true_tree, astrid_tree) + 1e-6
}
run_analysis(astrid_data, "ASTRID")

# ASTRAL-IV Stuff
astral4_data <- data.frame(
  shape = numeric(nrow(data_table)),
  rf_error = numeric(nrow(data_table)),
  gene_discordance = numeric(nrow(data_table))
)
for (i in seq_len(nrow(data_table))) {
  true_tree <- read.tree(data_table$true_tree[i])
  
  astral4_tree <- read.tree(data_table$astrid[i])
  astral4_tree$tip.label <- sub("_0_0$", "", astral4_tree$tip.label)
  astral4_data$shape[i] <- shape_calc(true_tree)
  astral4_data$gene_discordance[i] <- data_table$gene_tree_error[i]
  astral4_data$rf_error[i] <- normalized_rf(true_tree, astral4_tree) + 1e-6
}
run_analysis(astral4_data, "ASTRAL-IV")