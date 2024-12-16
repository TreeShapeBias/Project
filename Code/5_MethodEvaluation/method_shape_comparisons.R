library(ape)
library(treebalance)
library(dplyr)
library(ggplot2)

shape_calc <- function(tree) {
  shape <- collessI(tree, method="corrected")
  return(shape)
}

data_table <- read.csv("./data_table.csv")

model_data <- data.frame(
  true_tree = numeric(nrow(data_table)* 4),
  inferred_tree = numeric(nrow(data_table)* 4),
  inference_method = numeric(nrow(data_table)* 4)
)

for (i in seq_len(nrow(data_table))) {
  true_tree <- read.tree(data_table$true_tree[i])
  model_data$true_tree[(i*4 - 3):(i*4)] <- shape_calc(true_tree)
  
  treeqmc_tree <- read.tree(data_table$tree.qmc[i])
  model_data$inferred_tree[i*4 - 3] <- shape_calc(treeqmc_tree)
  model_data$inference_method[i*4 - 3] <- "Tree-QMC"
  
  astral.iii_tree <- read.tree(data_table$astral.iii[i])
  model_data$inferred_tree[i*4 - 2] <- shape_calc(astral.iii_tree)
  model_data$inference_method[i*4 - 2] <- "ASTRAL-III"
  
  caster_tree <- read.tree(data_table$caster[i])
  model_data$inferred_tree[i*4 - 1] <- shape_calc(caster_tree)
  model_data$inference_method[i*4 - 1] <- "CASTER"
  
  astral.iv_tree <- read.tree(data_table$astral.iv[i])
  model_data$inferred_tree[i*4] <- shape_calc(astral.iv_tree)
  model_data$inference_method[i*4] <- "ASTRAL-IV"
}


ggplot(data = model_data, mapping=aes(true_tree, inferred_tree)) + geom_point(mapping = aes(color = inference_method)) + 
  geom_abline(slope=1, intercept=0, color="black", linewidth=1, linetype="dashed") + ylim(0, 1) + xlim(0, 1) + 
  ylab("Inferred Tree") + xlab("True Tree") + theme_bw() + guides(color=guide_legend(title="Inference Method")) +
  geom_smooth(mapping = aes(color=inference_method), method = "lm")
