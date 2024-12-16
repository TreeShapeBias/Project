# main script for generating tree pool and selecting trees.
library(ape)
library(poweRbal)
library(dplyr)
library(stringr)
library(treebalance)

set.seed(4932)

# assign branch length using exponential distribution
assign_yule_branch_lengths <- function(tree, speciation_rate = 1.0) {
  if (is.null(tree$edge.length)) {
    branch_lengths <- rexp(n = nrow(tree$edge), rate = speciation_rate)
    tree$edge.length <- branch_lengths
  }
  return(tree)
}

make_ultrametric <- function(tree, desired_length = 1, lambda = 0.5) {
  ultrametric_tree <- tryCatch({
    chronos(tree, lambda = lambda)
  }, error = function(e) {
    message("chronos error: ", e$message)
    return(NULL)
  })
  
  if (is.null(ultrametric_tree)) {
    return(NULL)
  }
  
  if (!is.ultrametric(ultrametric_tree)) {
    message("chronos failed to produce an ultrametric tree.")
    return(NULL)
  }
  
  # Rescale the tree to the desired length
  current_max_length <- max(node.depth.edgelength(ultrametric_tree))
  scaling_factor <- desired_length / current_max_length
  ultrametric_tree$edge.length <- ultrametric_tree$edge.length * scaling_factor
  
  return(ultrametric_tree)
}

calculate_colless <- function(tree) {
  tryCatch({
    collessI(tree, method = "corrected")
  }, error = function(e) {
    NA
  })
}

generate_tree_pool <- function(pool_simphy_folder, num_generated_per_model = 100) {
  all_species_trees <- list()
  subfolder_map <- c()  # To track the origin of each tree
  
  subfolders <- list.dirs(pool_simphy_folder, recursive = FALSE)
  
  for (i in seq_along(subfolders)) {
    subfolder <- subfolders[i]
    species_tree_file <- file.path(subfolder, "s_tree.trees")
    
    if (file.exists(species_tree_file)) {
      tree <- tryCatch({
        read.tree(species_tree_file)
      }, error = function(e) {
        message("Error reading tree from ", species_tree_file, ": ", e$message)
        NULL
      })
      
      if (!is.null(tree)) {
        all_species_trees <- c(all_species_trees, list(tree))
        subfolder_map <- c(subfolder_map, paste0("pool_simphy_", i))
      }
    }
  }

models <- c("PDA", "ETM", "FordsAlpha", "Grow")
  
for (model in models) {
  for (j in 1:num_generated_per_model) {
    tree <- tryCatch({
      switch(model,
             "PDA" = poweRbal::genPDATree(100),
             "ETM" = poweRbal::genETMTree(100),
             "FordsAlpha" = poweRbal::genFordsAlphaTree(100,0.87),
             "Grow" = poweRbal::genGrowTree(100, use_built_in = "IF_sym", ZETA = 2))
    }, error = function(e) {
      message("Error generating tree with model ", model, ": ", e$message)
      NULL
    })
    
    if (!is.null(tree)) {
      all_species_trees <- c(all_species_trees, list(tree))
      subfolder_map <- c(subfolder_map, paste0("generated_", model, "_", j))
    }
  }
}

  return(list(trees = all_species_trees, map = subfolder_map))
}

categorize_trees <- function(all_species_trees) {
  colless_indices <- sapply(all_species_trees, calculate_colless)
  
  tree_data <- data.frame(
    TreeIndex = 1:length(all_species_trees),
    CollessIndex = colless_indices,
    stringsAsFactors = FALSE
  )
  
  # Remove trees with NA Colless index
  tree_data <- tree_data %>% filter(!is.na(CollessIndex))
  
  # Define Colless index ranges
  tree_data <- tree_data %>%
    mutate(Range = case_when(
      CollessIndex >= 0.0 & CollessIndex < 0.2 ~ "0.0-0.2",
      CollessIndex >= 0.2 & CollessIndex < 0.4 ~ "0.2-0.4",
      CollessIndex >= 0.4 & CollessIndex < 0.6 ~ "0.4-0.6",
      CollessIndex >= 0.6 & CollessIndex < 0.8 ~ "0.6-0.8",
      CollessIndex >= 0.8 & CollessIndex <= 1.0 ~ "0.8-1.0",
      TRUE ~ "Other"
    )) %>%
    filter(Range != "Other")
  
  return(tree_data)
}

select_trees <- function(tree_data, all_species_trees, subfolder_map, desired_trees_per_range = 4,
                         speciation_rate = 0.000001, desired_length = 2000000) {
  
  selected_trees <- list()
  selected_metadata <- data.frame()
  
  desired_ranges <- c("0.0-0.2","0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0")
  
  for (range in desired_ranges) {
    message("Selecting trees for range: ", range)
    
    available_trees <- tree_data %>% filter(Range == range)
    
    if (nrow(available_trees) == 0) {
      warning("No trees available in range ", range, ". Skipping this range.")
      next
    }
    
    available_indices <- available_trees$TreeIndex
    
    count <- 0
    
    while (count < desired_trees_per_range && length(available_indices) > 0) {

      selected_idx <- sample(available_indices, 1)
      tree <- all_species_trees[[selected_idx]]
      origin <- subfolder_map[selected_idx]
      tree <- assign_yule_branch_lengths(tree, speciation_rate)
      if (!is.ultrametric(tree)) {
        ultrametric_tree <- make_ultrametric(tree, desired_length = desired_length)
      
        if (is.null(ultrametric_tree)) {
          message("chronos failed for tree index ", selected_idx, ". Removing from pool.")
          available_indices <- available_indices[available_indices != selected_idx]
          next
        } else {
          tree <- ultrametric_tree
        }
      }
      
      # If successful, add to selected_trees
      selected_trees <- c(selected_trees, list(tree))
      selected_metadata <- rbind(selected_metadata, data.frame(
        TreeIndex = selected_idx,
        Range = range,
        Origin = origin,
        CollessIndex = tree_data$CollessIndex[tree_data$TreeIndex == selected_idx],
        stringsAsFactors = FALSE
      ))
      
      # Remove the selected index from available_indices to avoid reselection
      available_indices <- available_indices[available_indices != selected_idx]
      count <- count + 1
    }
    
    if (count < desired_trees_per_range) {
      warning("Only ", count, " trees could be selected for range ", range, 
              ". Desired: ", desired_trees_per_range, ".")
    }
  }
  
  return(list(selected_trees = selected_trees, metadata = selected_metadata))
}

save_selected_trees <- function(selected_trees, selected_metadata, output_base_folder) {
  dir.create(output_base_folder, showWarnings = FALSE, recursive = TRUE)
  
  for (i in 1:nrow(selected_metadata)) {
    tree <- selected_trees[[i]]
    range <- selected_metadata$Range[i]
    origin <- selected_metadata$Origin[i]
    tree_idx <- selected_metadata$TreeIndex[i]
    
    range_folder <- file.path(output_base_folder, range)
    dir.create(range_folder, showWarnings = FALSE, recursive = TRUE)
    
    if (str_starts(origin, "pool_simphy")) {
      # Extract subfolder number from origin
      subfolder_num <- str_extract(origin, "\\d+")
      filename <- paste0("s_tree_", subfolder_num, "_", tree_idx, ".trees")
    } else if (str_starts(origin, "generated")) {
      filename <- paste0(origin, ".trees")
    } else {
      filename <- paste0("tree_", tree_idx, ".trees")
    }
    
    # Save the tree
    write.tree(tree, file = file.path(range_folder, filename))
  }
  summary_csv_path <- file.path(output_base_folder, "selected_trees_summary.csv")
  write.csv(selected_metadata, summary_csv_path, row.names = FALSE)
  
  message("Selected trees and summary have been saved to ", output_base_folder)
}

# Define folder paths
pool_simphy_folder <- "path_to_trees_generated_using_SimPhy"  # Folder containing original species trees
selected_trees_output_folder <- "path_to_folder"  # Folder to save final selected trees

# Step 1: Generate the tree pool
message("=== Generating Tree Pool ===")
tree_pool <- generate_tree_pool(pool_simphy_folder)
all_species_trees <- tree_pool$trees
subfolder_map <- tree_pool$map

# Step 2: Calculate Colless index and categorize trees
message("=== Calculating Colless Index and Categorizing Trees ===")
tree_data <- categorize_trees(all_species_trees)

# Step 3: Select trees based on Colless index ranges with error handling
message("=== Selecting Trees Based on Colless Index Ranges ===")
selection_result <- select_trees(
  tree_data = tree_data,
  all_species_trees = all_species_trees,
  subfolder_map = subfolder_map,
  desired_trees_per_range = 4,
  speciation_rate = 0.000001,
  desired_length = 2000000
)

selected_trees <- selection_result$selected_trees
selected_metadata <- selection_result$metadata

# Step 4: Save selected trees and summary
message("=== Saving Selected Trees and Summary ===")
save_selected_trees(
  selected_trees = selected_trees,
  selected_metadata = selected_metadata,
  output_base_folder = selected_trees_output_folder
)
