import os
import sys
from ete3 import Tree

# Input arguments
input_folder = sys.argv[1]
reference_tree_file = os.path.join(input_folder, "s_tree.trees")


try:
    with open(reference_tree_file, 'r') as f:
        reference_tree = Tree(f.readline().strip(), format=1)
    print(f"Reference tree loaded from: {reference_tree_file}")
except FileNotFoundError:
    print(f"Error: Reference tree file {reference_tree_file} not found.")
    sys.exit(1)

for leaf in reference_tree:
    leaf.name = leaf.name + "_0_0"

rf_distances = []
for filename in os.listdir(input_folder):
    if filename.startswith("g_tree"):
        file_path = os.path.join(input_folder, filename)
        try:
            with open(file_path, 'r') as f:
                gene_tree = Tree(f.readline().strip(), format=1)
            # Compute normalized RF distance
            rf_distance, max_rf, *_ = reference_tree.robinson_foulds(gene_tree, unrooted_trees=True)
            normalized_rf = rf_distance / max_rf if max_rf > 0 else 0.0
            rf_distances.append(normalized_rf)
        except Exception as e:
            print(f"Error processing {file_path}: {e}")

if rf_distances:
    average_rf = sum(rf_distances) / len(rf_distances)
    print(f"Average RF error: {average_rf:.4f}")
else:
    print("No valid gene trees were found to compute RF error.")
