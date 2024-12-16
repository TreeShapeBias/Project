import os
import sys
import subprocess

input_folder = sys.argv[1]
combined_trees_file = input_folder + sys.argv[2]
astral_jar_path = "../../Astral/astral.5.7.8.jar"
output_species_tree = input_folder + sys.argv[3]


with open(combined_trees_file, 'w') as outfile:
	for filename in os.listdir(input_folder):
		if filename.startswith('g_tree'):
			file_path = os.path.join(input_folder, filename)
			with open(file_path, 'r') as infile:
				for line in infile:
					outfile.write(line.strip() + '\n')
print(f"Gene trees combined into: {combined_trees_file}")


try:
	subprocess.run(
		['java', '-jar', astral_jar_path, '-i', combined_trees_file, '-o', output_species_tree],
		check=True
	)
	print(f"Species tree saved to: {output_species_tree}")
except subprocess.CalledProcessError as e:
	print(f"Error running ASTRAL: {e}")
except FileNotFoundError:
	print("Java or ASTRAL JAR not found. Ensure both are correctly installed.")

