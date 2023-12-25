import os
import yaml

# Step 1: Find all .wdl files in the current directory and its subdirectories
wdl_files = [os.path.join(root, file)
             for root, dirs, files in os.walk(".")
             for file in files
             if file.endswith(".wdl")]

# Step 2: Initialize an empty list to store the names of the WDL files
tools = []

# Step 3: Iterate over the list of WDL files
for wdl_file in wdl_files:
    # Step 4: For each WDL file, add a new entry to the tools list
    tools.append({"name": wdl_file})

# Create the .dockstore.yml file
dockstore_yml = {
    "version": 1.2,
    "tools": tools
}

# Write the .dockstore.yml file
with open(".dockstore.yml", "w") as outfile:
    yaml.dump(dockstore_yml, outfile, default_flow_style=False)