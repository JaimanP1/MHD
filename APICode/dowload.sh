#!/bin/bash

# If you want to download and use the data locally instead of using the API. Adjust ts and file_num parameters as needed

# Output file
output_file="commands.txt"
> "$output_file"  # Clear the file if it exists

# List of variables
vars=("BX" "BY" "BZ" "CB2" "CT_BT" "VX" "VY" "VZ")

echo "vdccreate -dimension 361x361x361 -numts 30 -ncvars3d vx:vy:vz:cb2:ct_bt:bx:by:bz test.vdc" >> "$output_file"

# Loop through variables
for var in "${vars[@]}"; do
    for ((ts=0; ts<=30; ts++)); do
        file_num=$(printf "%03d" $((ts + 50)))
        echo "raw2vdc -ts $ts -varname ${var,,} test.vdc B3D.001.${var}.R.${file_num}" >> "$output_file"
    done
done

