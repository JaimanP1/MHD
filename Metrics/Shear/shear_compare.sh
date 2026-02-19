#!/bin/bash

# Directory to save PNG images
output_dir="./TwistingMotion/Old"
mkdir -p "$output_dir"

# Loop over time steps 001 to 100
for t in $(seq -w 1 100); do
    for var in VX VY VZ; do
        input_file="/project/cstr/Jaiman/sp25/VAPORdata/Test4/Merge/BOTTOM_${var}.${t}"
        output_file="${output_dir}/${var}_${t}.png"

        gnuplot <<EOF
set terminal pngcairo size 800,600
set output "${output_file}"
set pm3d
unset surface
set view 0,0,1,1
splot "${input_file}" using 1:2:3 with lines
EOF

    done
done

