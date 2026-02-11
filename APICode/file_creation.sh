#!/bin/bash

# List of variables
#variables=("BX" "BY" "BZ" "CB2" "CT_BT" "VX" "VY" "VZ")
variables=("RO")

# Number of timesteps
start_time=1
end_time=50

# File size and format parameters
data_size="361 361 361"
brick_size="360. 360. 360."
data_format="DOUBLE"
brick_origin="0. 0. 0."
byte_offset=4

# Loop through each variable
for var in "${variables[@]}"; do
    # Loop through each time step
    for ((t=$start_time; t<=$end_time; t++)); do
        # Format timestep with leading zeros (e.g., 001, 002, ..., 120)
        timestep=$(printf "%03d" $t)
        data_file="/project/cstr/Jaiman/sp26/Project1/Test2/VAPOR/B3D.001.${var}.${timestep}"
        bov_file="/project/cstr/Jaiman/sp26/Project1/Test2/API/BOV/${var}_${timestep}.bov"

        # Write the BOV file
        cat > $bov_file <<EOL
# TIME is a floating point value specifying the timestep being read in DATA_FILE

TIME: $t

# DATA_FILE points to a binary data file.  It can be a full file path, or a path relative to the BOV header.

DATA_FILE: $data_file

# The data size corresponds to NX,NY,NZ in the above example code.  It must contain three values

DATA_SIZE: $data_size

BRICK_SIZE: $brick_size

# Allowable values for DATA_FORMAT are: INT,FLOAT,DOUBLE

DATA_FORMAT: $data_format

# VARIABLE is a string that specifies the variable being read in DATA_FILE.  Must be alphanumeric (abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890_-)

VARIABLE: $var

# BRICK_ORIGIN lets you specify a new coordinate system origin for
# the mesh that will be created to suit your data.  It must contain three values.

BRICK_ORIGIN: $brick_origin

# Fortran has 4 byte header and 4 byte footer
BYTE_OFFSET: $byte_offset
EOL

       # echo "Generated $bov_file"
    done
done

echo "All BOV files generated successfully."

