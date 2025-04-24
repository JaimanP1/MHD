import numpy as np

import os



# File path template

file_template = '/project/cstr/Jaiman/sp25/Test2/VAPOR/Merge/B3D_MHD.{:03d}'



# Output file path

output_file_path = 'output.txt'



# Array dimensions and size

nx, ny, nz = 361, 361, 361

array_size_bytes = nx * ny * nz * 8  # Size of one array in bytes

header_size_bytes = 4



# Time step range

time_steps = range(1, 160)  # 001 to 160



# Prepare to write the output file

with open(output_file_path, 'w') as f_out:

    f_out.write("# TimeStep\tk\tValue\n")



    # Iterate over all time steps

    for t in time_steps:

        file_path = file_template.format(t)

        

        # Check if the file exists

        if not os.path.exists(file_path):

            print(f"Warning: File '{file_path}' not found. Skipping...")

            continue



        try:

            # Open the binary file and read the data

            with open(file_path, 'rb') as f:

                # Skip the 4-byte header

                f.seek(header_size_bytes)

                

                # Skip the first array

                f.seek(array_size_bytes, 1)  # Move forward by 1 array

                

                # Read the second array

                second_array = np.fromfile(f, dtype=np.float64, count=nx * ny * nz)

                

                # Reshape to 3D array, assuming column-major order (Fortran order)

                second_array = second_array.reshape((nx, ny, nz), order='F')

                

                # Extract values at (180, 180, k) for all k

                values = second_array[180, 180, :]



            # Find the first occurrence of a sign change

            smallest_k = -1

            smallest_value = 0.0

            for k in range(1, nz):

                if np.sign(values[k]) != np.sign(values[k - 1]) and values[k] != 0.0 and values[k - 1] != 0.0:

                    smallest_k = k

                    smallest_value = values[k]

                    break



            # Record the result if a sign change was found

            if smallest_k != -1:

                f_out.write(f"{t}\t{smallest_k}\t{smallest_value:.6e}\n")

                print(f"Sign change at timestep {t}: k = {smallest_k}, value = {smallest_value:.6e}")

            else:

                print(f"No sign change found for timestep {t}.")

        

        except Exception as e:

            print(f"Error processing file '{file_path}': {e}")



print(f"Sign change summary saved to '{output_file_path}' successfully.")


