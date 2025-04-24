import numpy as np

import os



# File path

file_path = '/project/cstr/Jaiman/sp25/Test2/VAPOR/Merge/B3D_MHD.001'



# Output file

output_file_path = 'n_z_output.txt'



# Array dimensions and sizes

nx, ny, nz = 361, 361, 361

array_size = nx * ny * nz

array_size_bytes = array_size * 8  # float64 = 8 bytes

header_size_bytes = 4



# Index to extract from

ix, iy = 180, 180



# Check if file exists

if not os.path.exists(file_path):

    print(f"File '{file_path}' not found.")

else:

    try:

        with open(file_path, 'rb') as f:

            # Skip header

            f.seek(header_size_bytes)



            # Read Bx

            Bx = np.fromfile(f, dtype=np.float64, count=array_size)

            Bx = Bx.reshape((nx, ny, nz), order='F')



            # Read By

            By = np.fromfile(f, dtype=np.float64, count=array_size)

            By = By.reshape((nx, ny, nz), order='F')



        # Extract Bx, By along z at (180, 180)

        Bx_line = Bx[ix, iy, :]

        By_line = By[ix, iy, :]



        # Compute B(z)

        B = np.sqrt(Bx_line**2 + By_line**2)



        # Compute dB/dz using finite differences

        dB_dz = np.zeros(nz)

        for z in range(nz):

            if z == 0:

                dB_dz[z] = (B[z + 1] - B[z])

            elif z == nz - 1:

                dB_dz[z] = (B[z] - B[z - 1])

            else:

                dB_dz[z] = (B[z + 1] - B[z - 1]) / 2



        # Compute n(z)

        n = np.zeros(nz)

        for z in range(nz):

            if B[z] != 0:

                n[z] = -z / B[z] * dB_dz[z]

            else:

                n[z] = 0.0  # avoid division by zero



        # Write results

        with open(output_file_path, 'w') as f_out:

            f_out.write("# z\tn(z)\n")

            for z in range(nz):

                f_out.write(f"{z}\t{n[z]:.6e}\n")



        print(f"n(z) successfully written to '{output_file_path}'.")



    except Exception as e:

        print(f"Error processing file '{file_path}': {e}")


