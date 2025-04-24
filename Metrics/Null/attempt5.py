import numpy as np



# File path

file_path = '/project/cstr/Jaiman/sp25/Test2/VAPOR/Merge/B3D_MHD.091'

output_file_path = 'sign_changes2.txt'



# Array dimensions and size

nx, ny, nz = 361, 361, 361

array_size_bytes = nx * ny * nz * 8  # Size of one array in bytes

header_size_bytes = 4



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



# Find locations where sign changes

sign_changes = []

for k in range(1, nz):

    if np.sign(values[k]) != np.sign(values[k - 1]) and values[k] != 0.0 and values[k - 1] != 0.0:

        sign_changes.append((k, values[k]))



# Save the sign change information to a file

with open(output_file_path, 'w') as f_out:

    f_out.write("# k\tvalue\n")

    for k, v in sign_changes:

        f_out.write(f"{k}\t{v:.6e}\n")



print(f"Sign changes recorded in '{output_file_path}' successfully.")


