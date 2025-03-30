import numpy as np



# Define parameters

filename = '/project/cstr/Jaiman/sp25/Test2/VAPOR/Merge/B3D_MHD.040'  # Replace with actual file name

nx, ny, nz = 361, 361, 361  # Replace with actual grid dimensions

num_vars = 6  # bx, by, bz, vx, vy, vz



# Open file in binary mode

with open(filename, 'rb') as f:

    var_data = []

    

    for i in range(num_vars):

        # Read 4-byte header (record length)

        f.read(4)

        

        # Read array data (double precision floats)

        arr = np.fromfile(f, dtype=np.float32, count=nx*ny*nz)

        

        # Read 4-byte trailer (record length)

        f.read(4)

        

        # Reshape to 3D array

        arr = arr.reshape((nx, ny, nz), order='F')  # Fortran order

        

        # Append to list

        var_data.append(arr)



# Extract arrays

bx_r, by_r, bz_r, vx_r, vy_r, vz_r = var_data



# Extract by_r(0,0,z)

by_r_180_180_z = by_r[180, 180, :]

#print("by_r(0,0,z):", by_r_0_0_z)



# Find the indices where by_r(0,0,z) = 0

z_indices = np.where(np.isclose(by_r_180_180_z, 0.0, atol=1e-8))[0]



# Print results

if len(z_indices) > 0:

    #print(f"by_r(0,0,z) = 0 at indices: {z_indices}")

    for z in z_indices:

        print(f"z = {z}, by_r(180,180,{z}) = {by_r_180_180_z[z]}")

else:

    print("No point found where by_r(180,180,z) = 0.")


