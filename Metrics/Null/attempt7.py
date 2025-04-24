import numpy as np

import os



# File path template

file_template = '/project/cstr/Jaiman/sp25/VAPORdata/Test2/VAPOR/Merge/B3D_MHD.{:03d}'



# Output file path

output_file_path = 'output.csv'



# Array dimensions and size

nx, ny, nz = 361, 361, 361

array_size_bytes = nx * ny * nz * 8

header_size_bytes = 4



# Time step range

time_steps = range(1, 160)



# Store all k-lists to determine the max number of sign changes

all_k_lists = []



for t in time_steps:

    file_path = file_template.format(t)

    if not os.path.exists(file_path):

        print(f"Warning: File '{file_path}' not found. Skipping...")

        all_k_lists.append([])  # Pad with empty list

        continue



    try:

        with open(file_path, 'rb') as f:

            f.seek(header_size_bytes)

            f.seek(array_size_bytes, 1)  # Skip first array

            second_array = np.fromfile(f, dtype=np.float64, count=nx * ny * nz)

            second_array = second_array.reshape((nx, ny, nz), order='F')

            values = second_array[180, 180, :]



        sign_change_ks = [

            k for k in range(1, nz)

            if np.sign(values[k]) != np.sign(values[k - 1])

            and values[k] != 0.0 and values[k - 1] != 0.0

        ]



        all_k_lists.append(sign_change_ks)

        print(f"TimeStep {t}: Sign changes at k = {','.join(map(str, sign_change_ks)) if sign_change_ks else 'None'}")



    except Exception as e:

        print(f"Error processing file '{file_path}': {e}")

        all_k_lists.append([])



# Determine the max number of sign changes across all timesteps

max_k_count = max(len(k_list) for k_list in all_k_lists)



# Write to CSV

with open(output_file_path, 'w') as f_out:

    # Header

    headers = ['TimeStep'] + [f'k{i+1}' for i in range(max_k_count)]

    f_out.write(','.join(headers) + '\n')



    for t_idx, k_list in enumerate(all_k_lists):

        padded_k_list = k_list + [''] * (max_k_count - len(k_list))  # Pad with empty strings

        row = [str(time_steps[t_idx])] + [str(k) for k in padded_k_list]

        f_out.write(','.join(row) + '\n')



print(f"Sign change CSV saved to '{output_file_path}' successfully.")


