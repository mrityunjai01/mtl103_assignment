from ilp import gomory
import numpy as np
import os
eps = 1e-6
file_list = ['1.in', '2.in', '3.in', '4.in', '5.in', '6.in']
# '3.in', '4.in', '5.in', '6.in', '7.in', '8.in', '9.in', '10.in']
expected_outputs = {
    '1.in': [2, 0],
    '2.in': [2, 1],
    '3.in': [0, 0, 3, 0],
    '4.in': [0, 4, 0, 0],
    '5.in': [0, 4, 0, 0, 16, 0, 0, 1],
    '6.in': [0, 0, 0, 0, 21, 0,0,0,0,0,0,0,0,2,0,0,0,0]
}
for file in file_list:
    obtained_output = gomory(os.path.join('examples', file))
    for ind, output in enumerate(obtained_output):
        if (np.abs(output - expected_outputs[file][ind]) > eps):
            # pass
            print(f"file {file}, \nexpected output - {expected_outputs[file]} \nobtained_output - {obtained_output} \ndiffer at index {ind}, where {output} != {expected_outputs[file][ind]}")

# print(f"no output means theres no mismatch")
