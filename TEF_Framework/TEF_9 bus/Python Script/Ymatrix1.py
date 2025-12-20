import psse3603
import psspy
import os
import re
import numpy as np
import scipy.io as sio

# ==============================================================================
# 1. SETUP
# ==============================================================================
work_dir = r"C:\Users\Acer\Desktop\final year project\energy function\My_Work\WSCC_Programs_MY"
result_dir=r"C:\Users\Acer\Desktop\final year project\energy function\My_Work\WSCC_Programs_MY\Result_Directory"
raw_file = os.path.join(work_dir, "IEEE9bus.raw")
txt_file = os.path.join(result_dir, "Ybus_Export.txt")
mat_file = os.path.join(result_dir, "Ybus_Export.mat")

psspy.psseinit(50)
psspy.read(0, raw_file)

# ==============================================================================
# 2. POWER FLOW
# ==============================================================================
print("Solving Power Flow...")
ierr = psspy.fnsl([0,0,0,1,1,0,99,0])
print(f"fnsl returned {ierr}")

# ==============================================================================
# 3. CONVERT
# ==============================================================================
print("Converting Generators...")
psspy.cong(0)

print("Converting Loads...")
psspy.conl(0,1,1,[0,0],[100,0,0,100])
psspy.conl(0,1,2,[0,0],[100,0,0,100])
psspy.conl(0,1,3,[0,0],[100,0,0,100])

psspy.ordr(0)
psspy.fact()
psspy.tysl(0)

# ==============================================================================
# 4. EXPORT Y-BUS AS TEXT
# ==============================================================================
print("Exporting Y-Bus...")
ierr_out = psspy.output_y_matrix(0,1,0,0,txt_file)
print(f"output_y_matrix returned {ierr_out}")

if ierr_out != 0:
    print("Error exporting Y-bus.")
    exit()

# ==============================================================================
# 5. PARSE THE TEXT FILE AND BUILD Y-BUS MATRIX
# ==============================================================================
print("Converting text file → numeric Y-bus...")

Yreal = {}
Yimag = {}
max_bus = 0

with open(txt_file, "r") as f:
    for line in f:
        parts = re.findall(r"[-+]?\d*\.\d+|\d+", line)
        if len(parts) == 4:
            i = int(parts[0])
            j = int(parts[1])
            real = float(parts[2])
            imag = float(parts[3])

            max_bus = max(max_bus, i, j)

            Yreal[(i, j)] = real
            Yimag[(i, j)] = imag

# Build matrix
Y = np.zeros((max_bus, max_bus), dtype=complex)
for (i, j), val in Yreal.items():
    Y[i-1, j-1] = val + 1j * Yimag[(i, j)]

# ==============================================================================
# 6. SAVE TO .MAT FILE
# ==============================================================================
print("Saving Y-bus to MATLAB format...")

sio.savemat(mat_file, {"Ybus": Y})

print(f"SUCCESS: Ybus.mat created at:\n   {mat_file}")
