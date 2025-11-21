import psse3603
"""
extract_ybus.py
Extracts full Y-bus from a PSSE .sav case using psspy.y_matrix()
Produces:
 - ybus_dense.csv  (real and imag combined as complex strings)
 - ybus_real.csv   (real part)
 - ybus_imag.csv   (imag part)
 - ybus.xlsx       (Excel with complex strings)
"""

import os
import sys
import numpy as np
import pandas as pd

# ---- Adjust these to your environment ----
CASE_FILE = r"C:\path\to\Converted_ieee9bus.sav"   # <- change to your .sav path
OUTPUT_DIR = os.getcwd()                           # or specify folder
EXCEL_FILENAME = os.path.join(OUTPUT_DIR, "ybus.xlsx")
CSV_COMPLEX = os.path.join(OUTPUT_DIR, "ybus_dense.csv")
CSV_REAL = os.path.join(OUTPUT_DIR, "ybus_real.csv")
CSV_IMAG = os.path.join(OUTPUT_DIR, "ybus_imag.csv")

# ---- PSSE imports ----
import psspy

# initialize PSSE
psspy.psseinit()

# load case
ierr = psspy.case(CASE_FILE)
if ierr != 0:
    raise RuntimeError(f"psspy.case() returned ierr={ierr}. Check file path and PSSE environment.")

# run a powerflow if needed
ierr = psspy.fnsl([0, 0, 0, 0, 0, 0, 0, 0])
if ierr != 0:
    print(f"Warning: psspy.fnsl returned ierr={ierr}. Proceeding anyway.")

# get sparse Y representation (rows, cols, G, B)
ierr, rows, cols, G, B = psspy.y_matrix()
if ierr != 0:
    raise RuntimeError(f"psspy.y_matrix() returned ierr={ierr}. Cannot extract Y-bus.")

rows = list(rows)
cols = list(cols)
G = list(G)
B = list(B)

# build unique bus ID set
unique_ids = sorted(set(rows) | set(cols))

# map bus IDs to dense indices
ierr, nbus_est = psspy.abuscount()
if ierr != 0:
    nbus_est = None

use_index_mapping = True
if nbus_est is not None:
    if max(unique_ids) <= nbus_est and min(unique_ids) >= 1:
        use_index_mapping = False

if not use_index_mapping:
    buslist = unique_ids
    size = len(buslist)
    id_to_idx_map = {bus: idx for idx, bus in enumerate(buslist)}
    def id_to_idx(i):
        return id_to_idx_map[i]
else:
    max_index = nbus_est if nbus_est is not None else max(unique_ids)
    size = max_index
    def id_to_idx(i): return int(i) - 1

# build dense Y matrix
Y = np.zeros((size, size), dtype=complex)

for r, c, g, b in zip(rows, cols, G, B):
    ri = id_to_idx(r)
    ci = id_to_idx(c)
    val = complex(g, b)
    Y[ri, ci] += val
    if ri != ci:
        Y[ci, ri] += val  # ensure symmetry

# create labels for DataFrame
if use_index_mapping:
    index_to_bus = buslist
else:
    index_to_bus = list(range(1, size+1))

labels = [str(b) for b in index_to_bus]

# prepare complex string and real/imag matrices
complex_str = np.vectorize(lambda z: f"{z.real:.12g}{'+' if z.imag>=0 else '-'}{abs(z.imag):.12g}j")(Y)
real_mat = Y.real
imag_mat = Y.imag

df_complex = pd.DataFrame(complex_str, index=labels, columns=labels)
df_real = pd.DataFrame(real_mat, index=labels, columns=labels)
df_imag = pd.DataFrame(imag_mat, index=labels, columns=labels)

# save to CSV
df_complex.to_csv(CSV_COMPLEX)
df_real.to_csv(CSV_REAL)
df_imag.to_csv(CSV_IMAG)

# save to Excel
with pd.ExcelWriter(EXCEL_FILENAME, engine="openpyxl") as writer:
    df_complex.to_excel(writer, sheet_name="Y_complex")
    df_real.to_excel(writer, sheet_name="Y_real")
    df_imag.to_excel(writer, sheet_name="Y_imag")

print(f"Done. Dense Y-bus size = {Y.shape}")
print(f"Files written:\n - {CSV_COMPLEX}\n - {CSV_REAL}\n - {CSV_IMAG}\n - {EXCEL_FILENAME}")
