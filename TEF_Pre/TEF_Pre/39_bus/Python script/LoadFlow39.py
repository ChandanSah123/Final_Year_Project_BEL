import psse3603
import psspy
import os
import pandas as pd
import math  # Added for conversion

# ==============================================================================
# 1. SETUP & CONFIGURATION
# ==============================================================================
work_dir = r"C:\Users\Acer\Desktop\final year project\energy function\TEF_Framework\TEF_Pre\39_bus\Python Script"
result_dir = r"C:\Users\Acer\Desktop\final year project\energy function\TEF_Framework\TEF_Pre\39_bus\Matlab script"
raw_file = os.path.join(work_dir, "ieee39bus1")
load_flow_file = os.path.join(result_dir, "LoadFlowResults.xlsx")

# Initialize PSS/E
psspy.psseinit(50)

# ==============================================================================
# 2. LOAD CASE & SOLVE
# ==============================================================================
print("Loading Case...")
psspy.read(0, raw_file)

print("Solving Load Flow...")
ierr_fnsl = psspy.fnsl([0, 0, 0, 1, 1, 0, 99, 0])

if ierr_fnsl > 1:
    print(f"Warning: Load flow solution code: {ierr_fnsl}")
else:
    print(" > Load Flow Solved Successfully.")

# ==============================================================================
# 3. EXTRACT BUS DATA (AND CONVERT RAD -> DEG)
# ==============================================================================
print("Extracting Bus Data...")
# Arguments: SID(-1), Flag(2), StringList
ierr_num, bus_nums = psspy.abusint(-1, 2, ['NUMBER'])
ierr_real, bus_data = psspy.abusreal(-1, 2, ['PU', 'ANGLE'])

bus_results = []
if ierr_num == 0:
    for num, pu, ang_rad in zip(bus_nums[0], bus_data[0], bus_data[1]):
        # --- FIX: CONVERT RADIANS TO DEGREES ---
        ang_deg = ang_rad * (180.0 / math.pi)
        bus_results.append([num, pu, ang_deg])
    print(f" > Extracted {len(bus_results)} buses.")
else:
    print(f"ERROR: Bus extraction failed (Code {ierr_num})")

# ==============================================================================
# 4. EXTRACT BRANCH DATA (FIXED ARGUMENTS)
# ==============================================================================
print("Extracting Branch Data...")
line_results = []

# Arguments for PSS/E 36 abrn___: 
# SID(-1), Owner(0), Ties(0), Type(3), Entry(2), StringList
# Type 3 = Select BOTH Lines and Transformers
try:
    l_int_err, l_ints = psspy.abrnint(-1, 0, 0, 3, 2, ['FROMNUMBER', 'TONUMBER'])
    l_chr_err, l_chrs = psspy.abrnchar(-1, 0, 0, 3, 2, ['ID'])
    l_real_err, l_flw = psspy.abrnreal(-1, 0, 0, 3, 2, ['P', 'Q'])

    if l_int_err == 0 and l_chr_err == 0 and l_real_err == 0:
        for frm, to, ckt, p, q in zip(l_ints[0], l_ints[1], l_chrs[0], l_flw[0], l_flw[1]):
            line_results.append([frm, to, ckt, p, q])
        print(f" > Extracted {len(line_results)} branches.")
    else:
        print(f"ERROR: Branch extraction failed. Codes: {l_int_err}")

except TypeError:
    print("CRITICAL ERROR: PSS/E version mismatch on arguments. Trying fallback (4 ints)...")
    # Fallback for older PSS/E versions just in case
    try:
        l_int_err, l_ints = psspy.abrnint(-1, 0, 0, 2, ['FROMNUMBER', 'TONUMBER'])
        l_chr_err, l_chrs = psspy.abrnchar(-1, 0, 0, 2, ['ID'])
        l_real_err, l_flw = psspy.abrnreal(-1, 0, 0, 2, ['P', 'Q'])
        if l_int_err == 0:
             for frm, to, ckt, p, q in zip(l_ints[0], l_ints[1], l_chrs[0], l_flw[0], l_flw[1]):
                line_results.append([frm, to, ckt, p, q])
    except:
        print("Fallback failed.")

# ==============================================================================
# 5. SAVE TO EXCEL
# ==============================================================================
try:
    with pd.ExcelWriter(load_flow_file) as writer:
        if bus_results:
            df_bus = pd.DataFrame(bus_results, columns=['Bus Number', 'Voltage (pu)', 'Angle (deg)'])
            df_bus.to_excel(writer, sheet_name='Bus Data', index=False)
        
        if line_results:
            df_line = pd.DataFrame(line_results, columns=['From Bus', 'To Bus', 'ID', 'P (MW)', 'Q (MVar)'])
            df_line.to_excel(writer, sheet_name='Line Flows', index=False)

    print(f"SUCCESS: Results saved to {load_flow_file}")

except Exception as e:
    print(f"CRITICAL ERROR: Could not save Excel file. Is it open? \nDetails: {e}")

print("--- Script Finished ---")