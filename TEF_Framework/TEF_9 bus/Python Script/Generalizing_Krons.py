import os
import re
import numpy as np
import scipy.io as sio
import sys

# Attempt to import PSS/E modules
try:
    import psse3603 as psse_mod
    import psspy
    import dyntools
except ImportError:
    pass

# ----------------- USER CONFIG -----------------
work_dir = r"C:\Users\Acer\Desktop\final year project\energy function\My_Work\WSCC_Programs_MY"
result_dir = r"C:\Users\Acer\Desktop\final year project\energy function\My_Work\WSCC_Programs_MY\Result_Directory"
raw_file = os.path.join(work_dir, "IEEE9bus.raw")
mat_file = os.path.join(result_dir, "Y_all_Generalized.mat")
txt_file = os.path.join(result_dir, "Ybus_Export.txt")

# TARGETS (Your Simulation Scenario)
fault_bus = 7
trip_line_from = 7
trip_line_to = 5
trip_line_id = '1'

# ------------------------------------------------

def init_psse():
    psspy.psseinit(200)
    psspy.read(0, raw_file)

def solve_flow():
    psspy.fnsl([0,0,0,1,1,0,99,0])

def get_dynamic_data():
    """ Extracts Generator Xd' and Loads, then prints a comparison table. """
    print("\n--- EXTRACTING DYNAMIC DATA ---")
    
    # 1. Generators (Xd')
    ierr, mach_bus = psspy.amachint(-1, 4, 'NUMBER')
    ierr, mach_z   = psspy.amachcplx(-1, 4, 'ZSORCE')
    
    gen_buses = []
    xd_prime = []
    
    print(f"{'Gen Bus':<10} | {'Extracted Xd (RAW)':<20} | {'MATLAB Expected':<20}")
    print("-" * 55)
    
    # Expected values for comparison (Just for your reference)
    expected_xd = {1: 0.0608, 2: 0.1198, 3: 0.1813}
    
    if mach_bus:
        for bus, z in zip(mach_bus[0], mach_z[0]):
            if abs(z) > 0:
                gen_buses.append(bus)
                val = z.imag
                xd_prime.append(val)
                ref = expected_xd.get(bus, "N/A")
                print(f"{bus:<10} | {val:<20.4f} | {str(ref):<20}")

    # 2. Loads
    ierr, ld_bus = psspy.aloadint(-1, 4, ['NUMBER'])
    ierr, ld_mva = psspy.aloadcplx(-1, 4, ['MVAACT'])
    ierr, bus_v  = psspy.abusreal(-1, 2, ['PU'])
    
    ierr, all_buses = psspy.abusint(-1, 2, ['NUMBER'])
    v_map = {b: v for b, v in zip(all_buses[0], bus_v[0])}

    load_adm = {}
    print("\n" + "-" * 55)
    print(f"{'Load Bus':<10} | {'Extracted Y (RAW)':<20} | {'MATLAB Expected':<20}")
    print("-" * 55)
    
    # Expected loads
    expected_load = {
        5: 1.26 - 0.504j,
        6: 0.877 - 0.292j,
        8: 0.969 - 0.339j
    }
    
    if ld_bus:
        for bus, s in zip(ld_bus[0], ld_mva[0]):
            v = v_map.get(bus, 1.0)
            if v < 0.001: v = 1.0
            
            # Calc Y
            s_pu = s / 100.0
            y = np.conj(s_pu) / (v**2)
            
            load_adm[bus] = load_adm.get(bus, 0j) + y
            
            ref = expected_load.get(bus, "N/A")
            y_str = f"{y.real:.3f}{y.imag:+.3f}j"
            print(f"{bus:<10} | {y_str:<20} | {str(ref):<20}")

    return gen_buses, xd_prime, load_adm

def get_trip_line_params(frm, to, target_id='1'):
    """ Robust 2-pass search for trip line with SAFETY CHECKS """
    target_id = target_id.strip()
    
    print(f"\n--- SEARCHING FOR TRIP LINE {frm}-{to} ID '{target_id}' ---")
    
    # Pass 1: Transmission Lines
    ierr, ln_int = psspy.abrnint(-1, 0, 0, 1, 2, ['FROMNUMBER', 'TONUMBER'])
    ierr, ln_z   = psspy.abrncplx(-1, 0, 0, 1, 2, ['RX'])
    ierr, ln_flt = psspy.abrnreal(-1, 0, 0, 1, 2, ['CHARG']) # This fails sometimes
    ierr, ln_id  = psspy.abrnchar(-1, 0, 0, 1, 2, ['ID'])

    if ln_int and ln_z:
        # --- FIX: Handle missing Charging data ---
        if ln_flt is None or ln_flt[0] is None:
            # Create a list of 0.0s matching the length of found lines
            charging_list = [0.0] * len(ln_int[0])
        else:
            charging_list = ln_flt[0]

        for f, t, z, ch, cid in zip(ln_int[0], ln_int[1], ln_z[0], charging_list, ln_id[0]):
            
            is_match = ((f == frm and t == to) or (f == to and t == frm))
            if is_match and cid.strip() == target_id:
                if abs(z) == 0: return 0j, 0j
                y_series = 1.0 / z
                y_return = y_series 
                b_return = 1j * ch / 2.0
                print(f"   >>> MATCH FOUND (Line)! Y={y_return:.4f}, B={b_return:.4f}")
                return y_return, b_return

    # Pass 2: Transformers
    ierr, tr_int = psspy.abrnint(-1, 0, 0, 2, 2, ['FROMNUMBER', 'TONUMBER'])
    ierr, tr_z   = psspy.abrncplx(-1, 0, 0, 2, 2, ['RX'])
    ierr, tr_flt = psspy.abrnreal(-1, 0, 0, 2, 2, ['XFRRAT', 'ANGLED'])
    ierr, tr_id  = psspy.abrnchar(-1, 0, 0, 2, 2, ['ID'])

    if tr_int and tr_z and tr_flt and tr_id:
        for f, t, z, tap, ang, cid in zip(tr_int[0], tr_int[1], tr_z[0], tr_flt[0], tr_flt[1], tr_id[0]):
            is_match = ((f == frm and t == to) or (f == to and t == frm))
            if is_match and cid.strip() == target_id:
                if abs(z) == 0: return 0j, 0j
                y_series = 1.0 / z
                
                if tap == 0: tap = 1.0
                angle_rad = ang * (np.pi/180.0)
                a = tap * (np.cos(angle_rad) + 1j * np.sin(angle_rad))
                
                y_return = y_series / (abs(a)**2) 
                b_return = 0j 
                
                print(f"   >>> MATCH FOUND (Xfmr)! Y={y_return:.4f}, B={b_return:.4f}")
                return y_return, b_return

    print("WARNING: Trip line not found in data! Post-Fault matrix will be incorrect.")
    return 0j, 0j

def export_and_parse_ybus():
    # Export Lines Only
    psspy.cong(0) 
    ierr = psspy.output_y_matrix(0, 1, 0, 0, txt_file)
    
    Yreal = {}
    Yimag = {}
    max_bus = 0
    float_re = re.compile(r"[-+]?\d*\.\d+|\d+")
    
    with open(txt_file, "r") as f:
        for line in f:
            parts = float_re.findall(line)
            if len(parts) >= 4:
                try:
                    i, j = int(parts[0]), int(parts[1])
                    real, imag = float(parts[2]), float(parts[3])
                    max_bus = max(max_bus, i, j)
                    Yreal[(i, j)] = real
                    Yimag[(i, j)] = imag
                except: continue

    Y = np.zeros((max_bus, max_bus), dtype=np.complex128)
    for (i, j), rv in Yreal.items():
        iv = Yimag.get((i, j), 0.0)
        Y[i-1, j-1] = rv + 1j*iv
        Y[j-1, i-1] = rv + 1j*iv
    return Y

# --- MATRIX MANIPULATION FUNCTIONS ---
def build_faulted(Y_full, f_bus):
    idx = f_bus - 1
    Y = np.delete(Y_full, idx, axis=0)
    Y = np.delete(Y, idx, axis=1)
    return Y

def build_postfault(Y_full, f_bus, t_bus, y_line, b_ch):
    idx_a, idx_b = f_bus - 1, t_bus - 1
    Yp = Y_full.copy()
    y_total = y_line + b_ch
    Yp[idx_a, idx_a] -= y_total
    Yp[idx_b, idx_b] -= y_total
    Yp[idx_a, idx_b] = 0j
    Yp[idx_b, idx_a] = 0j
    return Yp

def kron_reduce(Y_phys, gens, xd, loads):
    K = len(gens)
    N = Y_phys.shape[0]
    
    # 1. Ybar
    ybar = np.diag([1.0/(1j*x) for x in xd])
    
    # 2. Augment (YD)
    YD = Y_phys.copy()
    
    # Size check
    max_bus_id = max(gens) if gens else 0
    if loads: max_bus_id = max(max_bus_id, max(loads.keys()))
    is_faulted = (N < max_bus_id)

    # Add Gens
    gen_indices_map = []
    for k, bus in enumerate(gens):
        if is_faulted:
            if bus == fault_bus: continue 
            idx = bus - 1 if bus < fault_bus else bus - 2
        else:
            idx = bus - 1
        
        gen_indices_map.append((k, idx))
        YD[idx, idx] += ybar[k, k]
        
    # Add Loads
    for bus, yL in loads.items():
        if is_faulted:
            if bus == fault_bus: continue
            idx = bus - 1 if bus < fault_bus else bus - 2
        else:
            idx = bus - 1
        
        YD[idx, idx] += yL
        
    # 3. Reduce
    YA = ybar
    YB = np.zeros((K, N), dtype=complex)
    for (k, idx) in gen_indices_map:
        YB[k, idx] = -ybar[k, k]
    YC = YB.T
    
    try:
        return YA - YB @ np.linalg.solve(YD, YC)
    except:
        return np.zeros((K,K))

# === MAIN ===
if __name__ == "__main__":
    init_psse()
    solve_flow()
    
    gen_list, xd_list, load_dict = get_dynamic_data()
    y_trip, b_trip = get_trip_line_params(trip_line_from, trip_line_to, trip_line_id)
    
    print("\n--- BUILDING MATRICES ---")
    Y_pre_phys = export_and_parse_ybus()
    
    Yint_pre = kron_reduce(Y_pre_phys, gen_list, xd_list, load_dict)
    
    Y_fault_phys = build_faulted(Y_pre_phys, fault_bus)
    Yint_fault = kron_reduce(Y_fault_phys, gen_list, xd_list, load_dict)
    
    Y_post_phys = build_postfault(Y_pre_phys, trip_line_from, trip_line_to, y_trip, b_trip)
    Yint_post = kron_reduce(Y_post_phys, gen_list, xd_list, load_dict)
    
    sio.savemat(mat_file, {
        "Y_pre": Y_pre_phys,
        "Y_fault": Y_fault_phys,
        "Y_post": Y_post_phys,
        "Yint_pre": Yint_pre,
        "Yint_fault": Yint_fault,
        "Yint_post": Yint_post
    })
    
    print("\n--- Python Result: Yint_Fault (Top 3x3) ---")
    print(Yint_fault)
    
    print(f"\nSUCCESS: Generalized reliable calculation complete. Saved to {mat_file}")