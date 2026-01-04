import psse3603
import psspy
import dyntools
import os
import config

# ======================================================================
# 1. SETUP (centralized config)
# ======================================================================
work_dir = config.WORK_DIR
result_dir = config.RESULT_DIR
raw_file = config.raw_path("IEEE39bus1.raw")
dyr_file = config.DYR_FILE
out_file = config.out_path("IEEE9_CCT_Search.out")

# Initialize PSS/E
_i = psspy.getdefaultint()
_f = psspy.getdefaultreal()
psspy.psseinit(config.PSSE_INIT)

# ==============================================================================
# 2. SIMULATION ENGINE
# ==============================================================================
def run_simulation(clearing_duration):
    """
    Runs the simulation for a specific fault clearing duration.
    Returns: (is_stable, max_spread_found)
    """
    # --- Simulation Parameters ---
    fault_bus = 29
    from_bus = 29
    to_bus = 26
    line_id = '1 '
    
    t_fault_start = 1.0
    t_clear = t_fault_start + clearing_duration
    t_end = 4.0 

    # --- Clean Previous Run Output ---
    if os.path.exists(out_file):
        try:
            os.remove(out_file)
        except:
            pass

    # --- Load Data ---
    psspy.read(0, raw_file)
    psspy.dyre_new([1,1,1,1], dyr_file, "", "", "")

    # --- Solver Settings (use centralized timestep) ---
    psspy.dynamics_solution_param_2([_i]*8, [_f, _f, config.DYNAMICS_STEP, _f, _f, _f, _f, _f])
    psspy.fnsl([0,0,0,1,1,0,99,0])

    # --- Convert Generators ---
    psspy.cong(0) 
    psspy.conl(0, 1, 1, [0, 0], [0.0, 100.0, 0.0, 100.0])
    psspy.conl(0, 1, 2, [0, 0], [0.0, 100.0, 0.0, 100.0])
    psspy.conl(0, 1, 3, [0, 0], [0.0, 100.0, 0.0, 100.0])
    psspy.fact()
    psspy.tysl(0)

    # --- Setup Channels ---
    psspy.delete_all_plot_channels()
    # 1. Angles (Channels 1-10 for IEEE 39)
    psspy.chsb(0,1,[-1,-1,-1,1,1,0])
    # 2. Pmech (Channels 11-20)
    psspy.chsb(0,1,[-1,-1,-1,1,6,0])
    # 3. Pelec (Channels 21-30)
    psspy.chsb(0,1,[-1,-1,-1,1,2,0])
    # 4. Speed (Channels 31-40)
    psspy.chsb(0,1,[-1,-1,-1,1,7,0])

    # --- Run Simulation Sequence ---
    ierr = psspy.strt(0, out_file)
    if ierr != 0:
        return False, 999.0 # Fail

    # Steady State
    psspy.run(0, t_fault_start, 0, 1, 0)

    # Apply Fault
    psspy.dist_bus_fault(fault_bus, 1, 0.0, [0.0, -0.2E+10])
    psspy.run(0, t_clear, 0, 1, 0)

    # Clear Fault
   # psspy.dist_branch_trip(from_bus, to_bus, line_id)
    psspy.dist_clear_fault(1)
    psspy.run(0, t_end, 0, 1, 0)

    # --- Check Stability ---
    try:
        chnf_obj = dyntools.CHNF(out_file)
        _, chanid, chandata = chnf_obj.get_data()
    except:
        return False, 999.0
    
    if 'time' not in chandata or len(chandata['time']) == 0:
        return False, 999.0

    angle_keys = range(1, 11) 
    num_steps = len(chandata['time'])
    
    max_spread_over_run = 0.0
    is_stable = True
    
    for i in range(num_steps):
        angles = []
        for k in angle_keys:
            if k in chandata:
                angles.append(chandata[k][i])
        
        if angles:
            spread = max(angles) - min(angles)
            if spread > max_spread_over_run:
                max_spread_over_run = spread
            
            # Instability Threshold (180 deg)
            if spread > 180.0: 
                is_stable = False
                
    return is_stable, max_spread_over_run

# ==============================================================================
# 3. PHASE 1: BINARY SEARCH FOR CCT
# ==============================================================================
print("\n========================================")
print("   PHASE 1: FINDING CCT (Binary Search)")
print("========================================\n")

# Silence PSS/E output for the loop
psspy.report_output(6, "", []) 
psspy.progress_output(6, "", [])

t_min = 0
t_max = 0.400
precision = 0.00001 

while (t_max - t_min) > precision:
    t_test = (t_max + t_min) / 2
    
    # Run Simulation
    is_stable, spread = run_simulation(t_test)
    
    # Determine Status String for Output
    if is_stable:
        status = "STABLE  "
        t_min = t_test
    else:
        status = "UNSTABLE"
        t_max = t_test
    
    # OUTPUT: Print progress to terminal
    print(f"Testing t_cl = {t_test:.4f} s -> {status} (Max Angle: {spread:.2f} deg)")
        
CCT = t_min 

# Re-enable output
psspy.report_output(2, "", [])
psspy.progress_output(2, "", [])

# ==============================================================================
# 4. PHASE 2: RE-RUN AT CCT & TERMINAL OUTPUT
# ==============================================================================
print("\n========================================")
print(f"   PHASE 2: Verifying CCT @ {CCT:.4f}s")
print("========================================\n")

# Run exactly at CCT (Stable Limit) to get final data
final_stable, final_max_angle = run_simulation(CCT)

# --- TERMINAL OUTPUT SUMMARY ---
print("\n" + "#"*50)
print("             FINAL SIMULATION RESULTS")
print("#"*50)
print(f" CRITICAL CLEARING TIME (CCT) : {CCT:.4f} seconds")
print(f" MAXIMUM ROTOR ANGLE          : {final_max_angle:.2f} degrees")

if final_stable:
    print(" SYSTEM STATUS                : STABLE (Boundary Case)")
else:
    # This might happen if 0.001 precision falls exactly on edge
    print(" SYSTEM STATUS                : UNSTABLE (Marginal)")
print("#"*50 + "\n")


print("\nExporting CCT to MATLAB format...")
# Define output path
mat_file_path = os.path.join(result_dir, "CCT_TimeDomain.mat")

try:
    import scipy.io
    # Save as .mat file (Variable name in MATLAB will be 'CCT_TD')
    scipy.io.savemat(mat_file_path, {'CCT_TD': CCT})
    print(f" [SUCCESS] Saved .mat file: {mat_file_path}")
    print("           Load in MATLAB using: load('CCT_TimeDomain.mat')")

except ImportError:
    # Fallback for PSS/E environments without scipy
    txt_file_path = os.path.join(result_dir, "CCT_TimeDomain.txt")
    with open(txt_file_path, "w") as f:
        f.write(f"{CCT:.6f}")
    print(f" [INFO] 'scipy' not found. Saved as text file: {txt_file_path}")
    print("        Load in MATLAB using: CCT_TD = load('CCT_TimeDomain.txt')")