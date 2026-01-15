try:
    import psse3603  # type: ignore
    import psspy     # type: ignore
    import dyntools  # type: ignore
    pss_available = True
except Exception:
    pss_available = False
import matplotlib.pyplot as plt
import os
import config

# ==============================================================================
# 1. SETUP & CONFIGURATION
# ==============================================================================
# UPDATE THESE PATHS TO MATCH YOUR SYSTEM EXACTLY
work_dir = config.WORK_DIR
result_dir = config.RESULT_DIR
raw_file = config.RAW_FILE
dyr_file = config.DYR_FILE
out_file = config.OUT_FILE

# Initialize PSS/E (safe)
try:
    _i = psspy.getdefaultint()
    _f = psspy.getdefaultreal()
    psspy.psseinit(50)
except Exception:
    _i = 0
    _f = 0.0
#==============================================================================
# 2. SIMULATION ENGINE
# ==============================================================================
def run_simulation(clearing_duration):
    """
    Runs the simulation for a specific fault clearing duration.
    Returns: (is_stable, max_spread_found)
    """
    # --- Simulation Parameters ---
    fault_bus = config.FAULT_BUS
    from_bus = config.TRIP_LINE_FROM
    to_bus = config.TRIP_LINE_TO
    line_id = config.LINE_ID

    t_fault_start = config.T_FAULT_START
    t_clear = t_fault_start + clearing_duration
    t_end = config.T_END

    # --- Clean Previous Run Output ---
    if os.path.exists(out_file):
        try:
            os.remove(out_file)
        except:
            pass

    # --- Load Data ---
    psspy.read(0, raw_file)
    psspy.dyre_new([1,1,1,1], dyr_file, "", "", "")

    # --- Solver Settings (CRITICAL: 0.001s Step) ---
    psspy.dynamics_solution_param_2([_i]*8, [_f, _f, 0.001, _f, _f, _f, _f, _f])
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
    
    # 1. Rotor Angles (Channels 1-3 for 3 Generators)
    psspy.chsb(0,1,[-1,-1,-1,1,1,0])
    # 2. Mechanical Power (Channels 4-6)
    psspy.chsb(0,1,[-1,-1,-1,1,6,0])
    # 3. Electrical Power (Channels 7-9)
    psspy.chsb(0,1,[-1,-1,-1,1,2,0])
    # 4. Generator Speed (Channels 10-12)
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
    #psspy.dist_branch_trip(from_bus, to_bus, line_id)
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

    # FIX: IEEE 9 Bus has 3 Generators -> We only check Channels 1, 2, 3
    angle_keys = [1, 2, 3] 
    
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

t_min = 0.100  # Minimum possible clearing time
t_max = 0.500  # Maximum possible clearing time (adjust if needed)
precision = 0.0001 

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
    print(" SYSTEM STATUS                : UNSTABLE (Marginal)")
print("#"*50 + "\n")


# ==============================================================================
# 5. EXPORT CCT TO MATLAB
# ==============================================================================
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