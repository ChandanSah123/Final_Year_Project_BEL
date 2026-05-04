try:
    import psse3605  # type: ignore
    import psspy     # type: ignore
    import dyntools  # type: ignore
    pss_available = True
except Exception:
    pss_available = False
import os
import config
import scipy.io

# ==============================================================================
# 1. SETUP & CONFIGURATION
# ==============================================================================
# --- CRITERIA ---
# 360 degrees allows for wide swings (first swing stability) without 
# declaring instability prematurely.
MAX_ANGLE_THRESHOLD = 360

work_dir = config.WORK_DIR
result_dir = config.RESULT_DIR
raw_file = config.RAW_FILE
dyr_file = config.DYR_FILE
out_file = config.OUT_FILE

# Initialize PSS/E
print("Before Init")
try:
    _i = psspy.getdefaultint()
    _f = psspy.getdefaultreal()
    psspy.psseinit(50)
except Exception:
    _i = 0
    _f = 0.0

print("After Init")

print("--- Starting Robust CCT Search (Step-Forward & Refine) ---")

# ==============================================================================
# 2. SIMULATION ENGINE
# ==============================================================================
def run_simulation(clearing_duration):
    """
    Runs PSS/E simulation. Returns (is_stable, max_spread_found).
    """
    # --- Parameters ---
    fault_bus = config.FAULT_BUS
    from_bus = config.TRIP_LINE_FROM
    to_bus = config.TRIP_LINE_TO
    line_id = config.LINE_ID

    t_fault_start = config.T_FAULT_START
    t_clear = t_fault_start + clearing_duration
    # Simulate for 3 seconds AFTER fault to catch second-swing instability
    t_end = 5

    # --- Clean Previous Output ---
    if os.path.exists(out_file):
        try: os.remove(out_file)
        except: pass

    # --- Load Data ---
    psspy.read(0, raw_file)
    psspy.dyre_new([1,1,1,1], dyr_file, "", "", "")

    # --- Solver Settings ---
    # 0.002s (2ms) is a good balance for CCT search speed vs accuracy
    psspy.dynamics_solution_param_2([_i]*8, [_f, _f, 0.001, _f, _f, _f, _f, _f])
    psspy.fnsl([0,0,0,1,1,0,99,0])

    # --- Convert Generators ---
    psspy.cong(0)
    psspy.conl(0,1,1,[0,0],[0.0, 100.0,0.0, 100.0])
    psspy.conl(0,1,2,[0,0],[0.0, 100.0,0.0, 100.0])
    psspy.conl(0,1,3,[0,0],[0.0, 100.0,0.0, 100.0])
    psspy.fact()
    psspy.tysl(0)

    # --- Setup Channels ---
    psspy.delete_all_plot_channels()
    # Record Angles (Code 1) - Channels 1, 2, 3...
    psspy.chsb(0,1,[-1,-1,-1,1,1,0]) 

    # --- Run Sequence ---
    ierr = psspy.strt(0, out_file)
    if ierr != 0: 
        return False, 999.0 

    # Steady State
    psspy.run(0, t_fault_start, 0, 1, 0)
    
    # Fault
    psspy.dist_bus_fault(fault_bus, 1, 0.0, [0.0, -0.2E+10])
    #psspy.dist_branch_fault(from_bus, to_bus, line_id, 1, 0.0, [0.0, -0.2E+10])
    psspy.run(0, t_clear, 0, 1, 0)

    # Clear
    #psspy.dist_branch_trip(from_bus, to_bus, line_id)
    psspy.dist_clear_fault(1)
    psspy.run(0, t_end, 0, 1, 0)

    # --- Check Stability ---
    try:
        chnf_obj = dyntools.CHNF(out_file)
        _, _, chandata = chnf_obj.get_data()
    except:
        return False, 999.0
    
    if 'time' not in chandata: return False, 999.0

    # Get all integer keys (Angle channels)
    angle_keys = [k for k in chandata.keys() if isinstance(k, int)]
    
    num_steps = len(chandata['time'])
    max_spread_global = 0.0
    is_stable = True
    
    # Check every 5th step (sufficient resolution)
    for i in range(0, num_steps, 5): 
        current_angles = [chandata[k][i] for k in angle_keys]
        if not current_angles: continue
            
        spread = max(current_angles) - min(current_angles)
        if spread > max_spread_global:
            max_spread_global = spread
        
        # Robust Threshold Check
        if spread > MAX_ANGLE_THRESHOLD:
            is_stable = False
            break 
            
    return is_stable, max_spread_global

# ==============================================================================
# 3. ROBUST SEARCH LOGIC
# ==============================================================================
print("\n========================================")
print("   STARTING LINEAR STEP & REFINE SEARCH")
print("========================================\n")

# Silence PSS/E output during loop
psspy.report_output(6, "", []) 
psspy.progress_output(6, "", [])

# --- Search Parameters ---
current_t = 0.050      # Start Time (50ms)
step_size = 0.010      # Initial Step (50ms)
precision_limit = 0.001 # Target Precision (1ms)

last_stable_t = 0.0
found_cct = False

print(f"{'TEST TIME':<12} | {'STATUS':<10} | {'MAX ANGLE':<10}")
print("-" * 40)

while not found_cct:
    # Safety: Stop if CCT > 1.0s (Very unlikely for transient stability)
    if current_t > 1.0: 
        print("Safety limit reached (>1.0s). Simulation stopping.")
        last_stable_t = 1.0
        break

    stable, spread = run_simulation(current_t)
    
    if stable:
        # --- CASE STABLE ---
        print(f"{current_t:.4f} s     | STABLE     | {spread:.2f} deg")
        last_stable_t = current_t
        current_t += step_size
    else:
        # --- CASE UNSTABLE ---
        print(f"{current_t:.4f} s     | UNSTABLE   | {spread:.2f} deg")
        
        # Check if we are precise enough
        if step_size > precision_limit:
            # Not precise enough yet. Step back and refine.
            print(f"   >>> UNSTABLE. Backing up to {last_stable_t:.4f}s and refining step...")
            
            # Reduce step size (e.g., 0.05 -> 0.01 -> 0.001)
            step_size = step_size / 5.0 
            
            # Snap to nice numbers for readability
            if 0.009 < step_size < 0.011: step_size = 0.010
            if 0.0019 < step_size < 0.0021: step_size = 0.002
            if step_size < 0.001: step_size = 0.001
            
            # Restart search from last known stable point
            current_t = last_stable_t + step_size
        else:
            # We are at 1ms precision. This is the CCT.
            found_cct = True

CCT = last_stable_t

# Re-enable output
#psspy.report_output(2, "", [])
#psspy.progress_output(2, "", [])

# ==============================================================================
# 4. EXPORT RESULTS
# ==============================================================================
print("\n" + "#"*50)
print(f" ROBUST CCT FOUND: {CCT:.4f} seconds")
print("#"*50)

# Export to MATLAB
mat_file_path = os.path.join(result_dir, "CCT_TimeDomain.mat")
try:
    scipy.io.savemat(mat_file_path, {'CCT_TD': CCT})
    print(f" [SUCCESS] Saved .mat file: {mat_file_path}")
    print("           Load in MATLAB using: load('CCT_TimeDomain.mat')")
except ImportError:
    print(" [ERROR] Scipy not found. Cannot save .mat file.")