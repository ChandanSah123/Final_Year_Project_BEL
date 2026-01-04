import os
import sys
import textwrap
import scipy.io
import matplotlib.pyplot as plt
import config

try:
    import psse3603; import psspy; import dyntools
except ImportError:
    sys.path.append(r"C:\Program Files\PTI\PSSE35\PSSPY37")
    sys.path.append(r"C:\Program Files\PTI\PSSE35\PSSBIN")
    import psspy; import dyntools

_i = psspy.getdefaultint(); _f = psspy.getdefaultreal()
psspy.psseinit(1000)

# ==============================================================================
# 1. SETUP & CONFIGURATION (LINE FAULT)
# ==============================================================================
work_dir = config.WORK_DIR
result_dir = config.RESULT_DIR
raw_file = config.RAW_FILE
dyr_file = config.DYR_FILE
out_file = config.OUT_FILE.replace('.out', '_CCT_Line.out')

# ==============================================================================
# 2. SIMULATION ENGINE (LINE FAULT - TRIP)
# ==============================================================================
def run_simulation(clearing_duration):
    t_fault_start = config.T_FAULT_START
    t_clear = t_fault_start + clearing_duration
    t_end = config.T_END
    
    from_bus = config.TRIP_LINE_FROM
    to_bus   = config.TRIP_LINE_TO
    ckt_id   = config.LINE_ID

    if os.path.exists(out_file):
        try: os.remove(out_file)
        except: pass

    psspy.read(0, raw_file)
    psspy.dyre_new([1,1,1,1], dyr_file, "", "", "")
    psspy.dynamics_solution_param_2([_i]*8, [_f, _f, 0.001, _f, _f, _f, _f, _f])
    psspy.fnsl([0,0,0,1,1,0,99,0])

    psspy.cong(0)
    psspy.conl(0, 1, 1, [0, 0], [0.0, 100.0, 0.0, 100.0])
    psspy.conl(0, 1, 2, [0, 0], [0.0, 100.0, 0.0, 100.0])
    psspy.conl(0, 1, 3, [0, 0], [0.0, 100.0, 0.0, 100.0])
    psspy.fact(); psspy.tysl(0)

    psspy.delete_all_plot_channels()
    psspy.chsb(0,1,[-1,-1,-1,1,1,0]) # Angles
    psspy.chsb(0,1,[-1,-1,-1,1,6,0])
    psspy.chsb(0,1,[-1,-1,-1,1,2,0])
    psspy.chsb(0,1,[-1,-1,-1,1,7,0])

    ierr = psspy.strt(0, out_file)
    if ierr != 0: return False, 999.0

    psspy.run(0, t_fault_start, 0, 1, 0)
    
    # Fault at Line Terminal
    psspy.dist_bus_fault(from_bus, 1, 0.0, [0.0, -0.2E+10])
    
    psspy.run(0, t_clear, 0, 1, 0)
    
    # Trip Line & Clear
    psspy.dist_branch_trip(from_bus, to_bus, ckt_id)
    psspy.dist_clear_fault(1)
    
    psspy.run(0, t_end, 0, 1, 0)

    try:
        chnf_obj = dyntools.CHNF(out_file)
        _, chanid, chandata = chnf_obj.get_data()
    except: return False, 999.0

    if 'time' not in chandata: return False, 999.0
    
    angle_keys = [1, 2, 3]
    num_steps = len(chandata['time'])
    max_spread = 0.0
    is_stable = True
    
    for i in range(num_steps):
        angles = [chandata[k][i] for k in angle_keys if k in chandata]
        if angles:
            spread = max(angles) - min(angles)
            if spread > max_spread: max_spread = spread
            if spread > 180.0: is_stable = False
            
    return is_stable, max_spread

# ==============================================================================
# 3. BINARY SEARCH FOR CCT
# ==============================================================================
print("\n=== FINDING CCT: LINE FAULT (TRIP) ===")
psspy.report_output(6, "", []); psspy.progress_output(6, "", [])

t_min, t_max = 0.100, 0.500
precision = 0.0001

while (t_max - t_min) > precision:
    t_test = (t_max + t_min) / 2
    is_stable, spread = run_simulation(t_test)
    status = "STABLE  " if is_stable else "UNSTABLE"
    if is_stable: t_min = t_test
    else: t_max = t_test
    print(f"Testing t_cl = {t_test:.4f} s -> {status} (Spread: {spread:.2f})")

CCT = t_min
psspy.report_output(2, "", []); psspy.progress_output(2, "", [])

# ==============================================================================
# 4. FINAL VERIFICATION & EXPORT
# ==============================================================================
print(f"\nFINAL CCT: {CCT:.4f} seconds")
mat_file_path = os.path.join(result_dir, "CCT_TimeDomain.mat")
scipy.io.savemat(mat_file_path, {'CCT_TD': CCT})
print(f"Saved to {mat_file_path}")