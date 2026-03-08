import psse3603
import psspy
import dyntools
import pandas as pd
import os

# ==============================================================================
# 1. CONFIGURATION & DATA PREPARATION
# ==============================================================================
work_dir = r"C:\Users\Dell\OneDrive\Desktop\LSTM"
raw_file = os.path.join(work_dir, "ieee39bus1.raw")
dyr_file = os.path.join(work_dir, "ieee39buscls.dyr")
out_file = os.path.join(work_dir, "temp_run.out") 

FAULT_TIME = 1.0
SIM_END_TIME = 5.0

# List of branches to simulate faults on: (From_Bus, To_Bus, Circuit_ID)
# Based on the network data provided in your prompt
branch_list = [
    (1, 2, '1'), (1, 2, '2'), (1, 39, '1'), (1, 39, '2'),
    (2, 3, '1'), (2, 3, '2'), (2, 25, '1'), (2, 25, '2'),
    (3, 4, '1'), (3, 4, '2'), (3, 18, '1'), (3, 18, '2'),
    (4, 5, '1'), (4, 5, '2'), (4, 14, '1'), (4, 14, '2'),
    (5, 6, '1'), (5, 6, '2'), (5, 8, '1'), (5, 8, '2'),
    (6, 7, '1'), (6, 7, '2'), (6, 11, '1'), (6, 11, '2'),
    (7, 8, '1'), (7, 8, '2'), 
    (8, 9, '1'), (8, 9, '2'),
    (9, 39, '1'), (9, 39, '2'),
    (10, 11, '1'), (10, 11, '2'), (10, 13, '1'), (10, 13, '2'),
    (13, 14, '1'), (13, 14, '2'),
    (14, 15, '1'), (14, 15, '2'),
    (15, 16, '1'), (15, 16, '2'),
    (16, 17, '1'), (16, 17, '2'), (16, 19, '1'), (16, 19, '2'),
    (16, 21, '1'), (16, 21, '2'), (16, 24, '1'), (16, 24, '2'),
    (17, 18, '1'), (17, 18, '2'), (17, 27, '1'), (17, 27, '2'),
    (21, 22, '1'), (21, 22, '2'),
    (22, 23, '1'), (22, 23, '2'),
    (23, 24, '1'), (23, 24, '2'),
    (25, 26, '1'), (25, 26, '2'),
    (26, 27, '1'), (26, 27, '2'), (26, 28, '1'), (26, 28, '2'), (26, 29, '1'), (26, 29, '2'),
    (28, 29, '1'), (28, 29, '2')
]

headers = ['Time', 'FaultLocation', 'FaultDuration', 
           'Ang_G1', 'Ang_G2', 'Ang_G3', 'Ang_G4', 'Ang_G5', 'Ang_G6', 'Ang_G7', 'Ang_G8', 'Ang_G9', 'Ang_G10',
           'Spd_G1', 'Spd_G2', 'Spd_G3', 'Spd_G4', 'Spd_G5', 'Spd_G6', 'Spd_G7', 'Spd_G8', 'Spd_G9', 'Spd_G10',
           'Pe_G1', 'Pe_G2', 'Pe_G3', 'Pe_G4', 'Pe_G5', 'Pe_G6', 'Pe_G7', 'Pe_G8', 'Pe_G9', 'Pe_G10',
           'Vt_G1', 'Vt_G2', 'Vt_G3', 'Vt_G4', 'Vt_G5', 'Vt_G6', 'Vt_G7', 'Vt_G8', 'Vt_G9', 'Vt_G10']

# ==============================================================================
# 2. INITIALIZE PSS/E
# ==============================================================================
_i = psspy.getdefaultint()
_f = psspy.getdefaultreal()
psspy.psseinit(50)

# ==============================================================================
# 3. MAIN SIMULATION LOOP
# ==============================================================================
# Durations from 0.1 to 1.0
durations = [round(i * 0.1, 1) for i in range(1, 11)]

for FAULT_DURATION in durations:
    # 3.1. Create CSV File for this duration (e.g., Line_0.1.csv)
    # Using format {:.1g} removes trailing zeros if integer, but here we want {:.1f} for consistency (0.1, 1.0)
    filename = f"Line_{FAULT_DURATION:.1f}.csv"
    csv_file = os.path.join(work_dir, filename)
    
    # Initialize Master DataFrame for this duration
    df_init = pd.DataFrame(columns=headers)
    df_init.to_csv(csv_file, index=False)
    print(f"\n[NEW FILE] Created: {filename}")

    # 3.2. Loop through every branch in the list
    for br in branch_list:
        f_from = br[0]
        f_to = br[1]
        f_id = br[2] # Circuit ID
        
        # String identifier for the CSV
        location_str = f"{f_from}_{f_to}_Ckt{f_id.strip()}"
        
        print(f"   -> Simulating Fault: Line {f_from}-{f_to} (ID: {f_id}), Duration: {FAULT_DURATION}s")

        # --- A. SETUP & STEADY STATE ---
        psspy.read(0, raw_file)
        psspy.dyre_new([1,1,1,1], dyr_file, "", "", "")
        
        # Solution Parameters
        psspy.dynamics_solution_param_2([_i]*8, [_f, _f, 0.001, _f, _f, _f, _f, _f])
        psspy.fnsl([0,0,0,1,1,0,99,0])

        # --- B. CHANNEL SETUP ---
        psspy.delete_all_plot_channels()
        psspy.chsb(0,1,[-1,-1,-1,1,1,0]) # Angle
        psspy.chsb(0,1,[-1,-1,-1,1,7,0]) # Speed
        psspy.chsb(0,1,[-1,-1,-1,1,2,0]) # P_elec
        psspy.chsb(0,1,[-1,-1,-1,1,4,0]) # V_term

        # --- C. NETWORK CONVERSION (CRITICAL FIX) ---
        psspy.cong(0) 
        psspy.conl(0, 1, 1, [0, 0], [0.0, 100.0, 0.0, 100.0])
        psspy.conl(0, 1, 2, [0, 0], [0.0, 100.0, 0.0, 100.0])
        psspy.conl(0, 1, 3, [0, 0], [0.0, 100.0, 0.0, 100.0])
        psspy.fact()
        psspy.tysl(0)

        # --- D. DYNAMIC SIMULATION ---
        psspy.strt(0, out_file)
        
        # 1. Run to Fault Time
        psspy.run(0, FAULT_TIME, 0, 1, 0)
        
        # 2. Apply Branch Fault
        # dist_branch_fault(from, to, id, method, impedance_real, [R, X])
        psspy.dist_branch_fault(f_from, f_to, f_id, 1, 0.0, [0.0, -0.2E+10])
        
        # 3. Run for Duration
        psspy.run(0, FAULT_TIME + FAULT_DURATION, 0, 1, 0)
        
        # 4. Trip Branch & Clear Fault
        psspy.dist_branch_trip(f_from, f_to, f_id)
        psspy.dist_clear_fault(1)
        
        # 5. Run to End
        psspy.run(0, SIM_END_TIME, 0, 1, 0)

        # --- E. DATA EXTRACTION ---
        try:
            chnfobj = dyntools.CHNF(out_file)
            short_title, chanid, chandata = chnfobj.get_data()
            
            # Map Data
            time = chandata['time']
            # We assume 10 generators based on your previous scripts
            # PSS/E channels store data in dict keys 1, 2, 3... corresponding to chsb order
            
            data_dict = {
                'Time': time,
                'FaultLocation': [location_str] * len(time),
                'FaultDuration': [FAULT_DURATION] * len(time),
                # Angle (1-10)
                'Ang_G1': chandata[1], 'Ang_G2': chandata[2], 'Ang_G3': chandata[3],
                'Ang_G4': chandata[4], 'Ang_G5': chandata[5], 'Ang_G6': chandata[6],
                'Ang_G7': chandata[7], 'Ang_G8': chandata[8], 'Ang_G9': chandata[9], 'Ang_G10': chandata[10],
                # Speed (11-20)
                'Spd_G1': chandata[11], 'Spd_G2': chandata[12], 'Spd_G3': chandata[13],
                'Spd_G4': chandata[14], 'Spd_G5': chandata[15], 'Spd_G6': chandata[16],
                'Spd_G7': chandata[17], 'Spd_G8': chandata[18], 'Spd_G9': chandata[19], 'Spd_G10': chandata[20],
                # Pe (21-30)
                'Pe_G1': chandata[21], 'Pe_G2': chandata[22], 'Pe_G3': chandata[23],
                'Pe_G4': chandata[24], 'Pe_G5': chandata[25], 'Pe_G6': chandata[26],
                'Pe_G7': chandata[27], 'Pe_G8': chandata[28], 'Pe_G9': chandata[29], 'Pe_G10': chandata[30],
                # Vt (31-40)
                'Vt_G1': chandata[31], 'Vt_G2': chandata[32], 'Vt_G3': chandata[33],
                'Vt_G4': chandata[34], 'Vt_G5': chandata[35], 'Vt_G6': chandata[36],
                'Vt_G7': chandata[37], 'Vt_G8': chandata[38], 'Vt_G9': chandata[39], 'Vt_G10': chandata[40],
            }

            df_temp = pd.DataFrame(data_dict)
            
            # Append to current Duration CSV
            df_temp.to_csv(csv_file, mode='a', header=False, index=False)
            
        except Exception as e:
            print(f"Error extracting data for {location_str}: {e}")

    print(f"Finished processing all lines for Duration {FAULT_DURATION}s")

print("All simulations completed successfully.")