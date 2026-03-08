import psse3603  # type: ignore
import psspy  # type: ignore
import dyntools  # type: ignore
import pandas as pd
import numpy as np
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
            
            # Convert all channel data to numpy arrays and filter after fault clearing
            time = np.array(chandata['time'])
            fault_clear_time = FAULT_TIME + FAULT_DURATION
            mask = time >= fault_clear_time
            
            # Map Data with filtering
            time = time[mask]
            data_dict = {
                'Time': time,
                'FaultLocation': [location_str] * len(time),
                'FaultDuration': [FAULT_DURATION] * len(time),
                # Angle (1-10)
                'Ang_G1': np.array(chandata[1])[mask], 'Ang_G2': np.array(chandata[2])[mask], 'Ang_G3': np.array(chandata[3])[mask],
                'Ang_G4': np.array(chandata[4])[mask], 'Ang_G5': np.array(chandata[5])[mask], 'Ang_G6': np.array(chandata[6])[mask],
                'Ang_G7': np.array(chandata[7])[mask], 'Ang_G8': np.array(chandata[8])[mask], 'Ang_G9': np.array(chandata[9])[mask], 'Ang_G10': np.array(chandata[10])[mask],
                # Speed (11-20)
                'Spd_G1': np.array(chandata[11])[mask], 'Spd_G2': np.array(chandata[12])[mask], 'Spd_G3': np.array(chandata[13])[mask],
                'Spd_G4': np.array(chandata[14])[mask], 'Spd_G5': np.array(chandata[15])[mask], 'Spd_G6': np.array(chandata[16])[mask],
                'Spd_G7': np.array(chandata[17])[mask], 'Spd_G8': np.array(chandata[18])[mask], 'Spd_G9': np.array(chandata[19])[mask], 'Spd_G10': np.array(chandata[20])[mask],
                # Pe (21-30)
                'Pe_G1': np.array(chandata[21])[mask], 'Pe_G2': np.array(chandata[22])[mask], 'Pe_G3': np.array(chandata[23])[mask],
                'Pe_G4': np.array(chandata[24])[mask], 'Pe_G5': np.array(chandata[25])[mask], 'Pe_G6': np.array(chandata[26])[mask],
                'Pe_G7': np.array(chandata[27])[mask], 'Pe_G8': np.array(chandata[28])[mask], 'Pe_G9': np.array(chandata[29])[mask], 'Pe_G10': np.array(chandata[30])[mask],
                # Vt (31-40)
                'Vt_G1': np.array(chandata[31])[mask], 'Vt_G2': np.array(chandata[32])[mask], 'Vt_G3': np.array(chandata[33])[mask],
                'Vt_G4': np.array(chandata[34])[mask], 'Vt_G5': np.array(chandata[35])[mask], 'Vt_G6': np.array(chandata[36])[mask],
                'Vt_G7': np.array(chandata[37])[mask], 'Vt_G8': np.array(chandata[38])[mask], 'Vt_G9': np.array(chandata[39])[mask], 'Vt_G10': np.array(chandata[40])[mask],
            }

            df_temp = pd.DataFrame(data_dict)
            
            # Append to current Duration CSV
            df_temp.to_csv(csv_file, mode='a', header=False, index=False)
            
        except Exception as e:
            print(f"Error extracting data for {location_str}: {e}")

    print(f"Finished processing all lines for Duration {FAULT_DURATION}s")

print("All simulations completed successfully.")