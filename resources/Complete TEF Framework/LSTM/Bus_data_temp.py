import psse3603
import psspy
import dyntools
import pandas as pd
import os
FAULT_TIME = 1.0
SIM_END_TIME = 5

# Paths
work_dir = r"C:\Users\Dell\OneDrive\Desktop\LSTM"
raw_file = os.path.join(work_dir, "ieee39bus1.raw")
dyr_file = os.path.join(work_dir, "ieee39buscls.dyr")
out_file = os.path.join(work_dir, "temp_run.out") # Temporary file for each run

# Buses to simulate faults on (IEEE 9 Bus has buses 1 to 9)
fault_buses = [1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39]

# ==============================================================================
# 2. INITIALIZE PSS/E
# ==============================================================================
_i = psspy.getdefaultint()
_f = psspy.getdefaultreal()
psspy.psseinit(50)

# Initialize the CSV headers we will use for every duration
headers = ['Time', 'FaultLocation', 'FaultDuration', 
           'Ang_G1', 'Ang_G2', 'Ang_G3', 'Ang_G4', 'Ang_G5', 'Ang_G6', 'Ang_G7', 'Ang_G8', 'Ang_G9', 'Ang_G10',
           'Spd_G1', 'Spd_G2', 'Spd_G3', 'Spd_G4', 'Spd_G5', 'Spd_G6', 'Spd_G7', 'Spd_G8', 'Spd_G9', 'Spd_G10',
           'Pe_G1', 'Pe_G2', 'Pe_G3', 'Pe_G4', 'Pe_G5', 'Pe_G6', 'Pe_G7', 'Pe_G8', 'Pe_G9', 'Pe_G10',
           'Vt_G1', 'Vt_G2', 'Vt_G3', 'Vt_G4', 'Vt_G5', 'Vt_G6', 'Vt_G7', 'Vt_G8', 'Vt_G9', 'Vt_G10']

# ==============================================================================
# 3. MAIN LOOP
# ==============================================================================
# Create durations list 0.1, 0.2, ..., 1.0
durations = [round(i * 0.1, 1) for i in range(1, 11)]

for FAULT_DURATION in durations:
    # CSV file per clearing time
    csv_file = os.path.join(work_dir, f"Bus_{FAULT_DURATION:.1f}.csv")

    # Create empty CSV with headers for this duration
    df_init = pd.DataFrame(columns=headers)
    df_init.to_csv(csv_file, index=False)
    print(f"Created CSV for duration {FAULT_DURATION}s: {csv_file}")

    # Loop over each bus for this clearing time
    for f_bus in fault_buses:
        print(f"--- Simulating Fault at Bus {f_bus} (Duration: {FAULT_DURATION}s) ---")

        # --- A. RUN SIMULATION ---
        psspy.read(0, raw_file)
        psspy.dyre_new([1,1,1,1], dyr_file, "", "", "")

        # Solver Settings (0.001s step size as requested)
        psspy.dynamics_solution_param_2([_i]*8, [_f, _f, 0.001, _f, _f, _f, _f, _f])
        psspy.fnsl([0,0,0,1,1,0,99,0])

        # Network Conversion
        psspy.cong(0)
        psspy.conl(0, 1, 1, [0, 0], [0.0, 100.0, 0.0, 100.0])
        psspy.conl(0, 1, 2, [0, 0], [0.0, 100.0, 0.0, 100.0])
        psspy.conl(0, 1, 3, [0, 0], [0.0, 100.0, 0.0, 100.0])

        # Setup Channels (Only recording Gens 1..10 as before)
        psspy.delete_all_plot_channels()
        psspy.chsb(0,1,[-1,-1,-1,1,1,0]) # angle channel
        psspy.chsb(0,1,[-1,-1,-1,1,7,0]) # speed channel
        psspy.chsb(0,1,[-1,-1,-1,1,2,0]) # electrical power channel
        psspy.chsb(0,1,[-1,-1,-1,1,4,0]) # terminal voltage channel

        # Run Dynamics
        psspy.strt_2([0,0], out_file)
        psspy.run(0, FAULT_TIME, 0, 1, 1)

        # Apply Fault
        psspy.dist_bus_fault(f_bus, 1, 0.0, [0.0, -0.2E+10])
        psspy.run(0, FAULT_TIME + FAULT_DURATION, 0, 1, 1)

        # Clear Fault
        psspy.dist_clear_fault(1)
        psspy.run(0, SIM_END_TIME, 0, 1, 1)

        # --- B. EXTRACT DATA USING DYNTOOLS ---
        chnf_obj = dyntools.CHNF(out_file)
        short_title, chanid, chandata = chnf_obj.get_data()

        # Extract channel data (keys assumed as in original script)
        time_data = chandata['time']
        ang_g1 = chandata[1]
        ang_g2 = chandata[2]
        ang_g3 = chandata[3]
        ang_g4 = chandata[4]
        ang_g5 = chandata[5]
        ang_g6 = chandata[6]
        ang_g7 = chandata[7]
        ang_g8 = chandata[8]
        ang_g9 = chandata[9]
        ang_g10 = chandata[10]
        spd_g1 = chandata[11]
        spd_g2 = chandata[12]
        spd_g3 = chandata[13]
        spd_g4 = chandata[14]
        spd_g5 = chandata[15]
        spd_g6 = chandata[16]
        spd_g7 = chandata[17]
        spd_g8 = chandata[18]
        spd_g9 = chandata[19]
        spd_g10 = chandata[20]
        pe_g1 = chandata[21]
        pe_g2 = chandata[22]
        pe_g3 = chandata[23]
        pe_g4 = chandata[24]
        pe_g5 = chandata[25]
        pe_g6 = chandata[26]
        pe_g7 = chandata[27]
        pe_g8 = chandata[28]
        pe_g9 = chandata[29]
        pe_g10 = chandata[30]
        vt_g1 = chandata[31]
        vt_g2 = chandata[32]
        vt_g3 = chandata[33]
        vt_g4 = chandata[34]
        vt_g5 = chandata[35]
        vt_g6 = chandata[36]
        vt_g7 = chandata[37]
        vt_g8 = chandata[38]
        vt_g9 = chandata[39]
        vt_g10 = chandata[40]

        # --- C. PREPARE DATAFRAME ---
        df_temp = pd.DataFrame({
            'Time': time_data,
            'FaultLocation': f_bus,
            'FaultDuration': FAULT_DURATION,
            'Ang_G1': ang_g1,
            'Ang_G2': ang_g2,
            'Ang_G3': ang_g3,
            'Ang_G4': ang_g4,
            'Ang_G5': ang_g5,
            'Ang_G6': ang_g6,
            'Ang_G7': ang_g7,
            'Ang_G8': ang_g8,
            'Ang_G9': ang_g9,
            'Ang_G10': ang_g10,
            'Spd_G1': spd_g1,
            'Spd_G2': spd_g2,
            'Spd_G3': spd_g3,
            'Spd_G4': spd_g4,
            'Spd_G5': spd_g5,
            'Spd_G6': spd_g6,
            'Spd_G7': spd_g7,
            'Spd_G8': spd_g8,
            'Spd_G9': spd_g9,
            'Spd_G10': spd_g10,
            'Pe_G1': pe_g1,
            'Pe_G2': pe_g2,
            'Pe_G3': pe_g3,
            'Pe_G4': pe_g4,
            'Pe_G5': pe_g5,
            'Pe_G6': pe_g6,
            'Pe_G7': pe_g7,
            'Pe_G8': pe_g8,
            'Pe_G9': pe_g9,
            'Pe_G10': pe_g10,
            'Vt_G1': vt_g1,
            'Vt_G2': vt_g2,
            'Vt_G3': vt_g3,
            'Vt_G4': vt_g4,
            'Vt_G5': vt_g5,
            'Vt_G6': vt_g6,
            'Vt_G7': vt_g7,
            'Vt_G8': vt_g8,
            'Vt_G9': vt_g9,
            'Vt_G10': vt_g10
        })

        # --- D. APPEND TO CSV ---
        df_temp.to_csv(csv_file, mode='a', header=False, index=False)

    print(f"Completed duration {FAULT_DURATION}s; data saved to: {csv_file}\n")

print("All durations complete.")