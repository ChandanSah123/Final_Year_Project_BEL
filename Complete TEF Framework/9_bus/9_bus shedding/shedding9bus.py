# File:"C:\Users\Acer\Desktop\9_bus shedding\shedding9bus.py", generated on MON, JAN 12 2026  19:51, PSS(R)E Xplore release 36.03.01

import psse3605  # type: ignore
import psspy     # type: ignore
import dyntools 
import os
import config
import redirect
import matplotlib.pyplot as plt
import dyntools  # type: ignore
psspy.psseinit(50)
_i = psspy.getdefaultint()
_f = psspy.getdefaultreal()
psspy.psseinit(50)
redirect.psse2py()

P=[71.641, 163.0, 85.0]
Xd=[0.0608, 0.1198, 0.1813]
H=[23.64, 6.4, 3.01]
shg=2-1
shed=0.36
remain=1-shed
Premain=P[shg]*remain
Pshed=P[shg]*shed
Xdremain=Xd[shg]*(1+shed/remain)
Xdshed=Xd[shg]*(1+remain/shed)
Hremain=H[shg]*remain
Hshed=H[shg]*shed

machine=2
fault_bus=7
from_bus=7
to_bus=5
line_id="1"
t_fault_start=1.0
t_clear=0.245
trip_delay=0.1
t_end=5

psspy.read(0,r"""C:\Users\Acer\Desktop\Complete TEF Framework\9_bus\9_bus shedding\IEEE9busshed.raw""")
psspy.dyre_new_2([1,1,1,1],r"""C:\Users\Acer\Desktop\Complete TEF Framework\9_bus\9_bus shedding\ieee9busshed.dyr""")
psspy.machine_data_5(machine,r"""2""",[_i,_i,_i,_i,_i,_i,_i],[P[shg],6.654,300.0,-300.0,100.0,10.0,_f,_f,Xd[shg],_f,_f,_f,_f,_f,_f,_f,_f],["",""])
psspy.machine_chng_5(machine,r"""1""",[_i,_i,_i,_i,_i,_i,_i],[Premain,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f],["",""])
psspy.machine_chng_5(machine,r"""2""",[_i,_i,_i,_i,_i,_i,_i],[Pshed,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f],["",""])
psspy.machine_chng_5(machine,r"""1""",[_i,_i,_i,_i,_i,_i,_i],[_f,_f,_f,_f,_f,_f,_f,_f,Xdremain,_f,_f,_f,_f,_f,_f,_f,_f],["",""])
psspy.machine_chng_5(machine,r"""2""",[_i,_i,_i,_i,_i,_i,_i],[_f,_f,_f,_f,_f,_f,_f,_f,Xdshed,_f,_f,_f,_f,_f,_f,_f,_f],["",""])
psspy.add_plant_model(machine,r"""2""",1,r"""GENCLS""",0,"",0,[],[],2,[0.0,0.0])
psspy.change_plmod_con(machine,r"""1""",r"""GENCLS""",1,Hremain)
psspy.change_plmod_con(machine,r"""2""",r"""GENCLS""",1,Hshed)
psspy.dynamics_solution_param_2([_i,_i,_i,_i,_i,_i,_i,_i],[_f,_f, 0.001,_f,_f,_f,_f,_f])
psspy.fnsl([0,0,0,1,1,0,99,0])
psspy.cong(0)
psspy.conl(0,1,1,[0,0],[0.0, 100.0,0.0, 100.0])
psspy.conl(0,1,2,[0,0],[0.0, 100.0,0.0, 100.0])
psspy.conl(0,1,3,[0,0],[0.0, 100.0,0.0, 100.0])
psspy.chsb(0,1,[-1,-1,-1,1,1,0])
psspy.chsb(0,1,[-1,-1,-1,1,2,0])
psspy.chsb(0,1,[-1,-1,-1,1,6,0])
psspy.chsb(0,1,[-1,-1,-1,1,7,0])
psspy.strt_2([0,0],r"""C:\Users\Acer\Desktop\Complete TEF Framework\9_bus\9_bus shedding\9bussheding.out""")
psspy.run(0,t_fault_start,0,1,1)
psspy.dist_3phase_bus_fault(fault_bus,0,1,230.0,[0.0,-0.2E+10])
psspy.run(0,t_fault_start + t_clear,0,1,1)
psspy.dist_clear_fault(1)
#psspy.dist_branch_trip(from_bus,to_bus,line_id)
psspy.run(0,t_fault_start + t_clear + trip_delay,0,1,1)
psspy.dist_machine_trip(2,r"""2""")
psspy.run(0,t_end,0,1,1)

# ==============================================================================
# 1. SETUP & CONFIGURATION
# ==============================================================================
out_file = config.OUT_FILE

# Initialize PSS/E (safe)
try:
    _i = psspy.getdefaultint()
    _f = psspy.getdefaultreal()
    psspy.psseinit(50)
except Exception:
    _i = 0
    _f = 0.0


# ==============================================================================
# 3. PLOTTING
# ==============================================================================

if not os.path.exists(out_file):
    print(f"CRITICAL ERROR: Output file not found at {out_file}")
else:
    chnf_obj = dyntools.CHNF(out_file)
    short_title, chanid, chandata = chnf_obj.get_data()
    time_data = chandata['time']

    print("Generating Plots...")

    # Plot 1: Rotor Angles (Channels 1, 2, 3)
    plt.figure(figsize=(10, 6))
    for ch in range(1, 5):
        if ch in chandata:
            clean_label = str(chanid[ch]).strip()
            plt.plot(time_data, chandata[ch], label=f"Ch{ch}: {clean_label}")
    plt.title(f"Rotor Angles )")
    plt.ylabel("Angle (Deg)")
    plt.grid(True)
    plt.legend()
    plt.show()


print("--- Script Finished ---")