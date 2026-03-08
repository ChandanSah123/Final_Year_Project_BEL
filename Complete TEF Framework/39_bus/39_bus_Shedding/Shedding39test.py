# File:"C:\Users\Acer\Desktop\39_bus_Shedding\Shedding39test.py", generated on WED, JAN 21 2026  10:26, PSS(R)E Xplore release 36.03.01
import psse3603
import psspy
import os
import redirect
import matplotlib.pyplot as plt
import dyntools  # type: ignore
psspy.psseinit(50)
_i = psspy.getdefaultint()
_f = psspy.getdefaultreal()
psspy.psseinit(50)
redirect.psse2py()

P=[250.0, 521.6, 650.0, 632.0, 508.0, 650.0, 762.0, 540.0, 830.0, 1000]
Xd=[0.0250,0.0500,0.0450,0.0350,0.0890,0.0400,0.044,0.0450,0.0450,0.0040]
H=[42, 30.3, 35.8, 28.6, 26, 34.8, 26.4, 24.3, 34.5, 500]

shg=9-1
shed=0.36
remain=1-shed
Premain=P[shg]*remain
Pshed=P[shg]*shed
Xdremain=Xd[shg]*(1+shed/remain)
Xdshed=Xd[shg]*(1+remain/shed)
Hremain=H[shg]*remain
Hshed=H[shg]*shed

machine=38
fault_bus=29
from_bus=29
to_bus=26
line_id="1"
t_fault_start=1.0
t_clear=0.184
trip_delay=0.2
t_end=5

psspy.read(0,r"""C:\Users\Acer\Desktop\Complete TEF Framework\39_bus\39_bus_Shedding\ieee39shed.raw""")
psspy.dyre_new([1,1,1,1],r"""C:\Users\Acer\Desktop\Complete TEF Framework\39_bus\39_bus_Shedding\ieee39shed.dyr""")
psspy.machine_data_5(machine,r"""2""",[_i,_i,_i,_i,_i,_i,_i],[Pshed,_f,_f,_f,_f,_f,_f,_f,Xdshed,_f,_f,_f,_f,_f,_f,_f,_f],["",""])
psspy.machine_chng_5(machine,r"""1""",[_i,_i,_i,_i,_i,_i,_i],[Premain,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f],["",""])
psspy.machine_chng_5(machine,r"""2""",[_i,_i,_i,_i,_i,_i,_i],[Pshed,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f],["",""])
psspy.machine_chng_5(machine,r"""1""",[_i,_i,_i,_i,_i,_i,_i],[_f,_f,_f,_f,_f,_f,_f,_f,Xdremain,_f,_f,_f,_f,_f,_f,_f,_f],["",""])
psspy.machine_chng_5(machine,r"""2""",[_i,_i,_i,_i,_i,_i,_i],[_f,_f,_f,_f,_f,_f,_f,_f,Xdshed,_f,_f,_f,_f,_f,_f,_f,_f],["",""])
psspy.change_plmod_con(machine,r"""1""",r"""GENCLS""",1,Hremain)
psspy.add_plant_model(machine,r"""2""",1,r"""GENCLS""",0,"",0,[],[],2,[0.0,0.0])
psspy.change_plmod_con(machine,r"""2""",r"""GENCLS""",1,Hshed)
psspy.dynamics_solution_param_2([_i,_i,_i,_i,_i,_i,_i,_i],[_f,_f, 0.001,_f,_f,_f,_f,_f])
psspy.fnsl([0,0,0,1,1,0,99,0])
psspy.cong(0)
psspy.conl(0,1,1,[0,0],[0.0, 100.0,0.0, 100.0])
psspy.conl(0,1,2,[0,0],[0.0, 100.0,0.0, 100.0])
psspy.conl(0,1,3,[0,0],[0.0, 100.0,0.0, 100.0])

psspy.chsb(0,1,[-1,-1,-1,1,1,0])
psspy.chsb(0,1,[-1,-1,-1,1,7,0])
psspy.strt_2([0,0],r"""C:\Users\Acer\Desktop\Complete TEF Framework\39_bus\39_bus_Shedding\shedding39.out""")
psspy.run(0,t_fault_start,0,1,1)
psspy.dist_3phase_bus_fault(from_bus,0,1,0,[0.0,-0.2E+10])
#psspy.dist_bus_fault(fault_bus,1,0.0,[0.0,-0.2E+10])
psspy.run(0,t_fault_start+t_clear,0,1,1)
psspy.dist_clear_fault(1)
psspy.run(0,t_fault_start+t_clear+trip_delay,0,1,1)
psspy.dist_machine_trip(machine,r"""2""")
psspy.run(0,t_end,0,1,1)



out_file = r"""C:\Users\Acer\Desktop\Complete TEF Framework\39_bus\39_bus_Shedding\shedding39.out"""

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
    for ch in range(1, 12):
        if ch in chandata:
            clean_label = str(chanid[ch]).strip()
            plt.plot(time_data, chandata[ch], label=f"Ch{ch}: {clean_label}")
    plt.title(f"Rotor Angles )")
    plt.ylabel("Angle (Deg)")
    plt.grid(True)
    plt.legend()
    plt.show()


print("--- Script Finished ---")
