# File:"C:\Users\DELL\Desktop\PhD\energy function\NE\shedding.py", generated on MON, NOV 09 2015   0:27, release 33.05.02
# File:"C:\Users\DELL\Desktop\shedding.py", generated on WED, JUN 24 2015  19:54, release 33.05.02
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
shed=0.335
remain=1-shed
Premain=P[shg]*remain
Pshed=P[shg]*shed
Xdremain=Xd[shg]*(1+shed/remain)
Xdshed=Xd[shg]*(1+remain/shed)
Hremain=H[shg]*remain
Hshed=H[shg]*shed

#b1=15
#b2=3
b=29
m=38
t=1.184
delay=0.2

psspy.read(0,r"""C:\Users\Acer\Desktop\Complete TEF Framework\39_bus\39_bus_Shedding\ieee39shed.raw""")
psspy.dyre_new([1,1,1,1],r"""C:\Users\Acer\Desktop\Complete TEF Framework\39_bus\39_bus_Shedding\ieee39shed.dyr""")
#psspy.resq(r"""C:\Users\Acer\Desktop\Complete TEF Framework\39_bus\39_bus_Shedding\ieee39bus.seq""")
psspy.machine_data_2(m,r"""2""",[_i,_i,_i,_i,_i,_i],[ Pshed, 1.0, 9999,_f,_f,_f,_f,_f, Xdshed,_f,_f,_f,_f,_f,_f,_f,_f])
psspy.machine_chng_2(m,r"""1""",[_i,_i,_i,_i,_i,_i],[ Premain,_f,_f,_f,_f,_f,_f,_f, Xdremain,_f,_f,_f,_f,_f,_f,_f,_f])
psspy.machine_chng_2(m,r"""2""",[_i,_i,_i,_i,_i,_i],[_f,_f, 1.0,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f])
psspy.change_plmod_con(m,r"""1""",r"""GENCLS""",1, Hremain)
psspy.add_plant_model(m,r"""2""",1,r"""GENCLS""",0,"",0,[],[],2,[0.0,0.0])
psspy.change_plmod_con(m,r"""2""",r"""GENCLS""",1, Hshed)
psspy.dynamics_solution_param_2([_i,_i,_i,_i,_i,_i,_i,_i],[_f,_f, 0.001,_f,_f,_f,_f,_f])
psspy.fnsl([0,0,0,1,1,0,99,0])
psspy.nsol([0,0,0,1,1,0,99])
psspy.nsol([0,0,0,1,1,0,99])
psspy.nsol([0,0,0,1,1,0,99])
psspy.chsb(0,1,[-1,-1,-1,1,1,0])
psspy.cong(0)
psspy.conl(0,1,1,[0,0],[0.0, 100.0,0.0, 100.0])
psspy.conl(0,1,2,[0,0],[0.0, 100.0,0.0, 100.0])
psspy.conl(0,1,3,[0,0],[0.0, 100.0,0.0, 100.0])
#psspy.strt(0,r"""C:\Users\DELL\Desktop\NE.out""")
psspy.strt(0,r"""C:\Users\Acer\Desktop\Complete TEF Framework\39_bus\39_bus_Shedding\shedding39.out""")
psspy.run(0, 1.0,0,1,0)
psspy.dist_bus_fault(b,1,0.0,[0.0,-0.2E+10])
#psspy.dist_spcb_fault_2(b1,b2,r"""2""",[3,0,3,1,0,0,1],[ 0.5,0.0,0.0,0.0,0.0])

psspy.run(0, t,0,1,0)
psspy.dist_clear_fault(1)

#psspy.dist_branch_trip(b1,b2,r"""2""")

psspy.run(0, t+delay,0,1,0)
psspy.dist_machine_trip(m,r"""2""")

psspy.run(0, 5.0,0,1,0)



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


