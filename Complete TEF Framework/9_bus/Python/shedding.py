# File:"C:\Users\DELL\Desktop\PhD\energy function\NE\shedding.py", generated on MON, NOV 09 2015   0:27, release 33.05.02
# File:"C:\Users\DELL\Desktop\shedding.py", generated on WED, JUN 24 2015  19:54, release 33.05.02
import psse3605
import psspy
import pssplot
_i = psspy.getdefaultint()
_f = psspy.getdefaultreal()
psspy.psseinit(50)
P=[250.0, 521.6, 650.0, 632.0, 508.0, 650.0, 762.0, 540.0, 830.0, 1000]
Xd=[0.0250,0.0500,0.0450,0.0350,0.0890,0.0400,0.044,0.0450,0.0450,0.0040]
H=[42, 30.3, 35.8, 28.6, 26, 34.8, 26.4, 24.3, 34.5, 500]

shg=7-1
shed=0.24
remain=1-shed
Premain=P[shg]*remain
Pshed=P[shg]*shed
Xdremain=Xd[shg]*(1+shed/remain)
Xdshed=Xd[shg]*(1+remain/shed)
Hremain=H[shg]*remain
Hshed=H[shg]*shed

#b1=15
#b2=3
b=2
m=36
t=1.33
delay=0.1

psspy.read(0,r"""C:\Users\DELL\Desktop\PhD\energy function\NE\ieee39shed.raw""")
psspy.dyre_new([1,1,1,1],r"""C:\Users\DELL\Desktop\PhD\energy function\NE\ieee39shed.dyr""","","","")
psspy.resq(r"""C:\Users\DELL\Desktop\PhD\energy function\NE\ieee39bus.seq""")
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
psspy.strt(0,r"""C:\Users\DELL\Desktop\NE.out""")
psspy.run(0, 1.0,0,1,0)
psspy.dist_bus_fault(b,1,0.0,[0.0,-0.2E+10])
#psspy.dist_spcb_fault_2(b1,b2,r"""2""",[3,0,3,1,0,0,1],[ 0.5,0.0,0.0,0.0,0.0])

psspy.run(0, t,0,1,0)
psspy.dist_clear_fault(1)

#psspy.dist_branch_trip(b1,b2,r"""2""")

psspy.run(0, t+delay,0,1,0)
psspy.dist_machine_trip(m,r"""2""")

psspy.run(0, 5.0,0,1,0)
pssplot.openchandatafile(r"""C:\Users\DELL\Desktop\NE.out""")
pssplot.dragdropplotdata(r"""NE""",r"""1 - ANGL    30[GEN 10      11.000]1""")
pssplot.dragdropplotdata(r"""NE""",r"""2 - ANGL    31[GEN 2       11.000]1""")
pssplot.dragdropplotdata(r"""NE""",r"""3 - ANGL    32[GEN 3       11.000]1""")
pssplot.dragdropplotdata(r"""NE""",r"""4 - ANGL    33[GEN 4       11.000]1""")
pssplot.dragdropplotdata(r"""NE""",r"""5 - ANGL    34[GEN 5       11.000]1""")
pssplot.dragdropplotdata(r"""NE""",r"""6 - ANGL    35[GEN 6       11.000]1""")
pssplot.dragdropplotdata(r"""NE""",r"""7 - ANGL    36[GEN 7       11.000]1""")
pssplot.dragdropplotdata(r"""NE""",r"""8 - ANGL    36[GEN 7       11.000]2""")
pssplot.dragdropplotdata(r"""NE""",r"""9 - ANGL    37[GEN 8       11.000]1""")
pssplot.dragdropplotdata(r"""NE""",r"""10 - ANGL    38[GEN 9       11.000]1""")
pssplot.dragdropplotdata(r"""NE""",r"""11 - ANGL    39[GEN 1       230.00]1""")
