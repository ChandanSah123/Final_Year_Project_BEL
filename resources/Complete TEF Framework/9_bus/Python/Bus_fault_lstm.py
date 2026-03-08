import psse3603
import psspy
import pssplot
_i = psspy.getdefaultint()
_f = psspy.getdefaultreal()
psspy.psseinit(50)
fault_bus=7
fault_time=1
fault_clear=0.2
end=5
raw_file=r"""C:\Users\Acer\Desktop\TEF_Pre\9_bus\Python\IEEE9bus.raw"""
dyr_file=r"""C:\Users\Acer\Desktop\TEF_Pre\9_bus\Python\ieee9bus.dyr"""
out_file=r"""C:\Users\Acer\Desktop\TEF_Pre\9_bus\Python\data_gen.out"""
psspy.read(0,r"""C:\Users\Acer\Desktop\TEF_Pre\9_bus\Python\IEEE9bus.raw""")
psspy.dyre_new([1,1,1,1],r"""C:\Users\Acer\Desktop\TEF_Pre\9_bus\Python\ieee9bus.dyr""")
psspy.dynamics_solution_param_2([_i]*8, [_f, _f, 0.001, _f, _f, _f, _f, _f])
psspy.fnsl([0,0,0,1,1,0,99,0])
psspy.cong(2)
psspy.conl(0,1,1,[0,0],[100.0,0.0,0.0,100.0])
psspy.conl(0,1,2,[0,0],[100.0,0.0,0.0,100.0])
psspy.conl(0,1,3,[0,0],[100.0,0.0,0.0,100.0])
psspy.chsb(0,1,[-1,-1,-1,1,1,0])
psspy.chsb(0,1,[-1,-1,-1,1,7,0])
psspy.strt_2([0,0],r"""C:\Users\Acer\Desktop\TEF_Pre\9_bus\Python\data_gen.out""")
psspy.run(0,fault_time,0,1,1)
psspy.dist_bus_fault(fault_bus, 1, 0.0, [0.0, -0.2E+10])
psspy.run(0,fault_time+fault_clear,0,1,1)
psspy.dist_clear_fault(1)
psspy.run(0,end,0,1,1)
pssplot.openchandatafile(r"""C:\Users\Acer\Desktop\TEF_Pre\9_bus\Python\data_gen.out""")
pssplot.channelfileexcelexport(r"""C:\Users\Acer\Desktop\TEF_Pre\9_bus\Python\data_gen.out""")
