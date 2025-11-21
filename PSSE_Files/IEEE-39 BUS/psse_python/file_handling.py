
import psse3603
import psspy
psspy.psseinit()
savfile='savnw.sav'
psspy.psseinit()
psspy.case(savfile)
err,(genbuses)=psspy.amachint(-1,1,'NUMBER')
err,(genid,)=psspy.amachchar(-1,1,'ID')
generators_zip=zip(genbuses,genid)
generators=list(generators_zip)
print(generators)

