import copy
import numpy as np


l = Line("lp test")

if 0:
	l.add(OBJET2())
	for x in xrange(2):
		l.add(DRIFT("d1", XL=10))
		l.add(QUADRUPO("d1", XL=10, B_0=5 ))
		l.add(CHANGREF("d1", XCE=0, YCE=0, ALE=-360/8))
		l.add(CHANGREF("d1", XCE=0, YCE=1, ALE=0))
		l.add(DRIFT("d1", XL=10))

if 1:
	emma = Line('emma')
	xpas = (2,20,2)
	cells = 42
	angle = 360/cells
	d_offset = 34.048 * mm
	f_offset = 7.514 * mm
	#lengths
	ld = 210 * mm
	sd = 50 * mm
	fq = 58.782 * mm
	dq = 75.699 * mm
	# quad radius
	fr = 370 * mm
	dr = 530 * mm
	fb = -6.695 * fr * T
	db = 4.704 * dr * T

	for x in xrange(42):
		emma.add(DRIFT("start", XL=0* cm_))
		emma.add(FAISCNL("start", FNAME='zgoubi.fai'))
		emma.add(DRIFT('ld', XL=ld*cm_/2))
		emma.add(CHANGREF(ALE=-angle))
		emma.add(CHANGREF(YCE=d_offset*cm_))
		emma.add(QUADRUPO('defoc', XL=dq*cm_, R_0=dr*cm_, B_0=db*kgauss_, XPAS=xpas, KPOS=1))
		emma.add(CHANGREF(YCE=-d_offset*cm_))
		emma.add(DRIFT('sd', XL=sd*cm_))
		emma.add(CHANGREF(YCE=f_offset*cm_))
		emma.add(QUADRUPO('foc', XL=fq*cm_, R_0=fr*cm_, B_0=fb*kgauss_, XPAS=xpas, KPOS=1))
		emma.add(FAISCNL("ffoc", FNAME='zgoubi.fai'))
		emma.add(CHANGREF(YCE=-f_offset*cm_))
		emma.add(DRIFT('ld', XL=ld*cm_/2))
		emma.add(DRIFT("end", XL=0* cm_))
		emma.add(FAISCNL("end", FNAME='zgoubi.fai'))

	l = emma


labels = set()
n = 0
for e in l.elements():
	if hasattr(e, 'label1'):
		lab = e.label1.strip()
		if lab != "":
			if lab in labels:
				lab = lab + str(n)
				n += 1
				e.set(label1=lab)
			
			labels.add(lab)
			


tline = Line('test_line')
ob = OBJET2()
rigidity = - ke_to_rigidity(10e6,ELECTRON_MASS)
ob.set(BORO=rigidity)
tline.add(ob)
tline.add(ELECTRON())
tline.add(l)
tline.add(REBELOTE(NPASS=9, K=99))
tline.add(END())


closed_orbit =  find_closed_orbit_range(tline, init_YTZP=[0,0,0,0], max_iterations=100)
ob.clear()
Y,T,Z,P = closed_orbit
#ob.add(Y=Y, T=T, Z=Z, P=P, D=1)
#for dy in np.linspace(-1,1,10):
#	ob.add(Y=Y+dy, T=T, Z=Z, P=P, D=1)

ob.add(Y=Y+0.1, T=T, Z=Z, P=P, D=1)

print closed_orbit
tline.full_tracking()
res= tline.run()
ftrack = res.get_all('fai')
ptrack = res.get_all('plt')

from zgoubi.lab_plot import LabPlot


lp = LabPlot(l)
#lp.add_tracks(ftrack)
lp.add_tracks(ftrack, ptrack)

#exit()
lp.draw()
#lp.save("emma.pdf")
#lp.save("emma.svg")
lp.show()


