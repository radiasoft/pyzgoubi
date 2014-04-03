import copy
import numpy as np


l = Line("lp test")

if 0:
	# changref
	#l.add(OBJET2())
	for x in xrange(2):
		l.add(DRIFT("d1", XL=10))
		l.add(QUADRUPO("d1", XL=10, B_0=5, KPOS=1))
		l.add(CHANGREF("d1", XCE=0, YCE=0, ALE=-360/8))
		l.add(CHANGREF("d1", XCE=0, YCE=1, ALE=0))
		l.add(DRIFT("d1", XL=10))
if 0:
	# BENDS
	for x in xrange(1):
		l.add(DRIFT("start", XL=0))
		l.add(FAISCNL("start", FNAME='zgoubi.fai'))
		l.add(DRIFT("d1", XL=5))
		l.add(BEND("b1", XPAS=(10,10,10), XL=5, B1=-4, KPOS=3, W_E = radians(0), W_S = radians(0)))
		l.add(DRIFT("d2", XL=5))
		l.add(DRIFT("d2", XL=5))
		l.add(BEND("b2", XPAS=(10,10,10), XL=5, B1=4, KPOS=1, W_E = radians(0), W_S = radians(-20)))
		l.add(DRIFT("d3", XL=5))
		l.add(DRIFT("d4", XL=5))
		l.add(BEND("b3", XPAS=(10,10,10), XL=5, B1=-4, KPOS=1, W_E = radians(0), W_S = radians(0)))
		l.add(DRIFT("d5", XL=5))
		l.add(DRIFT("end",XL=0))
		l.add(FAISCNL("end", FNAME='zgoubi.fai'))
if 0:
	# DIPOLE
	ref_rid = -ke_to_rigidity(10e6, ELECTRON_MASS)
	dipole_angle = 2*pi/36
	dipole_field = -0.1*T
	dipole_bend_radius  = abs(ref_rid / (dipole_field *kgauss_ ) *cm)
	ad = abs(degrees(dipole_angle))
	for x in xrange(1):
		l.add(DRIFT("start", XL=0))
		l.add(FAISCNL("start", FNAME='zgoubi.fai'))
		l.add(DRIFT("d1", XL=5))
		#l.add(BEND("b1", XPAS=(10,10,10), XL=5, B1=-4, KPOS=3, W_E = radians(0), W_S = radians(0)))
		l.add(DIPOLE("b1", AT=ad, RM=dipole_bend_radius*cm_, ACN=ad/2, B_0=dipole_field*kgauss_, N=2,
			OMEGA_E=ad/2, OMEGA_S=-ad/2,
			R1_E=1e9, R2_E=1e9, R1_S=1e9, R2_S=1e9, R1_L=1e9, R2_L=1e9,R3=1e9,
			U1_E=-1e9, U2_E=1e9, U1_S=-1e9, U2_S=1e9, U1_L=-1e9, U2_L=1e9,
			IORDRE=2, Resol=10,
			LAM_E=1,NCE=2,CE_0=0.1,CE_1=2, XI_E=1,SHIFT_E=1,
			LAM_S=0.1,
			KPOS=2, RE=dipole_bend_radius*cm_, RS=dipole_bend_radius*cm_,
			XPAS=1
			))


		l.add(DRIFT("d2", XL=5))
		l.add(DRIFT("d3", XL=5))
		l.add(DIPOLE("b2", AT=ad, RM=dipole_bend_radius*cm_, ACN=ad/2, B_0=-dipole_field*kgauss_, N=-2,
			OMEGA_E=ad/2, OMEGA_S=-ad/2,
			R1_E=1e9, R2_E=1e9, R1_S=1e9, R2_S=1e9, R1_L=1e9, R2_L=1e9,R3=1e9,
			U1_E=-1e9, U2_E=1e9, U1_S=-1e9, U2_S=1e9, U1_L=-1e9, U2_L=1e9,
			IORDRE=2, Resol=10,
			LAM_E=1,NCE=2,CE_0=0.1,CE_1=2, XI_E=1,SHIFT_E=1,
			LAM_S=0.1,
			KPOS=2, RE=dipole_bend_radius*cm_, RS=dipole_bend_radius*cm_,
			XPAS=1
			))

		l.add(DRIFT("d4", XL=5))
		l.add(FAISCNL("end", FNAME='zgoubi.fai'))

if 1:
	# EMMA, quads and change refs
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

	l = emma * 2


from collections import Counter
import copy

def uniquify_labels(line):
	"make names unique by appending digits"
	max_label = 10
	label1s = Counter()
	dup_names = set()

	# dont need to rename if already unique
	for e in line.elements():
		if hasattr(e, 'label1'):
			lab1 = e.label1.strip()
			if lab1 != "":
				label1s[lab1] += 1
				if label1s[lab1] > 1:
					dup_names.add(lab1)
	
	label1s = Counter() # reset counter
	# new line with copy of elements
	# avoids issues if a line contains the same element instance multiple times
	new_line = Line(line.name)
	for e in line.elements():
		new_e = copy.deepcopy(e)
		if hasattr(new_e, 'label1'):
			lab1 = new_e.label1.strip()
			if lab1 != "" and lab1 in dup_names:
				new_lab1 = lab1 + str(label1s[lab1])
				if len(new_lab1) > max_label:
					raise ValueError("Renaming '%s' to '%s' exceeds label length %s"%(lab1, new_lab1, max_label))
				label1s[lab1] += 1
				new_e.set(label1 = new_lab1)
		new_line.add(new_e)
	return new_line



l = uniquify_labels(l)


tline = Line('test_line')
ob = OBJET2()
rigidity = - ke_to_rigidity(10e6,ELECTRON_MASS)
ob.set(BORO=rigidity)
ob.add(D=1)
tline.add(ob)
tline.add(ELECTRON())
tline.add(l)
tline.add(REBELOTE(NPASS=0, K=99))
tline.add(END())



closed_orbit = [0,0,0,0]
#closed_orbit =  find_closed_orbit_range(tline, init_YTZP=[0,0,0,0], max_iterations=100)
ob.clear()
Y,T,Z,P = closed_orbit
#ob.add(Y=Y, T=T, Z=Z, P=P, D=1)
if 0:
	for dy in np.linspace(-4,4,10):
		ob.add(Y=Y+dy, T=T, Z=Z, P=P, D=1)

if 1:
	for dd in np.linspace(0.1,1.5,10):
		ob.add(Y=Y, T=T, Z=Z, P=P, D=dd)

# for magnet field scan
if 0:
	for dy in np.linspace(-30,30,100):
		ob.add(Y=Y+dy, T=T, Z=Z, P=P, D=1)

#ob.add(Y=Y+1, T=T, Z=Z, P=P, D=1)
#ob.add(Y=Y+1, T=T, Z=Z, P=P, D=2)
#ob.add(Y=Y+1, T=T, Z=Z, P=P, D=0.8)


print closed_orbit
print tline
#tline.full_tracking()
tline.full_tracking(drift_to_multi=True)
print tline
res= tline.run(xterm=0)
ftrack = res.get_all('fai')
ptrack = res.get_all('plt')

fh = open("zgoubi.plt.csv", "w")
fh.write(",".join(ptrack.dtype.names) + "\n") 
for p in ptrack:
	fh.write(",".join([str(x) for x in p]) +"\n")
#print ptrack['element_label1']

from zgoubi.lab_plot import LabPlot

print l

lp = LabPlot(l, boro=-ke_to_rigidity(10e6,ELECTRON_MASS), style={"track":{"linewidth":2, "color":"g"}})
#lp.add_tracks(ftrack=ftrack)
#lp.add_tracks(ptrack=ptrack)
lp.add_tracks(ftrack, ptrack)

#exit()
lp.draw()
#lp.draw(draw_field_points=True, field_component='z')
#lp.draw(draw_tracks=True, draw_field_midplane=True, field_component='z')
#lp.save("emma.pdf")
#lp.save("emma.svg")
lp.show()


