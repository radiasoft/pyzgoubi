# -*- coding: utf-8 -*-
print "running emma example"

emma = Line('emma')

xpas = (20,20,20)

cells = 42
angle = 360/cells
d_offset = -34.048 * mm
f_offset = -7.514 * mm

#lengths
ld = 210 * mm
sd = 50 * mm

fq = 58.782 * mm
dq = 75.699 * mm

# quad radius
fr = 37 * mm
dr = 53 * mm

fb = -6.695 * fr * T
db = 4.704 * dr * T

ob = OBJET2()
emma.add(ob)

emma.add(ELECTRON())


emma.add(DRIFT('ld', XL=ld*cm_/2))
emma.add(CHANGREF(ALE=angle))

emma.add(CHANGREF(YCE=d_offset*cm_))
emma.add(QUADRUPO('defoc', XL=dq*cm_, R_0=dr*cm_, B_0=db*kgauss_, XPAS=xpas, IL=2))
emma.add(CHANGREF(YCE=-d_offset*cm_))

emma.add(DRIFT('sd', XL=sd*cm_))

emma.add(CHANGREF(YCE=f_offset*cm_))
emma.add(QUADRUPO('foc', XL=fq*cm_, R_0=fr*cm_, B_0=fb*kgauss_, XPAS=xpas, IL=2))
emma.add(CHANGREF(YCE=-f_offset*cm_))

emma.add(DRIFT('ld', XL=ld*cm_/2))

emma.add(FAISCNL(FNAME='zgoubi.fai'))

#add REBELOTE like this to allow a later modification 
reb=REBELOTE(K=99, NPASS=1)
emma.add(reb)


emma.add(END())

rigidity = ke_to_rigidity(10e6, 0.51099892e6)
ob.set(BORO=-rigidity)
ob.add(Y=0, T=0, D=1)

r = emma.run(xterm = False)

#first step is to find closed orbit
closedorb_YTZP = find_closed_orbit(emma, init_YTZP=[1,2,3,4], tol=1e-6)

#Change from OBJET2 to OBJET5 so that MATRIX can compute transfer matrix etc
ob5 = OBJET5()
emma.replace(ob, ob5)
ob5.set(BORO=-rigidity)
ob5.set(PY=1e-4,PT=1e-3,PZ=1e-4,PP=1e-3,PX=1e-3,PD=1e-3)
ob5.set(YR=closedorb_YTZP[0],TR=closedorb_YTZP[1],ZR=closedorb_YTZP[2],PR=closedorb_YTZP[3],DR=1.0)
matrix=MATRIX(IORD=1,IFOC=11)
emma.replace(reb,matrix)
#run zgoubi to find tune etc.
r = emma.run(xterm = False)

#find tune calculated by MATRIX over this periodic cell
tune = r.get_tune()

#get twiss parameters at end of cell, returns [beta_y,alpha_y,gamma_y,disp_y,disp_py,beta_z,alpha_z,gamma_z,disp_z,disp_pz]
twissparam = r.get_twiss_parameters()
betayz = [twissparam[0],twissparam[5]]
alphayz = [twissparam[1],twissparam[6]]
gammayz = [twissparam[2],twissparam[7]]

print "beta Y,Z at end of cell ",betayz
print "alpha Y,Z at end of cell ",alphayz
print "gamma Y,Z at end of cell ",gammayz

#get_twiss_profiles has format [s_coord, label, mu_y, beta_y, alpha_y, gamma_y, disp_y, disp_py, mu_z,beta_z, alpha_z, gamma_z, disp_z, disp_pz]
twiss_profiles = get_twiss_profiles(emma,'twiss_profiles.txt')

#Note - Could specify twiss parameters at beginning of cell
#twiss_profiles = get_twiss_profiles(emma,'twiss_profiles.txt',input_twiss_parameters = twissparam)

#extract s coordinate and beta_y from twiss_profiles and create figure beta_y_profile.png
s = twiss_profiles[0]
beta_y = twiss_profiles[3]
beta_z = twiss_profiles[9]
plot_data_xy_multi(s,[beta_y,beta_z], 'beta_profiles', labels=["beta profile","s [m]","beta_y [m]"],style=['k+','b.'],legend=['beta_y','beta_z'])

#plot horizontal phase advance mu_y
mu_y = twiss_profiles[2]
plot_data_xy_multi(s, mu_y, 'mu_y_profile', labels=["mu_y profile","s [m]","mu_y [rad]"],style='b+')

#plot dispersion in the horizontal plane
disp_y = twiss_profiles[6]
plot_data_xy_multi(s, disp_y, 'disp_profile', labels=["D profile","s [m]","D [m]"],style='b+')
