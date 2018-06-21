emma_cell = Line('emma')
xpas = (10,10,10)
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
#field
fb = -6.695 * fr * T
db = 4.704 * dr * T


emma_cell.add(DRIFT('ld', XL=ld/cm/2))
emma_cell.add(CHANGREF(ALE=angle))
emma_cell.add(CHANGREF(YCE=d_offset/cm))
emma_cell.add(QUADRUPO('defoc', XL=dq/cm, R_0=dr/cm, B_0=db/kgauss, XPAS=xpas, KPOS=1))
emma_cell.add(CHANGREF(YCE=-d_offset/cm))
emma_cell.add(DRIFT('sd', XL=sd/cm))
emma_cell.add(CHANGREF(YCE=f_offset/cm))
emma_cell.add(QUADRUPO('foc', XL=fq/cm, R_0=fr/cm, B_0=fb/kgauss, XPAS=xpas, KPOS=1))
emma_cell.add(CHANGREF(YCE=-f_offset/cm))
emma_cell.add(DRIFT('ld', XL=ld/cm/2))


data = gcp.get_cell_properties(cell=emma_cell, min_ke=10e6, max_ke=20e6, ke_steps=11, particle='e')

mkdir_p("plots")

gcp.plot_cell_tracks(cell=emma_cell, data=data, particle='e', output_file='plots/47-gcp-emma-plot_1.pdf')

gcp.plot_cell_tracks(cell=emma_cell, data=data, particle='e', output_file='plots/47-gcp-emma-plot_2.pdf', draw_field_midplane=True, min_y=-10, max_y=10, y_steps=100)

