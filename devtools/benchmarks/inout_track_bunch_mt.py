"""Compare times with different numbers of threads
Low overheads - Lots of steps in zgoubi, so zgoubi time dominates
High overheads - Few steps, so pyzgoubi time dominates

"""
import time
mass = PROTON_MASS
energy = 1e6
l_max_cpu = 3 # log of max cpu 3 -> 2**3 = 8
n_particles = 1e3
rep = 5

output = ""

b_orig = Bunch.gen_halo_x_xp_y_yp(n_particles, 1e-3, 1e-3, 4, 5, 1e-3, 2e-2, ke=energy, mass=mass, charge=1)
b_orig.particles()[0]['D'] = 1.1

for low_overhead in [True, False]:
	line_seg = Line("lineseg")

	line_seg.add(PROTON())
	if low_overhead:
		line_seg.add(QUADRUPO(XL=1, R_0=5, B_0=0.1, XPAS=(50,5000,50)))
	else:
		line_seg.add(QUADRUPO(XL=1, R_0=5, B_0=0.1, XPAS=(1,1,1)))

	itimes = []
	#t0 = time.time()
	for rep_x in xrange(rep):
		ti0 = time.time()
		st_bunch = line_seg.track_bunch_mt(b_orig, n_threads=1, max_particles=1e3)
		ti1 = time.time()
		itimes.append(ti1 - ti0)
	#t1 = time.time()
	#st_time = (t1-t0)/rep
	st_time = min(itimes)
	st_end = st_bunch.particles()[['Y', 'P', 'Z', 'T', 'D']].view(float).reshape([-1, 5])

	n_threads = []
	times = []
	for x in xrange(1,l_max_cpu + 1):
		n_thread = 2 ** x
		itimes = []
		#t0 = time.time()
		for rep_x in xrange(rep):
			ti0 = time.time()
			mt_bunch = line_seg.track_bunch_mt(b_orig, n_threads=n_thread, max_particles=1e3)
			ti1 = time.time()
			itimes.append(ti1 - ti0)
		#t1 = time.time()
		#mt_time = (t1-t0) / rep
		#times.append(mt_time)
		times.append(min(itimes))
		n_threads.append(n_thread)

		mt_end = mt_bunch.particles()[['Y', 'P', 'Z', 'T', 'D']].view(float).reshape([-1, 5])
		errors = abs((mt_end - st_end) / numpy.maximum(mt_end, st_end))
		print "%r" % errors[0][0]
		print "mean errors in YTZPD: single vs ", n_thread ," thread"
		print errors.mean(0)
		assert(numpy.all(errors.mean(0) < [1e-16, 2e-16, 1e-16, 2e-16, 1e-16])), "error to big"


	if low_overhead:
		output += "Low overhead (lots of time in zgoubi)\n"
	else:
		output += "\nHigh overhead (short time spend in zgoubi)\n"
	output += "threads\ttime\tspeedup\n"
	output += "1\t%.4f\t1\n" % st_time

	for x in xrange(len(n_threads)):
		output +=  "%d\t%.4f\t%.4f\n" % (n_threads[x], times[x], st_time/times[x])

print output
