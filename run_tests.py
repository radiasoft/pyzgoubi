#!/usr/bin/env python
import tempfile
import shutil
import os
import sys
import subprocess
import time
import datetime


log_file_path = "test-" + datetime.datetime.today().strftime("%Y%m%d-%H%M") + ".log"
log = open(log_file_path, "w")
print "writing test log to:", log_file_path

#install pyzoubi to temp folder
install_dir = tempfile.mkdtemp(prefix='/tmp/pyzgoubi_test_inst_')
print "installing to", install_dir
print >> log, "installing to", install_dir
#install_res = os.system("python setup.py install --prefix=%s"%install_dir)
#subprocess.Popen(["./setup.py", "clean", "--all"])
clean_proc = subprocess.Popen(["./setup.py", "clean", "--all"], stdout=log, stderr=subprocess.STDOUT)
clean_proc.wait()
install_res = subprocess.Popen(["./setup.py", "install", "--prefix=%s"%install_dir], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

for line in install_res.communicate()[0].split('\n'):
	print >>log, line
	if line.startswith('alias'):
		pyzgoubi_cmd = line.partition('=')[2]
		pyzgoubi_cmd = pyzgoubi_cmd.strip('"')
		
if install_res.returncode != 0:
	print "ERROR: install failed"
	print >> log, "ERROR: install failed"
	print "If there were permission errors, try running 'sudo ./setup clean --all' and 'sudo rm install.log'"
	print >> log,"If there were permission errors, try running 'sudo ./setup clean --all' and 'sudo rm install.log'"
	sys.exit(1)

print
print >>log
log.flush()
proc = subprocess.Popen(pyzgoubi_cmd+" --version", shell=True, stderr=subprocess.STDOUT, stdout=log)
proc.wait()

# move to another temp dir for running tests
run_dir = tempfile.mkdtemp(prefix='/tmp/pyzgoubi_test_run_')
print "running tests from", run_dir
print >> log,"running tests from", run_dir
os.chdir(run_dir)

test_dir = os.path.join(install_dir, 'share', 'pyzgoubi','test')

number_of_tests = len(os.listdir(test_dir))
tests_run = 0
tests_sucess = []
tests_fail = []
tot_time = 0

for test_file in os.listdir(test_dir):
	print
	print >>log, "\n", "="*40
	log.flush()
	if len(sys.argv) > 1:
		if test_file not in sys.argv:
			print "skipping", test_file
			print >> log, "skipping", test_file
			continue
	full_test_file = os.path.join(test_dir, test_file)
	tests_run += 1
	print "running test %s, %d of %d"%(test_file, tests_run, number_of_tests)
	print >> log,"running test %s, %d of %d"%(test_file, tests_run, number_of_tests)
	log.flush()
	command = pyzgoubi_cmd + " " + full_test_file
	t0 = time.time()
#	result = os.system(command)
	print command
	proc = subprocess.Popen(command, shell=True, stderr=subprocess.STDOUT, stdout=log)
	proc.wait()
	result = proc.returncode
	t1 = time.time()
	t = t1 - t0
	log.flush()
	print
	print >>log, "\n", "="*40
	log.flush()
	if result == 0:
		print "PASS:", test_file
		print >> log,"PASS:", test_file
		tests_sucess.append(test_file)
	else:
		print "FAIL:", test_file
		print >> log,"FAIL:", test_file
		tests_fail.append(test_file)
	print "Took %s sec"%t
	print >> log,"Took %s sec"%t
	tot_time += t


print "\nSummary:"
print >> log,"\nSummary:"
print "Ran %d tests"%number_of_tests
print >> log,"Ran %d tests"%number_of_tests
print "Pass %d"%len(tests_sucess)
print >> log,"Pass %d"%len(tests_sucess)
print "Fail %d"%len(tests_fail)
print >> log,"Fail %d"%len(tests_fail)
if len(tests_fail) != 0:
	print "Failed tests:"
	print >> log,"Failed tests:"
	for t in tests_fail:
		print t

print "Took %s sec"%tot_time
print >> log,"Took %s sec"%tot_time



shutil.rmtree(install_dir)
shutil.rmtree(run_dir)

print "Test written log to:", log_file_path
