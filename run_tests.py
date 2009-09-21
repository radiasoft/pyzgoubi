#!/usr/bin/env python
import tempfile
import shutil
import os
import sys
import subprocess
import time

#install pyzoubi to temp folder
install_dir = tempfile.mkdtemp(prefix='/tmp/pyzgoubi_test_inst_')
print "installing to", install_dir
#install_res = os.system("python setup.py install --prefix=%s"%install_dir)
install_res = subprocess.Popen(["./setup.py", "install", "--prefix=%s"%install_dir], stdout=subprocess.PIPE)

for line in install_res.communicate()[0].split('\n'):
	if line.startswith('alias'):
		pyzgoubi_cmd = line.partition('=')[2]
		pyzgoubi_cmd = pyzgoubi_cmd.strip('"')
		
if install_res.returncode != 0:
	print "ERROR: install failed"
	sys.exit(1)



# move to another temp dir for running tests
run_dir = tempfile.mkdtemp(prefix='/tmp/pyzgoubi_test_run_')
print "running tests from", run_dir
os.chdir(run_dir)

test_dir = os.path.join(install_dir, 'share', 'pyzgoubi','test')

number_of_tests = len(os.listdir(test_dir))
tests_run = 0
tests_sucess = []
tests_fail = []
tot_time = 0

for test_file in os.listdir(test_dir):
	full_test_file = os.path.join(test_dir, test_file)
	tests_run += 1
	print "running test %s, %d of %d"%(test_file, tests_run, number_of_tests)
	command = pyzgoubi_cmd + " " + full_test_file
	t0 = time.time()
	result = os.system(command)
	t1 = time.time()
	t = t1 - t0
	if result == 0:
		print "PASS:", test_file
		tests_sucess.append(test_file)
	else:
		print "FAIL:", test_file
		tests_fail.append(test_file)
	print "Took %s sec"%t
	tot_time += t


print ""
print "Summary:"
print "Ran %d tests"%number_of_tests
print "Pass %d"%len(tests_sucess)
print "Fail %d"%len(tests_fail)
if len(tests_fail) != 0:
	print "Failed tests:"
	for t in tests_fail:
		print t

print "Took %s sec"%tot_time



shutil.rmtree(install_dir)
shutil.rmtree(run_dir)
