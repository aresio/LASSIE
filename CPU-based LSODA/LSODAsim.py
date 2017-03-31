#!/usr/bin/python
import sys
sys.path.append("src")

from LSODA import *
from  GenerateODEs import *
import numpy as np
import os
import time


percorso = ""
percorso_output =""
dump = False

if __name__ == '__main__':

	t0 = time.time()

	if len(sys.argv)>1:
		percorso = sys.argv[1]
	if len(sys.argv)>2:
		percorso_output = sys.argv[2]
	if len(sys.argv)>3:
		dump = True if sys.argv[3]=="-v" else False

	if dump:
		print " * Opening project in folder", percorso
		print " * Saving output file in folder", percorso_output

	try:
		os.makedirs(percorso_output)
	except:
		if dump: print "WARNING:", percorso_output, "already exists"

	# generation of ODEs and Jacobian
	gen = GenerateODEs()
	subspecies = gen.open_files(folder=percorso, output=percorso+"ODE.py", dump=dump, use_cmatrix=False)
		
	# Solver 
	l = LSODER()
	l.load_data_nofile(gen.functions_strings_for_lsoda)
	
	l.time_instants = gen.time_instants
	l.atolvector = gen.atolvector
	l.rtol = gen.rtol
	l.max_steps = gen.max_integration_steps

	ret = l.run(subspecies)

	new_matrix = [[]]*len(l.time_instants)

	for n, a in enumerate(new_matrix):
		new_matrix[n] = [ l.time_instants[n] ]
		new_matrix[n].extend(ret[n])

	np.savetxt(percorso_output + os.sep + "ODESol", new_matrix, fmt="%.8f", delimiter = "\t")
	print " * Dynamics saved in folder", percorso_output+os.sep+"ODESol"

	t1 = time.time()
	total = t1-t0
	print "Simulation time:", total
