#!/usr/bin/python
import sys

try:
	sys.path.append("src")
	from LSODA import *
	from  GenerateODEs import *
except:
	sys.path.append("lsodaSIM/src")
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
	dump = False
	fit  = False

	if len(sys.argv)>1:
		percorso = sys.argv[1]
	if len(sys.argv)>2:
		percorso_output = sys.argv[2]
	# if len(sys.argv)>3:
	# 	dump = True if sys.argv[3]=="-v" else False
	if len(sys.argv)>3:
		for i in range(3, len(sys.argv)):
			if sys.argv[i]=="-v":
				dump = True
			elif sys.argv[i]=="-f":
				fit = True
	if fit:
		dump = False
		try:
			ts_matrix = np.loadtxt(percorso+os.sep+"ts_matrix")
		except:
			print " * There is not ts_matrix file containing DTTS in", percorso


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
	
	if gen.time_instants[0] != 0:
		tmp = np.zeros(len(gen.time_instants)+1)
		tmp[0] = 0
		tmp[1:] = gen.time_instants
		l.time_instants = tmp
	else:
		l.time_instants = gen.time_instants

	l.atolvector = gen.atolvector
	l.rtol = gen.rtol
	l.max_steps = gen.max_integration_steps

	ret = l.run()

	new_matrix = [[]]*len(l.time_instants)

	for n, a in enumerate(new_matrix):
		new_matrix[n] = [ l.time_instants[n] ]
		new_matrix[n].extend(ret[n])

	finalMatrix = np.zeros((len(gen.time_instants), len(subspecies) + 1))

	if gen.time_instants[0] != 0:
		for i in xrange(len(gen.time_instants)):
			finalMatrix[i][0] = new_matrix[i+1][0]

		for i in xrange(len(gen.time_instants)):
			for j in range(0, len(subspecies)):
				idx = subspecies[j]
				finalMatrix[i][j+1] = new_matrix[i+1][idx+1]
	else:
		for i in xrange(len(gen.time_instants)):
			finalMatrix[i][0] = new_matrix[i][0]
		for i in xrange(len(gen.time_instants)):
			for j in range(0, len(subspecies)):
				idx = subspecies[j]
				finalMatrix[i][j+1] = new_matrix[i][idx+1]

	if fit:
		acc = sys.float_info.max
		try:
			acc = 0
			for i in range(0, finalMatrix.shape[0]):
				for j in range(1, finalMatrix.shape[1]):
					valueTS  = ts_matrix[i,j]
					valueSIM = finalMatrix[i,j]
					if valueTS != 0:
						if ts_matrix[i,j] < 1e-16:
							valueTS = 1e-16
						
						if finalMatrix[i,j] < 1e-16:
							valueSIM = 1e-16
						acc += abs(valueTS - valueSIM) / valueTS
		except:
			acc = sys.float_info.max

		print acc
	
	else:
		np.savetxt(percorso_output + os.sep + "ODESol", finalMatrix, fmt="%.8e", delimiter = "\t")
	
	t1 = time.time()
	total = t1-t0

	if dump:
		print "* Dynamics saved in folder", percorso_output+os.sep+"ODESol"
		print "Simulation time:", total
