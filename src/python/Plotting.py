__author__ = 'andreatangherloni'
import numpy as np
import math
import matplotlib
import matplotlib.pylab as plt
from pylab import *
from matplotlib import rc
from matplotlib import rcParams
import matplotlib.ticker as mtick
import os


def main():

	rc('xtick', labelsize=24)
	rc('ytick', labelsize=24)
	rcParams['font.family'] = 'serif'
	rcParams['font.sans-serif'] = ['Times New Roman']
	rcParams['xtick.major.pad']='16'
	rcParams['ytick.major.pad']='16'

	folder = sys.argv[2]
	path = folder + "/output/Solution"
	
	with open(path, 'r') as f:
		odes = np.genfromtxt(f, skip_header=0)

	time = np.zeros(odes.shape[0],  dtype=np.float64)

	f = open(sys.argv[1], "r")
	subspec = []
	
	for line in f:
		subspec.append(line)
	f.close()
	
	subspecies = []
	for i in range(0, len(subspec)):
		string = subspec[i].split("\n")
		subspecies.append(int(string[0]))


	ode1 = np.transpose(odes)
	ode = ode1[1::][:]
	time = ode1[0][:]

	alphabet1 = sys.argv[3]
	alphabet = []

	if alphabet1 == "no":
		for i in range(0, len(subspec)):
			alphabet.append("X" + str(i))
	else:
		f = open(alphabet1, "r")
		string = []
		for line in f:
			string.append(line)
	f.close()
	alphabet = string[0].split("\t")

	verb = sys.argv[4]
	verbose = int(verb)
	if(verbose):
		print("\nPlotting solution....");
	
	plt.set_cmap('winter')
	cmap = plt.get_cmap()
	plt.figure(figsize=(16, 9))
	plt.ylabel('Molecular concentration', fontsize = 26)
	plt.xlabel('Time [s]', fontsize = 26)

	# colors = [0.0, 0.2, 0.4, 0.6]
	
	# for i in range(0, ode.shape[0]):
	# 	#s = alphabet[subspecies[i]]
	# 	s = "$X_" + str(i+1) + "$"
	# 	c = cmap(colors[i])
	# 	plt.plot(time, ode[i][:], color=c, linewidth = 2.0, label=s)

	# path = folder + "/image/all.pdf"
	# lg = plt.legend(loc = 1, numpoints=1, borderaxespad=0, fontsize = 22)
	# lg.draw_frame(False)
	# plt.savefig(path, bbox_inches = 'tight', additional_artists=lg, dpi=300)
	
	for i in range(0, ode.shape[0]):
		s = alphabet[subspecies[i]]
		path = folder + "/image/" + s + ".pdf"
		plt.figure(255)
		plt.figure(figsize=(25, 15))
		plt.figure(figsize=(16, 9))
		plt.ylabel('Molecular concentration', fontsize = 26)
		plt.xlabel('Time [s]', fontsize = 26)
		dim_y = np.arange(min(ode[i][:]), max(ode[i][:]))
		plt.plot(time, ode[i][:], 'b', linewidth = 2.0, label=s)
		#plt.savefig(path, bbox_inches = 'tight', additional_artists=plt.legend(loc = 0, numpoints=1, borderaxespad=0.5,  fontsize = 18), dpi=92)
		lg = plt.legend(loc = 0, numpoints=1, borderaxespad=0, fontsize = 22)
		lg.draw_frame(False)
		plt.savefig(path, bbox_inches = 'tight', additional_artists=lg, dpi=300)
		plt.close('all')

if __name__ == "__main__":
	main()
