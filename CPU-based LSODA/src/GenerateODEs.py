import sys
import scipy.misc
import math
import numpy as np
from pylab import *
from UserString import MutableString

def fattoriale(n):
	if n==1:
		return 1
	else:
		return n*fattoriale(n-1)


class GenerateODEs():


	def __init__(self):

		self.DEFAULT_FOLDER = ""
		self.DEFAULT_OUTPUT = ""
		self.DEFAULT_NAME   = "ODE"
		self.N_A = 6.022e23
		self.c_moles_molec = 0
		self.conv = 0
		self.time_instants = []

		self.prevent_functions_to_file = True

	def convertion_constant(self, k, a1):

		a = np.zeros((len(a1), len(a1[0])), dtype=np.int32)
		k1 = np.zeros(len(k), dtype=np.float64)
		k2 = []
		for i in range(0, len(a1)):
			for j in range(0, len(a1[0])):
				a[i][j] = int(a1[i][j])

		for i in range(0, len(k)):
			value = sum(a[i][:])
			if value >= 3:
				print("There are one or more reactions with higher order to the second")
				exit(-1)
			else:
				om = 0
				for j in range(0, len(a[0])):
					if(a[i][j] >= 2):
						conv =  self.c_moles_molec
						k1[i] =  (k[i]*conv)/2.0
						om = 1
						break
					if om == 0:
						if sum(a[i][:]) > 1:
							conv = self.c_moles_molec
							k1[i] =  k[i]*conv
						else:
							k1[i] = k[i]

			k2.append(k1[i])
		return k2

	def open_files(self, folder = "", output = "", row=0, dump=True, use_cmatrix=True):

		self.functions_strings_for_lsoda = ''

		if folder!="":
			self.DEFAULT_FOLDER = folder
		if output!="":
			self.DEFAULT_OUTPUT = output

		if self.DEFAULT_FOLDER[-2:]!="//": self.DEFAULT_FOLDER+="//"

		try:
			with open(self.DEFAULT_FOLDER+"modelkind") as f:
				string = f.readlines()
				string1 = string[0].split("\n")
				if string1[0] == 'deterministic':
					self.conv = 0
				elif string1[0] == 'stochastic':
					self.conv = 1
		except:
			self.conv = 0

		if self.conv == 1:
			with open(self.DEFAULT_FOLDER+"volume") as f:
				volume = float(f.readline())
				self.c_moles_molec = volume * self.N_A;
		else:
			volume = 1
			self.c_moles_molec = volume

		with open(self.DEFAULT_FOLDER+"left_side") as f:
			linee = f.readlines()
			matrice_sinistra =  map ( lambda x: x.strip("\n").split("\t"), linee )

		with open(self.DEFAULT_FOLDER+"right_side") as f:
			linee = f.readlines()
			matrice_destra =  map ( lambda x: x.strip("\n").split("\t"), linee )
		matr_sx_interi = map(int, np.array(matrice_sinistra).flatten())
		matr_dx_interi = map(int, np.array(matrice_destra).flatten())

		differenza = map( lambda (x,y): y-x, zip(matr_sx_interi, matr_dx_interi) )

		matr_differenza = np.reshape(differenza, (len(matrice_sinistra), len(matrice_sinistra[0]) )).tolist()

		if use_cmatrix:
			with open(self.DEFAULT_FOLDER+"c_matrix") as f:
				linee = f.readlines()
				parametri = map( lambda x: float(x), linee[row].strip("\n").split() )
		else:
			with open(self.DEFAULT_FOLDER+"c_vector") as f:
				parametri = []
				for linea in f:
					parametri.append(float(linea))
		
		with open(self.DEFAULT_FOLDER+"M_0") as f:
			linee = f.readlines()
			temp = linee[0].strip("\n").split("\t")
			quantita_iniziali = map( lambda x: float(x), temp)
		try:
			with open(self.DEFAULT_FOLDER+"M_feed") as f:
				feeds = map( float, f.readline().strip("\n").split() )

		except IOError:
			if dump: print "WARNING: no feed file found"
			feeds = [0]*len(matrice_sinistra[0])
			pass

		# read time instants for sampling
		self.time_instants = loadtxt(self.DEFAULT_FOLDER+"t_vector")
		time_max = self.time_instants[-1]

		try:
			with open(self.DEFAULT_FOLDER+"alphabet") as f:
				linee = f.readlines()
				alfabeto = linee[0].strip("\n").split("\t")
		except:
			alfabeto = []
			for i in xrange(len(quantita_iniziali)):
				alfabeto.append("X" + str(i))

		try:
			with open(self.DEFAULT_FOLDER+"cs_vector") as f:
				subspecies = []
				for linea in f:
					subspecies.append(int(linea))
		except:
			subspecies = []
			for i in xrange(len(alfabeto)):
				subspecies.append(i)

		elenco_costanti = [ "k"+str(i) for i in range(len(matrice_sinistra)) ]

		quantita_iniziali = map(lambda x: float(x)/self.c_moles_molec, quantita_iniziali)
		feeds = map(lambda x: float(x)/self.c_moles_molec, feeds)
		for i in range(0, len(feeds)):
			if feeds[i] > 0:
				quantita_iniziali[i] = feeds[i]

		try:
			with open(self.DEFAULT_FOLDER+"/max_step") as f:
				self.max_integration_steps = int(f.readline())
		except:
			self.max_integration_steps=10000

		try:
			self.atolvector = loadtxt(self.DEFAULT_FOLDER+"/atol_vector")
		except:
			self.atolvector = [1e-12]*len(alfabeto)

		try:
			self.rtol = loadtxt(self.DEFAULT_FOLDER+"/rtol")
		except:
			self.rtol = 1e-6

		if dump:
			print "Stoichiometry size:", len(matrice_sinistra), "X", len(alfabeto)
			print "Absolute tolerance:", self.atolvector
			print "Relative tolerance:", self.rtol
			print "Maximum number of steps:", self.max_integration_steps

		if self.prevent_functions_to_file:

			listaIntermedia = []
			self.functions_strings_for_lsoda = ''
			listaIntermedia.append("# differential equation output created by stoc2det\n")

			listaIntermedia.append("\ny0 = "+str(quantita_iniziali))

			listaIntermedia.append("\ntEnd = "+str(time_max))

			listaIntermedia.append("\nalfabeto = "+str(alfabeto)+"\n")

			listaIntermedia.append("\nvolume = "+str(volume)+"\n")

			if self.conv == 1:
				parametri1 = self.convertion_constant(parametri, matrice_sinistra)
				listaIntermedia.append("\nparametri = "+str(parametri1)+"\n\n")
			else:
				listaIntermedia.append("\nparametri = "+str(parametri)+"\n\n")

			listaIntermedia.append("def "+self.DEFAULT_NAME+" (y,t,p):\n")

			listaIntermedia.append("\n\tdy=zeros(["+str(len(alfabeto))+"])\n")

			service_array = [[]]*len(alfabeto)
			for specie in range(len(alfabeto)):

				sub_array = []

				if dump: print ""
				listaIntermedia.append("\n")

				if dump: print "d"+alfabeto[specie]+"/dt=",
				listaIntermedia.append("\t"+"dy["+str(specie)+"]=")

				qualcosa = False

				if feeds[specie]!=0:
					listaIntermedia.append("0")
					if dump: print "0",
					continue


				for riga in range(len(elenco_costanti)):

					if matrice_sinistra[riga][specie]!='0':

						if matr_differenza[riga][specie]<0:

							if dump: print "-"+elenco_costanti[riga],
							listaIntermedia.append(str(matr_differenza[riga][specie])+"*p["+str(riga)+"]")
							qualcosa = True
							sub_array.append("+p["+str(riga)+"]*"+str(matr_differenza[riga][specie]))

							for subs, spec in enumerate(alfabeto):
								if matrice_sinistra[riga][subs]=='1':
									if dump: print "* ["+alfabeto[subs]+"]",
									listaIntermedia.append("*y["+str(subs)+"]")
									sub_array.append(str(subs))

								elif matrice_sinistra[riga][subs]=='2':
									if dump: print "* ["+alfabeto[subs]+"]^2",
									listaIntermedia.append("*y["+str(subs)+"]*y["+str(subs)+"]")
									sub_array.append(str(subs))
									sub_array.append(str(subs))

				for riga in range(len(elenco_costanti)):

					if matrice_destra[riga][specie]!='0':

						if matr_differenza[riga][specie]>0:

							if dump: print "+"+elenco_costanti[riga],
							listaIntermedia.append("+"+str(matr_differenza[riga][specie])+"*p["+str(riga)+"]")
							qualcosa = True
							sub_array.append("+p["+str(riga)+"]*"+str(matr_differenza[riga][specie]))

							for subs, spec in enumerate(alfabeto):
								if matrice_sinistra[riga][subs]=='1':
									if dump: print "* ["+alfabeto[subs]+"]",
									listaIntermedia.append("*y["+str(subs)+"]")								
									sub_array.append(str(subs))

								elif matrice_sinistra[riga][subs]=='2':
									if dump: print "* ["+alfabeto[subs]+"]^2",
									listaIntermedia.append("*y["+str(subs)+"]*y["+str(subs)+"]")
									sub_array.append(str(subs))
									sub_array.append(str(subs))


				if not qualcosa:
					listaIntermedia.append("0")

				sub_array.append("!")
				service_array[specie] = sub_array

			listaIntermedia.append("\n\n\treturn dy\n\n\n")

			varray = [""] * len(alfabeto)
			matrice_jacobiana = []
			for i in range(len(alfabeto)):
				matrice_jacobiana.append( varray[:] )

			# per tutte le F (righe)
			for n,fun in enumerate(service_array):
				appoggio = fun
				for m, s in enumerate(alfabeto):
					collect = []
					for token in appoggio:
						if token[0] == "-" or token[0] == "+" or token[0] == "!":
							if len(collect)>0:
								if collect.count(str(m))==1:
									collect.remove(str(m))
									temp = [collect[0]]
									temp2 = map( lambda x: "y["+x+"]", collect[1:] )
									temp.extend(temp2)
									matrice_jacobiana[n][m] = matrice_jacobiana[n][m]  + "*".join(temp)
								elif collect.count(str(m))==2:
									collect.remove(str(m))
									collect.remove(str(m))
									collect.append("2*" + str(m))
									temp = [collect[0]]
									for sc in collect[1:]:
										if sc[0:2]=="2*":
											temp.append("2*y["+str(m)+"]")
										else:
											temp.append("y["+str(m)+"]")

									matrice_jacobiana[n][m] = matrice_jacobiana[n][m]  + "*".join(temp)
							collect = list()
							collect.append(token)
						else:
							if token!="!":
								collect.append(token)

					if matrice_jacobiana[n][m]=="":
						matrice_jacobiana[n][m] = "0"

			listaIntermedia.append("def JAC(y,t,p):\n")

			listaIntermedia.append("\n\tdy=zeros(["+str(len(alfabeto))+","+str(len(alfabeto))+"])\n")

			for n, r in enumerate(matrice_jacobiana):
				for m, c in enumerate(r):
					listaIntermedia.append("\tdy["+str(n)+"]["+str(m)+"] ="+c+" \n")

			listaIntermedia.append("\treturn dy\n")

			self.functions_strings_for_lsoda = "".join(listaIntermedia)
		else:			

			with open(self.DEFAULT_OUTPUT, 'w') as f:

				f.write("# differential equation output created by stoc2det\n")

				f.write("\ny0 = "+str(quantita_iniziali))

				f.write("\ntEnd = "+str(time_max))

				f.write("\nalfabeto = "+str(alfabeto)+"\n")

				f.write("\nvolume = "+str(volume)+"\n")

				if self.conv == 1:
					parametri1 = self.convertion_constant(parametri, matrice_sinistra)
					f.write("\nparametri = "+str(parametri1)+"\n\n")
				else:
					f.write("\nparametri = "+str(parametri)+"\n\n")

				f.write("def "+self.DEFAULT_NAME+" (y,t,p):\n")

				f.write("\n\tdy=zeros(["+str(len(alfabeto))+"])\n")

				service_array = [[]]*len(alfabeto)



				for specie in range(len(alfabeto)):

					sub_array = []

					if dump: print ""
					f.write("\n")

					if dump: print "d"+alfabeto[specie]+"/dt=",
					f.write("\t"+"dy["+str(specie)+"]=")

					qualcosa = False

					if feeds[specie]!=0:
						f.write("0")
						if dump: print "0",
						continue


					for riga in range(len(elenco_costanti)):

						if matrice_sinistra[riga][specie]!='0':

							if matr_differenza[riga][specie]<0:

								if dump: print "-"+elenco_costanti[riga],
								f.write(str(matr_differenza[riga][specie])+"*p["+str(riga)+"]")
								qualcosa = True
								sub_array.append("+p["+str(riga)+"]*"+str(matr_differenza[riga][specie]))

								for subs, spec in enumerate(alfabeto):
									if matrice_sinistra[riga][subs]=='1':
										if dump: print "* ["+alfabeto[subs]+"]",
										f.write("*y["+str(subs)+"]")
										sub_array.append(str(subs))

									elif matrice_sinistra[riga][subs]=='2':
										if dump: print "* ["+alfabeto[subs]+"]^2",
										f.write("*y["+str(subs)+"]*y["+str(subs)+"]")
										sub_array.append(str(subs))
										sub_array.append(str(subs))


					for riga in range(len(elenco_costanti)):

						if matrice_destra[riga][specie]!='0':

							if matr_differenza[riga][specie]>0:

								if dump: print "+"+elenco_costanti[riga],
								f.write("+"+str(matr_differenza[riga][specie])+"*p["+str(riga)+"]")
								qualcosa = True
								sub_array.append("+p["+str(riga)+"]*"+str(matr_differenza[riga][specie]))

								for subs, spec in enumerate(alfabeto):
									if matrice_sinistra[riga][subs]=='1':
										if dump: print "* ["+alfabeto[subs]+"]",
										f.write("*y["+str(subs)+"]")
										sub_array.append(str(subs))

									elif matrice_sinistra[riga][subs]=='2':
										if dump: print "* ["+alfabeto[subs]+"]^2",
										f.write("*y["+str(subs)+"]*y["+str(subs)+"]")
										sub_array.append(str(subs))
										sub_array.append(str(subs))

					if not qualcosa:
						f.write("0")

					sub_array.append("!")
					service_array[specie] = sub_array


				f.write("\n\n\treturn dy\n\n\n")

				varray = [""] * len(alfabeto)
				matrice_jacobiana = []
				for i in range(len(alfabeto)):
					matrice_jacobiana.append( varray[:] )

				# per tutte le F (righe)
				for n,fun in enumerate(service_array):
					appoggio = fun
					for m, s in enumerate(alfabeto):
						collect = []
						for token in appoggio:
							if token[0] == "-" or token[0] == "+" or token[0] == "!":
								if len(collect)>0:
									if collect.count(str(m))==1:
										collect.remove(str(m))
										temp = [collect[0]]
										temp2 = map( lambda x: "y["+x+"]", collect[1:] )
										temp.extend(temp2)
										matrice_jacobiana[n][m] = matrice_jacobiana[n][m]  + "*".join(temp)
									elif collect.count(str(m))==2:
										collect.remove(str(m))
										collect.remove(str(m))
										collect.append("2*" + str(m))
										temp = [collect[0]]
										for sc in collect[1:]:
											if sc[0:2]=="2*":
												temp.append("2*y["+str(m)+"]")
											else:
												temp.append("y["+str(m)+"]")
										matrice_jacobiana[n][m] = matrice_jacobiana[n][m]  + "*".join(temp)
								collect = list()
								collect.append(token)
							else:
								if token!="!":
									collect.append(token)

						if matrice_jacobiana[n][m]=="":
							matrice_jacobiana[n][m] = "0"

				f.write("def JAC(y,t,p):\n")

				f.write("\n\tdy=zeros(["+str(len(alfabeto))+","+str(len(alfabeto))+"])\n")

				for n, r in enumerate(matrice_jacobiana):
					for m, c in enumerate(r):
						f.write( "\tdy["+str(n)+"]["+str(m)+"] ="+c+" \n" )

				f.write("\treturn dy\n")

		return subspecies
