#from gev import GEV
from PyQt4 import QtGui, QtCore, uic
from numpy import *
import sys
from subprocess import check_output, Popen, call
from platform import platform
from pylab import *
from SBML2BSW import SBMLloader

import resources_rc

class MyDialog(QtGui.QDialog):
	def __init__(self):
		super(MyDialog, self).__init__()
		uic.loadUi('about-lassie.ui', self)




class MyWindow(QtGui.QMainWindow):

	def __init__(self):
		super(MyWindow, self).__init__()
		uic.loadUi('lassie-gui.ui', self)

		self.input_valid = False
		self.output_valid = False

		self.show()
		self.figure = figure()
		self._about = MyDialog()
		self._about.setModal(True)
		#self.populate_model_data("C:\\Users\\aresio\\Documents\\PythonCode\\ModelliTest\\MM")
		
	def show_help(self):
		self._about.show()

	def import_SBML(self):
		file = str(QtGui.QFileDialog.getOpenFileName(self, "Select Directory"))
		if file!="":
			SL = SBMLloader(file)  
    		SL.convert_to_biosimware("converted_sbml")
    		self._open_model("converted_sbml")

	def select_all(self):
		for x in xrange(self.speciestable.rowCount()):
			self.speciestable.item(x,2).setCheckState(QtCore.Qt.Checked)
	

	def deselect_all(self):
		for x in xrange(self.speciestable.rowCount()):
			self.speciestable.item(x,2).setCheckState(QtCore.Qt.Unchecked)

	def is_valid_input_directory(self, path):
		return True

	def show_warning(self, text):
		x  = QtGui.QMessageBox()
		x.setWindowTitle("Warning")
		x.setText(text)
		x.setIcon(QtGui.QMessageBox.Warning)
		x.exec_()

	def populate_model_data(self, path):

		self.input_valid = False

		try:
			A = loadtxt(path+"/left_side") 
		except:
			self.show_warning("Cannot open model directory")
			self.disable_simulation(); return
		reactions = len(A)
		species = len(A[0])		

		try:
			with open(path+"/alphabet") as fi:
				species_names = fi.readline().split()
		except:				
			print "Cannot open", path+"/alphabet"
			species_names = ["S_"+str(x) for x in xrange(species)]
		print species_names

		try:
			species_amounts = loadtxt(path+"/M_0")
			self.populate_species(species_names, species_amounts)
		except:
			self.show_warning("Error loading initial quantities")
			self.disable_simulation(); return

		try:
			left = loadtxt(path+"/left_side")
			right = loadtxt(path+"/right_side")
			constants = loadtxt(path+"/c_vector")
			self.populate_reactions(species_names, left, right, constants)
		except:
			self.show_warning("Error loading reactions")
			self.disable_simulation(); return

		self.input_valid = True

		if self.is_everything_ready(): self.enable_simulation()

	def is_everything_ready(self):
		return self.input_valid and self.output_valid 


	def enable_simulation(self): self.simulatebutton.setEnabled(True)
	def disable_simulation(self): self.simulatebutton.setEnabled(False)

	def populate_species(self, names, amounts):
		self.speciestable.setRowCount(len(names))
		self.speciestable.setColumnCount(3)
		for i in xrange(len(names)):
			n = names[i]
			a = amounts[i] 
			w_n = QtGui.QTableWidgetItem(n)
			w_a = QtGui.QTableWidgetItem(str(a))

			w_c = QtGui.QTableWidgetItem()
			w_c.setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled)
			w_c.setCheckState(QtCore.Qt.Checked)       
			w_c.setTextAlignment(QtCore.Qt.AlignCenter)
			
			self.speciestable.setItem(i, 0, w_n)
			self.speciestable.setItem(i, 1, w_a)
			self.speciestable.setItem(i, 2, w_c)

	def populate_reactions(self, names, left, right, K):
		self.reactionstable.setRowCount(len(K))
		self.reactionstable.setColumnCount(2)
		left_side_reactions = []
		right_side_reactions = []
		for i, reactants in enumerate(left):
			react_list = []
			for j in xrange(len(reactants)):
				if reactants[j]>0: 
					if reactants[j]==1:
						react_list.append(names[j])
					else:
						react_list.append(str(reactants[j])+"*"+names[j])
			left_side_reactions.append( "+".join(react_list) )
		for i, products in enumerate(right):
			prod_list = []
			for j in xrange(len(products)):
				if products[j]>0: 
					if products[j]==1:
						prod_list.append(names[j])
					else:
						prod_list.append(str(products[j])+"*"+names[j])
			right_side_reactions.append( "+".join(prod_list) )
		for k in xrange(len(left_side_reactions)):
		#for a,b in zip(left_side_reactions, right_side_reactions):
			a,b=left_side_reactions[k], right_side_reactions[k]
			reazione = a+"->"+b		
			w_r = QtGui.QTableWidgetItem(reazione)
			w_k = QtGui.QTableWidgetItem(str(K[k]))
			self.reactionstable.setItem(k, 0, w_r)
			self.reactionstable.setItem(k, 1, w_k)

	def open_model(self):
		last_dir = "."
		try:
			with open(".last_dir") as fi:
				last_dir = fi.readline()
		except:
			pass
		direct = str(QtGui.QFileDialog.getExistingDirectory(self, "Select Directory", last_dir))
		if direct != "":
			self._open_model(direct)

	def _open_model(self, direct):
		#if self.is_valid_input_directory(direct):
		self.inputpathlabel.setText(direct)
		self.populate_model_data(direct)
		with open(".last_dir", "w") as fo:
			fo.write(direct)

	def choose_output_dir(self):
		last_dir = "."
		try:
			with open(".last_outdir") as fi:
				last_dir = fi.readline()
		except:
			pass
		direct = str(QtGui.QFileDialog.getExistingDirectory(self, "Select Directory", last_dir))		
		if direct != "":
			if self.is_valid_input_directory(direct):
				self.outputpathlabel.setText(direct)
				#self.populate_model_data(direct)
				with open(".last_outdir", "w") as fo:
					fo.write(direct)
				self.output_valid = True
		else:
			self.output_valid = False

		if self.is_everything_ready(): self.enable_simulation()

	def simulate(self):
		# step 0: create files
		t_vector = linspace(0, float(self.timemax.value()), int(self.samples.value()))
		in_dir  = str(self.inputpathlabel.text())
		out_dir = str(self.outputpathlabel.text())
		savetxt(in_dir+"/t_vector", t_vector)

		cs_vector = []
		for x in xrange(self.speciestable.rowCount()):
			if self.speciestable.item(x,2).checkState() == QtCore.Qt.Checked:
				cs_vector.append(x)
		savetxt(in_dir+"/cs_vector", cs_vector, fmt="%d")

		with open(in_dir+"/modelkind", "w") as fo:
			fo.write("deterministic")
		
		# step 1: launch simulation (platform dependent)
		if "Windows" in platform():
			print " * Windows detected, using specific binary"
			binary = "lassieWin.exe"			

		command = [binary, "-double", in_dir, out_dir]
		print " ".join(command)
		ret = call(command)

		# step 2: plot results
		results = loadtxt(out_dir+"/output/Solution")
		self.figure.clf()
		ax = self.figure.add_subplot(1,1,1)
		y = 0
		for x in xrange(self.speciestable.rowCount()):
			if self.speciestable.item(x,2).checkState() == QtCore.Qt.Checked:
				ax.plot(results.T[0], results.T[y+1], label=str(self.speciestable.item(x,0).text()))
				y += 1
		ax.set_xlabel("Time")
		ax.set_ylabel("Amount")
		self.figure.tight_layout()
		ax.legend()
		self.figure.show()

if __name__ == '__main__':
	

	app    = QtGui.QApplication(sys.argv)
	window = MyWindow()
	sys.exit(app.exec_())
