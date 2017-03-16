#from gev import GEV
from PyQt4 import QtGui, QtCore, uic
from numpy import *
import sys
from subprocess import check_output, Popen, call
from platform import platform
from pylab import *
try:
	from SBML2BSW import SBMLloader
	use_sbml = True
except:
	print "WARNING: libsbml not found"
	use_sbml = False

import resources_rc

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
matplotlib.rcParams.update({'font.size': 8})


class MyDialog(QtGui.QDialog):
	def __init__(self):
		super(MyDialog, self).__init__()
		uic.loadUi('about-lassie.ui', self)


class OptionsDialog(QtGui.QDialog):
	def __init__(self, mainref=None):
		super(OptionsDialog, self).__init__()
		uic.loadUi('options.ui', self)
		self._mainref=mainref

	def select_file(self):
		file = str(QtGui.QFileDialog.getOpenFileName(self, "Select Directory"))
		if file!="":
			self.binary_path.setText(file)

	def accept(self):
		self._mainref._binary = self.binary_path.text()
		print "Settando:", self._mainref._binary
		self.close()


class GraphCanvas(FigureCanvas):

    def __init__(self, parent=None, dpi=100):

        self.fig = Figure(dpi=dpi)                                
        FigureCanvas.__init__(self, self.fig)
        FigureCanvas.setSizePolicy(self,
                                   QtGui.QSizePolicy.Expanding,
                                   QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def drawGraph(self, title="", data={}, speciestable=None, results=[]):
		if len(results)==0: return
		self.axes = self.fig.add_subplot(111)        
		self.axes.cla()         					
		for x in xrange(speciestable.rowCount()):
			if speciestable.item(x,2).checkState() == QtCore.Qt.Checked:
				self.axes.plot(results.T[0], results.T[x+1], "o-", markersize=3, label=str(speciestable.item(x,0).text()))
				print "Plotting", str(speciestable.item(x,0).text())
		self.axes.set_xlabel("Time")
		self.axes.set_ylabel("Amount")		
		self.axes.legend()
		self.fig.tight_layout()
		self.draw()
		self.background = self.fig.canvas.copy_from_bbox(self.axes.bbox)

    def onMove(self, event):
        

        # cursor moves on the canvas
        if event.inaxes:

            # restore the clean background
            self.fig.canvas.restore_region(self.background)
            ymin, ymax = self.axes.get_ylim()
            x = event.xdata - 1

            # draw each vertical line
            for line in self.verticalLines:
                line.set_xdata((x,))
                line.set_ydata((ymin, ymax))
                self.axes.draw_artist(line)
                x += 1

            self.fig.canvas.blit(self.axes.bbox)



class MyWindow(QtGui.QMainWindow):

	def __init__(self):
		super(MyWindow, self).__init__()
		uic.loadUi('lassie-gui.ui', self)

		self.input_valid = False
		self.output_valid = False

		self.show()
		self.figure = figure()
		self._options = OptionsDialog(mainref=self)
		self._options.setModal(True)
		self._about = MyDialog()
		self._about.setModal(True)
		self._results = []

		self._plot_canvas = GraphCanvas()
		self.plotlayout.layout().addWidget(self._plot_canvas)
		self.statusBar().showMessage("LASSIE ready, please open a model")
		#self.populate_model_data("C:\\Users\\aresio\\Documents\\PythonCode\\ModelliTest\\MM")


		# step 1: launch simulation (platform dependent)
		if "Windows" in platform():
			print " * Windows detected, using specific binary"
			self._binary = "lassieWin.exe"
		elif "Linux" in platform():
			print " * Linux detected, using specific binary"
			self._binary = "./lassie"
		else:
			print " * OSX detected, using specific binary"
			self._binary = "./lassie"

		try:
			with open(".binary") as fi:
				self._binary = fi.readline()
		except:
			pass		

	def save_files(self):
		file = str(QtGui.QFileDialog.getSaveFileName(self, "Select output file"))
		if file!="":
			savetxt(file, self._results)

	def show_help(self):
		self._about.show()

	def show_options(self):
		self._options.binary_path.setText(self._binary)
		self._options.show()
		print "Settato:", self._binary
		with open(".binary", "w") as fo:
			fo.write(self._binary)			

	def import_SBML(self):
		if not use_sbml:
			print "libsbml not found, aborting"
			return
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
		#print species_names

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
			
		self.speciestable.resizeColumnsToContents()
		self.speciestable.itemClicked.connect(self.handle_render)
		
	def handle_render(self, item):
		row, column = item.row(), item.column()
		if column==2:	# column of species to show
			self._plot_canvas.drawGraph(title="Bla", speciestable=self.speciestable, results=self._results)

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
		self.reactionstable.resizeColumnsToContents()

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
		self.inputpathlabel.setText(direct)
		self.populate_model_data(direct)
		if self.input_valid:
			print " * Model loaded"
			with open(".last_dir", "w") as fo:
				fo.write(direct)			
			self.input_greenlight.setEnabled(True)			
			self.statusBar().showMessage( "Species: "+str(self.speciestable.rowCount())+", Reactions: "+str(self.reactionstable.rowCount()) )

		else:
			print " * Model not loaded"
			#self.inputpathlabel.setStyleSheet('border: 3px solid red')
			self.input_greenlight.setEnabled(False)

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
				#self.outputpathlabel.setStyleSheet('border: 3px solid lime')
				self.output_greenlight.setEnabled(True)
		else:
			self.output_valid = False
			#self.outputpathlabel.setStyleSheet('border: 3px solid red')
			self.output_greenlight.setEnabled(False)

		if self.is_everything_ready(): self.enable_simulation()

	def simulate(self):
		# step 0: create files

		try:
			close('all')
		except:
			print "Not figure open"

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

		command = [str(self._binary), "-double", in_dir, out_dir]
		print " ".join(command)
		ret = call(command)

		# step 2: plot results
		self._results = loadtxt(out_dir+"/output/Solution")
		self._plot_canvas.drawGraph(title="Bla", speciestable=self.speciestable, results=self._results)
		self.save_to_file.setEnabled(True)

if __name__ == '__main__':


	app    = QtGui.QApplication(sys.argv)
	window = MyWindow()
	sys.exit(app.exec_())
