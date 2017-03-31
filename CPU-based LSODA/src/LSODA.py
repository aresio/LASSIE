import time
import sys
import math
import scipy
import collections
from scipy import *
from scipy import linspace
from scipy.integrate import ode
from scipy.integrate import odeint
import numpy as np


class LSODER():

	def __init__(self):
		self.PATH = ""
		self.t0 = 0
		self.points = 0
		self.parametri = 0
		self.time_max = 0
		self.time_instants = []
		self.max_steps = 0
		self.atolvector = None
		self.rtol = None


	def load_data(self, pa):
		self.PATH = pa
		exec open(self.PATH).read() in globals()
		self.parametri = parametri
		self.time_max = tEnd

	def load_data_nofile(self, stringa):
		exec stringa in globals()
		self.parametri = parametri
		self.time_max = tEnd


	def run(self, subspecies):

		ret1 = odeint(ODE, y0, self.time_instants , Dfun=JAC, args=(self.parametri,), mxhnil=0, printmessg=False, mxstep = self.max_steps, atol=self.atolvector, rtol=self.rtol )
			
		return ret1