from scipy.special import erfc
import numpy as np

from cmath import pi
from cmath import exp
import cmath
import math

from ase import *
from ase.visualize import view
from ase.geometry import Cell

def energy_function(r,pos):
	# Calculate electrostatic short range
	esum = 0
	# shifts = self.get_shifts(self.real_cut_off, atoms.get_cell())
	count = 1
	for ioni in range(0,len(r),3):
		for ionj in range(0,len(r),3):
			if ioni != ionj:  # skip in case it's the same ion in original unit cell
				print(count)
				rij = np.array([r[ioni]-r[ionj], r[ioni+1]-r[ionj+1], \
													r[ioni+2]-r[ionj+2]])
				dist = np.linalg.norm(rij)
				esum += (self.get_charges_mult(ioni, ionj) *
									 math.erfc(self.alpha*dist)/(2*dist))
			# take care of the rest lattice (+ Ln)
			for shift in shifts:
				rij = np.array([r[ioni]+shift[0]-r[ionj], \
								r[ioni+1]+shift[1]-r[ionj+1], \
								r[ioni+2]+shift[2]-r[ionj+2]])
				dist = np.linalg.norm(rij)
				esum += (self.get_charges_mult(ioni, ionj) *
									 math.erfc(self.alpha*dist)/(2*dist))
