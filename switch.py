from scipy.optimize import minimize
from ase.io import read as aread
import numpy as np
import sys

sys.path.append('../Switch_Implementation/gradients_implementation')
import potential as pot
from energy_function import *


charge_dict = {
	'O' : -2.,
	'Sr':  2.,
	'Ti':  4.,
	'Cl': -1.,
	'Na':  1.,
	'S' : -2.,
	'Zn':  2.
}


def flatten(positions):
	vector = []
	for sublist in positions:
		for element in sublist:
			vector.append(element)
	return vector

def energy_function_wrapper(r, N, *args):
	assert r.shape == (3*N,)
	pos = r.reshape(15,3)

	Cpot = pot.Coulomb()
	Cpot.set_parameters(*args[1:6])

	Bpot = pot.Buckingham()
	Bpot.set_parameters(*args[6:])
	energy = buckingham_coulomb(pos, Cpot, Bpot, N, *args[0])
	print(energy)

if __name__ == '__main__':
	atoms = aread("../../Data/random/RandomStart_Sr3Ti3O9/1.cif")
	r = np.array(flatten(atoms.positions))
	vects = np.array(atoms.get_cell())
	volume = abs(np.linalg.det(vects))
	N = len(atoms.get_positions())
	accuracy = 0.00001  # Demanded accuracy of terms 
	alpha = N**(1/6) * pi**(1/2) / volume**(1/3)
	real_cut = (-np.log(accuracy))**(1/2)/alpha
	recip_cut = 2*alpha*(-np.log(accuracy))**(1/2)


	energy_function_wrapper(r, len(atoms.positions), alpha, real_cut,\
								recip_cut, atoms.get_chemical_symbols(),\
								charge_dict)
	# res = minimize(funct, atoms.positions, method='BFGS', \
	# 	               options={'gtol': 1e-3, 'disp': True})