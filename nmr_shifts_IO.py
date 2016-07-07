#!/usr/bin/python
from itertools import cycle
from Tkinter import Tk
from tkFileDialog import askopenfilename

def get_shielding(atom_numbers, log_file_name):
	shielding = []

	with open(log_file_name, 'r') as f:
		lines = f.readlines()
		for i in range(0, len(lines)):
			line = lines[i]
			for i in range(0,len(atom_numbers)):
				if str(atom_numbers[i]) + "  " + "H    I"in line:
					shielding.append(line.split()[4])
	return shielding

# ZPVE not used for any calculations?
def calc_G(H, S, T):
	return H-(T*(S/1000))

def calc_Nj(delG, T):
	e = 2.71828182845904
	return e**(-delG/T/0.0019872)

def calc_shift(shielding, TMS):
	return -(shielding - TMS)

def print_table(G_ds, delGs, Njs, N, pops, shifts, Hr, Hs, shift_diff):
	print "--------------------------------------------------"
	print "\t", "D = 20", "\t \t", "D = 21", "\t \t", "D = 22"
	print "G", "\t", G_ds[0], "\t", G_ds[1], "\t", G_ds[2]
	print "delG", "\t", delGs[0], "\t \t", delGs[1], "\t", delGs[2]
	print "Nj", "\t", Njs[0], "\t \t", Njs[1], "\t", Njs[2]
	print "N", "\t", N
	print "Pops", "\t", pops[0], "\t", pops[1], "\t", pops[2]
	print "--------------------------------------------------"
	print "\t", "H = 20", "\t \t", "H = 21", "\t \t", "H = 22"
	print "shift", "\t", shifts[0], "\t \t", shifts[1], "\t \t", shifts[2]
	print "H_r", "\t", Hr, "ppm"
	print "H_s", "\t", Hs, "ppm"
	print "H_r - H-s =", shift_diff, " ppm"
	print "--------------------------------------------------"

def main():
	Tk().withdraw()

	T = 298.15

	G_d20 = calc_G(130.8150, 88.8790, T)
	G_d21 = calc_G(130.8930, 88.9100, T)
	G_d22 = calc_G(130.8290, 88.8840, T)

	G_ds = [G_d20, G_d21, G_d22]

	delGs = []
	for x in G_ds:
		delGs.append(x - G_ds[0])

	Njs = [1] # Nj_d20 is 1, start with it inside the list
	for x in delGs[1:]:
		Njs.append(calc_Nj(x, T))

	N = sum(Njs)

	pops = []
	for x in Njs:
		pops.append(x/N)

	TMS = 32.426

	shifts = []
	shifts.append(calc_shift(30.0125, TMS))
	shifts.append(calc_shift(30.9039, TMS))
	shifts.append(calc_shift(30.5504, TMS))

	shifts_cycle = cycle(shifts)
	next(shifts_cycle) # start one index up (shift_h21)

	Rvals = []
	Svals = []

	# flip pop numbers of d22 & d21; needed to do this for iteration below
	pops[1], pops[2] = pops[2], pops[1]
	
	for x in pops:
		Rvals.append(x * next(shifts_cycle))
		Svals.append(x * next(shifts_cycle))

	H_r = sum(Rvals)
	H_s = sum(Svals)
	shift_diff = H_r - H_s

	# print_table(G_ds, delGs, Njs, N, pops, shifts, H_r, H_s, shift_diff)
	print "Please enter the atom numbers you would like to get the shielding of, seperated by a space."
	atom_nums = raw_input().split()
	print "Please select the Gaussian log file containing the NMR calculations"
	fname = askopenfilename()

	shielding = get_shielding(atom_nums, fname)
	shielding = map(float, shielding)
	a = [calc_shift(x, TMS) for x in shielding]
	


if __name__ == '__main__':
	main()