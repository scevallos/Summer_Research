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

def calc_shifts(shielding_vals, TMS):
	shifts = []
	for s in shielding_vals:
		shifts.append(-(s-TMS))
	return shifts

def calc_delGs(Gs):
	delGs = []
	for x in Gs:
		delGs.append(x - Gs[0])
	return delGs

def calc_Njs_N(delGs, T):
	e = 2.71828182845904
	Njs = [1] # Nj_d20 is 1, start with it inside the list
	for x in delGs[1:]:
		Njs.append(e**(-x/T/0.0019872))
	N = sum(Njs)
	return Njs, N

def calc_pops(Njs, N):
	pops = []
	for x in Njs:
		pops.append(x/N)
	return pops

def calc_shift_diff(pops, shifts):
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

	return H_r, H_s, shift_diff

def calc_nmr_shifts(H_vals, S_vals, T, shielding_vals, TMS):
	Gs = []
	for h,s in zip(H_vals, S_vals):
		Gs.append(calc_G(h, s, T))
	delGs = calc_delGs(Gs)
	Njs,N = calc_Njs_N(delGs, T)
	pops = calc_pops(Njs, N)
	H_r, H_s, shift_diff = calc_shift_diff(pops, calc_shifts(shielding_vals, TMS))

	return H_r, H_s, shift_diff

def manual_shift_diff():
	print "Please enter the temperature (K)."
	T = float(raw_input()) # 298.15

	print "Please enter the H values, separated by a space.."
	hvals = raw_input().split() #not sure about this

	print "Please enter the S values, separated by a space.."
	svals = raw_input().split() #not sure about this

	print "Calculating NMR Shifts..."

	r, s, d = calc_nmr_shifts([130.8150, 130.8930, 130.8290], [88.8790, 88.9100, 88.8840], T, [30.0125, 30.9039, 30.5504], 32.426)

	# print_table(G_ds, delGs, Njs, N, pops, shifts, H_r, H_s, shift_diff)
	print "--------------------------------------------------"
	print "H_r", "\t", r, "ppm"
	print "H_s", "\t", s, "ppm"
	print "H_r - H-s =", d, " ppm"
	print "--------------------------------------------------"

def auto_shift_diff():
	Tk().withdraw()
	
	print "Please enter the atom numbers you would like to get the shifts of, seperated by a space."
	atom_nums = raw_input().split()
	print "Please select the Gaussian log file containing the NMR calculations"
	fname = askopenfilename()

	TMS = 32.426

	shielding = get_shielding(atom_nums, fname)
	shielding = map(float, shielding)
	shifts_IO = map(str, calc_shifts(shielding, TMS))

	print "Shift (ppm)"
	print "--------------------------------------------------"
	print "Atom No.", "\t", "\t".join(atom_nums)
	print "Shifts", "\t \t", "\t".join(shifts_IO)
	print "--------------------------------------------------"

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
	
	# manual_shift_diff()
	auto_shift_diff()

if __name__ == '__main__':
	main()