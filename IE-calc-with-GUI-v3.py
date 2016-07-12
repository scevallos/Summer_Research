#!/usr/bin/python

import numpy as np
import math as mt
import sys as sys
import os
import subprocess as sub
from operator import itemgetter
import csv
import datetime
import tkMessageBox
import matplotlib.pyplot as pyplot
from Tkinter import *
from ttk import *
import itertools
import tkFileDialog

def find_log_files(path):
	os.chdir(path)
	log_files = []
	for files in os.listdir('.'):
		if files.endswith("log"):
			log_files.append(str(files))
	return log_files

def create_KIE_tab(notebook):
	# Sets up frame
	t1 = Frame(notebook, name='kie')
	t1.grid()

	# Blank1 = Label(t1).grid(row=0)
	# Blank4 = Label(t1).grid(row=9)

	# Select directory to get log files from
	# dirname = tkFileDialog.askdirectory()

	# Find log files in chosen dir
	# options = find_log_files(dirname)

	#Preparing drop-down menu variables
	variable1 = StringVar(t1)
	variable1.set("")
	variable2 = StringVar(t1)
	variable2.set("")
	variable3 = StringVar(t1)
	variable3.set("")
	variable4 = StringVar(t1)
	variable4.set("")
	options = ["1", "2", "3", "4"]


	L8 = Label(t1, text = "File names:").grid(row=10)
	L9 = Label(t1, text = "Ground state, unlabeled:").grid(row=11)
	GS_unlabeled = apply(OptionMenu, (t1, variable1) + tuple(options))
	GS_unlabeled.grid(row=11, column=1)
	L10 = Label(t1, text = "Ground state, labeled:").grid(row=12)
	GS_labeled = apply(OptionMenu, (t1, variable2) + tuple(options))
	GS_labeled.grid(row=12, column=1)
	L11 = Label(t1, text = "Transition state, unlabeled:").grid(row=13)
	TS_unlabeled = apply(OptionMenu, (t1, variable3) + tuple(options))
	TS_unlabeled.grid(row=13, column=1)
	L12 = Label(t1, text = "Transition state, labeled:").grid(row=14)
	TS_labeled = apply(OptionMenu, (t1, variable4) + tuple(options))
	TS_labeled.grid(row=14, column=1)

	Blank6 = Label(t1).grid(row=15)
	Blank7 = Label(t1).grid(row=16)

	B = Button(t1, text ="Calculate Isotope Effect").grid(columnspan=2)

	Blank8 = Label(t1).grid(row=18)

	return t1

def create_EIE_tab(notebook):
	# Sets up frame
	t2 = Frame(notebook, name="eie")
	t2.grid()

	# EIE GUI goes here

	return t2

def create_NMR_tab(notebook):
	# Sets up frame
	t3 = Frame(notebook, name="nmr")
	t3.grid()

	la_nums = Label(t3, text="Atom Numbers:")
	la_nums.grid()
	a_nums = Entry(t3)
	a_nums.grid()

	nmrb = Button(t3, text="Compute NMR Shifts")
	nmrb.grid()

	return t3

def create_top_frame(pane):
	# Set up top panel (stuff not changing via tabs)
	topf = Frame(pane)
	topf.grid()

	# copy-pasted Alex's code of GUI setup
	Blank = Label(topf).grid(row=0)

	Label(topf, text = "Min. Temp. (K):").grid(row=1)
	min_temp = Entry(topf, width = 10)
	min_temp.grid(row=1, column=1)
	min_temp.insert(0, "298.15")

	Label(topf, text = "Max. Temp. (K):").grid(row=2)
	max_temp = Entry(topf, width = 10)
	max_temp.grid(row=2, column=1)
	max_temp.insert(0, "298.15")

	Label(topf, text = "Step between temperatures:").grid(row=3)
	step_size = Entry(topf, width=10)
	step_size.grid(row=3, column=1)
	step_size.insert(0, "1.0")

	# Setup graph option
	var = IntVar()
	var.set(2)

	# Diplay graph yes, no buttons
	Label(topf, text = "Display graph (temp vs IE?)").grid(row=4)
	Radiobutton(topf, text="Yes", variable=var, value=1).grid(row=4, column=1)
	Radiobutton(topf, text="No", variable=var, value=2).grid(row=4, column=2)

	Blank2 = Label(topf).grid(row=5)

	# Scaling factor entry
	Label(topf, text = "Scaling factor:").grid(row=8)
	scale_factor = Entry(topf, width = 10)
	scale_factor.grid(row=8, column=1)
	scale_factor.insert(0, "1.0")

	return topf

def create_bottom_frame(pane):
	# Setup bottom panel (Notebook, tab changing)
	botf = Frame(pane)
	botf.grid()

	# Notebook widget allows for tabbed functionality
	nb = Notebook(botf, name="nb")
	nb.grid()

	# Creating 3 tabs via respective methods & adding them to notebook
	t1 = create_KIE_tab(nb)
	nb.add(t1, text="KIE")

	t2 = create_EIE_tab(nb)
	nb.add(t2, text="EIE")

	t3 = create_NMR_tab(nb)
	nb.add(t3, text="NMR")

	return botf
	

def GUI():

	# Setup main window & set title
	root = Tk()
	root.title("IE Calculations")

	# Create main frame, holds PanedWindow, holding everything
	mf = Frame(root, name="mainframe")
	mf.grid()

	# Create pane, split into top & bottom panels
	pane = PanedWindow(root, orient=VERTICAL)
	pane.grid()

	# Create top frame
	top_frame = create_top_frame(pane)

	# Adding top frame to the pane
	pane.add(top_frame)

	# Create bottom frame
	bottom_frame = create_bottom_frame(pane)

	# Adding bottom frame to the pane
	pane.add(bottom_frame)

	# mainloop
	root.mainloop()

# =====================================================================================
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

	# print "Please enter the H values, separated by a space.."
	# hvals = raw_input().split() #not sure about this

	# print "Please enter the S values, separated by a space.."
	# svals = raw_input().split() #not sure about this

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
# =====================================================================================

def main():
	GUI()

if __name__ == '__main__':
	main()