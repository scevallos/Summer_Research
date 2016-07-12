#!/usr/bin/python

import numpy as np
import math as mt
import sys as sys
import os
import subprocess as sub
import itertools
from operator import itemgetter
import csv
import datetime
from Tkinter import *
import tkMessageBox
import matplotlib.pyplot as pyplot

def GUI():
    window = Tk()
    Blank1 = Label(window).grid(row=0)

    L1 = Label(window, text = "Min. Temp. (K):").grid(row=1)
    min_temp = Entry(window, width = 10)
    min_temp.grid(row=1, column=1)
    min_temp.insert(0, "298.15")

    L2 = Label(window, text = "Max. Temp. (K):").grid(row=2)
    max_temp = Entry(window, width = 10)
    max_temp.grid(row=2, column=1)
    max_temp.insert(0, "298.15")

    L3 = Label(window, text = "Step between temperatures:").grid(row=3)
    step_size = Entry(window, width=10)
    step_size.grid(row=3, column=1)
    step_size.insert(0, "1.0")

    var= IntVar()
    var.set(2)

    L4 = Label(window, text = "Display graph (temp vs IE?)").grid(row=4)
    Option1 = Radiobutton(window, text="Yes", variable=var, value=1)
    Option1.grid(row=4, column=1)
    Option2 = Radiobutton(window, text="No", variable=var, value=2)
    Option2.grid(row=4, column=2)

    Blank2 = Label(window).grid(row=5)

    L6 = Label(window, text = "Scaling factor:").grid(row=8)
    scale_factor = Entry(window, width = 10)
    scale_factor.grid(row=8, column=1)
    scale_factor.insert(0, "1.0")

    Blank4 = Label(window).grid(row=9)

    #Preparing drop-down menu variables
    variable1 = StringVar(window)
    variable1.set("")
    variable2 = StringVar(window)
    variable2.set("")
    variable3 = StringVar(window)
    variable3.set("")
    variable4 = StringVar(window)
    variable4.set("")
    options = find_log_files()

    if len(options) == 0:
        print "No log files found in root directory!"
        quit()

    L8 = Label(window, text = "File names:").grid(row=10)
    L9 = Label(window, text = "Ground state, unlabeled:").grid(row=11)
    GS_unlabeled = apply(OptionMenu, (window, variable1) + tuple(options))
    GS_unlabeled.grid(row=11, column=1)
    L10 = Label(window, text = "Ground state, labeled:").grid(row=12)
    GS_labeled = apply(OptionMenu, (window, variable2) + tuple(options))
    GS_labeled.grid(row=12, column=1)
    L11 = Label(window, text = "Transition state, unlabeled:").grid(row=13)
    TS_unlabeled = apply(OptionMenu, (window, variable3) + tuple(options))
    TS_unlabeled.grid(row=13, column=1)
    L12 = Label(window, text = "Transition state, labeled:").grid(row=14)
    TS_labeled = apply(OptionMenu, (window, variable4) + tuple(options))
    TS_labeled.grid(row=14, column=1)

    Blank6 = Label(window).grid(row=15)
    Blank7 = Label(window).grid(row=16)

    B = Button(window, text ="Calculate Isotope Effect", command = lambda: calculate_IE(variable1.get(), variable2.get(), variable3.get(), variable4.get(), max_temp.get(), min_temp.get(), scale_factor.get(), step_size.get(), var.get())).grid(columnspan=2)

    Blank8 = Label(window).grid(row=18)

    window.mainloop()

def find_log_files():
    log_files = []
    for files in os.listdir('.'):
        if files.endswith("log"):
            log_files.append(str(files))
    return log_files

def build_dictionary():
    """
        Input: none
        Output: dictionary with all elements that have multiple stable isotopes. Dictionary has only the most abundant isotope
        Used to determine the KIE in the correct order (light/heavy)
    """
    elements = {'H':1, 'He':4, 'Li':7, 'B':11, 'C':12, 'N':14, 'O':16, 'Ne':20, 'Mg':24, 'Si':28, 'S':32, 'Cl':35, 'Ar':40, 'K':39, 'Ca':40, 'Ti':48, 'Cr':52, 'Fe':56, 'Ni':58, 'Cu':63, 'Zn':64, 'Ga':69, 'Ge':74, 'Se':80, 'Br':79, 'Kr':84, 'Sr':88, 'Zr':90, 'Mo':98, 'Ru':102, 'Pd':106, 'Ag':107, 'Cd':114, 'Sn':120, 'Sb':121, 'Te':126, 'Xe':132, 'Ba':138, 'Ce':140, 'Nd':142, 'Sm':152, 'Eu':153, 'Gd':158, 'Dy':164, 'Er':166, 'Yb':174, 'Hf':180, 'W':184, 'Os':192, 'Ir':193, 'Pt':195, 'Hg':202, 'Tl':205, 'Pb':208}
    return elements

def try_int(x):
    """
        Input: list with strings
        Output: list with integers and strings
        If the index can be converted to an integer, it will be. Otherwise it will stay a string
    """
    try:
        return int(x)
    except ValueError:
        return x

def no_input():
    """
        Input: None. Def called only if txt file is not given in command line
        Output: Text file for user to fill out
        Populates a text file that the user then fills out to enter with the code
    """
    """ Input population for checkpoint files
    ERROR_OUT = open("INPUT_FCHK.txt","w")
    ERROR_OUT.write("Please fill out the following text file. Words in brackets are the standard response to the question" + "\n")
    ERROR_OUT.write("Ground state file (can be multiple):" + "\n\n")
    ERROR_OUT.write("Transition state file (can be multiple):" + "\n\n")
    ERROR_OUT.write("Do you want to write hyperchem files? [No]" + "\n\n")
    ERROR_OUT.write("Temperature (in Kelvin) (can be multiple):" + "\n\n")
    ERROR_OUT.write("Pressure (in atm) (only one):" + "\n\n")
    ERROR_OUT.write("Scaling factor for frequencies during thermochemistry: (only one)" + "\n\n")
    ERROR_OUT.write("Project out gradient direction [No]:" + "\n\n")
    ERROR_OUT.write("Atom symbol - ground state atom number - transition state atom number - isotope" + "\n")
    ERROR_OUT.write("Put at least one space between each input. DO NOT use any other characters to signify spaces!" + "\n\n")
    ERROR_OUT.close()
    """

    filename = datetime.datetime.now().strftime("%m%d%Y-%H%M%S")
    filename = "INPUT-", filename, ".txt"
    filename = ''.join(filename) 

    ERROR_OUT = open(filename,"w")
    ERROR_OUT.write("Please fill out the following text file." + "\n")
    ERROR_OUT.write("Unlabeled transition state file:" +"\n\n")
    ERROR_OUT.write("Labeled transition state file:" + "\n\n")
    ERROR_OUT.write("Unlabeled ground state file:" + "\n\n")
    ERROR_OUT.write("Labeled ground state file:" + "\n\n")
    ERROR_OUT.write("Temperature(in Kelvin):" + "\n\n")
    ERROR_OUT.write("Scaling factor:" + "\n\n")
    ERROR_OUT.close()

def freqchk_method():
    """
    Using a checkpoint file to calculate IE
    """
    y = 0
    csvfile = open('KIE_Output.csv', 'wb')
    outputwriter = csv.writer(csvfile, dialect='excel')
    outputwriter.writerow(["GS_file"] + ["TS_File"] + ["Temperature"] + ["Pressure"] + ["Scale_factor"] + ["Chem_symbol"] + ["GS_number"] + ["TS_number"] + ["Labeled_isotope"] + ["KIE"] + ["KIE_tunneling"])
    heading = ["GS_file", "TS_File", "Temperature", "Pressure", "Scale_factor", "Chem_symbol", "GS_number", "TS_number", "Labeled_isotope", "KIE", "KIE_tunneling"]
    x=0
    for each_GS in GS_chkpt_file:
        for each_TS in TS_chkpt_file:
            for each_item in isotope_changes: #split the changes to the isotope
                for each_temp in temp:
                    if len(each_GS) > x:
                        x = len(each_GS)
                    if len(each_TS) > x:
                        x = len(each_TS)
                        y += 1
                    chem_sym = []
                    gs_num = []
                    ts_num = []
                    isotope_mass = []
                    output = []
                    i = 0
                    while i < len(each_item):
                        chem_sym.append(each_item[i])
                        gs_num.append(each_item[i+1])
                        ts_num.append(each_item[i+2])
                        isotope_mass.append(each_item[i+3])
                        i += 4
                
                    # run freqchk for TS without a marker
                    run_freqchk_TS_no_marker(each_TS, hyperchem_files, each_temp, pressure, scale_factor, gradient_direction)
                    #run freqchk for TS with a marker
                    run_freqchk_TS_marker(each_TS, hyperchem_files, each_temp, pressure, scale_factor, gradient_direction, ts_num, isotope_mass, number_atoms_TS)
                    #run freqchk for GS without a marker
                    run_freqchk_GS_no_marker(each_GS, hyperchem_files, each_temp, pressure, scale_factor, gradient_direction)
                    #run freqchk for GS with a marker
                    run_freqchk_GS_marker(each_GS, hyperchem_files, each_temp, pressure, scale_factor, gradient_direction, gs_num, isotope_mass, number_atoms_GS)

                    #get frequencies from .txt file
                    frequency_TS_natural = map(float, get_frequencies("freq_TS_no_marker.txt"))
                    frequency_TS_isotope = map(float, get_frequencies("freq_TS_marker.txt"))
                    frequency_GS_natural = map(float, get_frequencies("freq_GS_no_marker.txt"))
                    frequency_GS_isotope = map(float, get_frequencies("freq_GS_marker.txt"))

                    #remove freqchk output files
                    os.system("rm freq_TS_no_marker.txt freq_TS_marker.txt freq_GS_no_marker.txt freq_GS_marker.txt")

                    # create array with u values
                    u_TS_natural = np.array(calc_u(frequency_TS_natural, each_temp, scale_factor))
                    u_TS_isotope = np.array(calc_u(frequency_TS_isotope, each_temp, scale_factor))
                    u_GS_natural = np.array(calc_u(frequency_GS_natural, each_temp, scale_factor))
                    u_GS_isotope = np.array(calc_u(frequency_GS_isotope, each_temp, scale_factor))
                    if u_TS_natural[0] < 0:
                        u_neg_TS_natural= u_TS_natural[0] #negative u value
                        u_TS_natural = u_TS_natural[1:] #allows calculation on all u values that are positive
                        u_neg_TS_isotope = u_TS_isotope[0] #negative u value
                        u_TS_isotope = u_TS_isotope[1:] #allows calculation on all u values that are positive
                    else:
                        u_neg_TS_natural =[]
                        u_neg_TS_isotope = []

                    # create array with exp(u/2) values
                    exp_TS_natural = np.array(exp_u_half(u_TS_natural))
                    exp_TS_isotope = np.array(exp_u_half(u_TS_isotope))
                    exp_GS_natural = np.array(exp_u_half(u_GS_natural))
                    exp_GS_isotope = np.array(exp_u_half(u_GS_isotope))

                    # create array with 1-exp(-u) values
                    one_minus_exp_TS_natural = np.array(calc_one_minus_exp(u_TS_natural))
                    one_minus_exp_TS_isotope = np.array(calc_one_minus_exp(u_TS_isotope))
                    one_minus_exp_GS_natural = np.array(calc_one_minus_exp(u_GS_natural))
                    one_minus_exp_GS_isotope = np.array(calc_one_minus_exp(u_GS_isotope))

                    # create array with prod values
                    prod_TS_natural = np.array(calc_prod(u_TS_natural, exp_TS_natural, one_minus_exp_TS_natural))
                    prod_TS_isotope = np.array(calc_prod(u_TS_isotope, exp_TS_isotope, one_minus_exp_TS_isotope))
                    prod_GS_natural = np.array(calc_prod(u_GS_natural, exp_GS_natural, one_minus_exp_GS_natural))
                    prod_GS_isotope = np.array(calc_prod(u_GS_isotope, exp_GS_isotope, one_minus_exp_GS_isotope))
               
                    # calculate FTS
                    if u_neg_TS_natural:
                        FTS_TS_natural = calc_FTS_TS(prod_TS_natural, u_neg_TS_natural)
                        FTS_TS_isotope = calc_FTS_TS(prod_TS_isotope, u_neg_TS_isotope)
                        FTS_GS_natural = calc_FTS(prod_GS_natural)
                        FTS_GS_isotope = calc_FTS(prod_GS_isotope)
                    else:
                        FTS_TS_natural = calc_FTS(prod_TS_natural)
                        FTS_TS_isotope = calc_FTS(prod_TS_isotope)
                        FTS_GS_natural = calc_FTS(prod_GS_natural)
                        FTS_GS_isotope = calc_FTS(prod_GS_isotope)
                
                    # calcualte qt for TS
                    if u_neg_TS_natural:
                        qt_TS_natural = calc_qt(u_neg_TS_natural)
                        qt_TS_isotope = calc_qt(u_neg_TS_isotope)
                    else:
                        qt_TS_natural = calc_qt(u_TS_natural[0])
                        qt_TS_isotope = calc_qt(u_TS_isotope[0])

                    # build dictionary with elements and get the mass of the
                    # elements being used
                    elements = {'H':1, 'He':4, 'Li':7, 'B':11, 'C':12, 'N':14, 'O':16, 'Ne':20, 'Mg':24, 'Si':28, 'S':32, 'Cl':35, 'Ar':40, 'K':39, 'Ca':40, 'Ti':48, 'Cr':52, 'Fe':56, 'Ni':58, 'Cu':63, 'Zn':64, 'Ga':69, 'Ge':74, 'Se':80, 'Br':79, 'Kr':84, 'Sr':88, 'Zr':90, 'Mo':98, 'Ru':102, 'Pd':106, 'Ag':107, 'Cd':114, 'Sn':120, 'Sb':121, 'Te':126, 'Xe':132, 'Ba':138, 'Ce':140, 'Nd':142, 'Sm':152, 'Eu':153, 'Gd':158, 'Dy':164, 'Er':166, 'Yb':174, 'Hf':180, 'W':184, 'Os':192, 'Ir':193, 'Pt':195, 'Hg':202, 'Tl':205, 'Pb':208}
                    temp_sym = ' '.join(chem_sym[0])
                    temp_isotope_mass = isotope_mass[0]
                    temp_isotope_mass = "".join(repr(temp_isotope_mass))
                    temp_isotope_mass = int(temp_isotope_mass)

                    # calculate KIE
                    for a in elements.keys():
                        if a == ''.join(temp_sym):
                            if elements[a] > temp_isotope_mass:
                                KIE1 = (FTS_TS_isotope / FTS_TS_natural)
                                KIE2 = (FTS_GS_natural / FTS_GS_isotope)
                                KIE3 = KIE2 * KIE1
                                KIE = 1 / KIE3
                            else:
                                KIE1 = (FTS_TS_isotope / FTS_TS_natural)
                                KIE2 = (FTS_GS_natural / FTS_GS_isotope)
                                KIE = KIE2 * KIE1

                    #calculate KIE with tunneling
                    KIE_tunneling = KIE * qt_TS_natural / qt_TS_isotope

                    #convert to strings
                    gs_num = map(str, gs_num)
                    ts_num = map(str, ts_num)
                    isotope_mass = map(str, isotope_mass)
                    KIE = str(KIE)
                    KIE_tunneling = str(KIE_tunneling)

                    each_temp = ''.join(each_temp)
                    pressure = ''.join(pressure)
                    scale_factor = ''.join(scale_factor)
                    chem_sym = ', '.join(chem_sym)
                    gs_num = ', '.join(gs_num)
                    ts_num = ', '.join(ts_num)
                    isotope_mass = ', '.join(isotope_mass)
                    
                    output.append(each_GS)
                    output.append(each_TS)
                    output.append(each_temp)
                    output.append(pressure)
                    output.append(scale_factor)
                    output.append(chem_sym)
                    output.append(gs_num)
                    output.append(ts_num)
                    output.append(isotope_mass)
                    output.append(KIE)
                    output.append(KIE_tunneling)
                    master_output.append(output)

                    outputwriter.writerow([each_GS] + [each_TS] + [each_temp] + [pressure] + [scale_factor] + [chem_sym] + [gs_num] + [ts_num] + [isotope_mass] + [KIE] + [KIE_tunneling])
        csvfile.close()

        print '      '.join(heading)
        for each_entry in master_output:
            print '       '.join(each_entry)
        
        print "\n", "All KIE's have been calculated!", "\n"

def parse_txt_fchk(text_file):
    """
        Input: Text file (see def no_input for format)
        Output: Each input for the four separate freqchks
        Preparing temperary text files that will be fed into the 'freqchk' command
    """
    IN_TEXT = open(text_file)
    all_lines = IN_TEXT.readlines()

    lines = []
    for i in range(0, len(all_lines)):
        if all_lines[i] == "Ground state file (can be multiple):\n":
            Q1 = i
        if all_lines[i] == "Transition state file (can be multiple):\n":
            Q2 = i
       	if all_lines[i] == "Do you want to write hyperchem files? [No]\n":
            Q3 = i
        if all_lines[i] == "Temperature (in Kelvin) (can be multiple):\n":
    	    Q4 = i
        if all_lines[i] == "Pressure (in atm) (only one):\n":
            Q5 = i
        if all_lines[i] == "Scaling factor for frequencies during thermochemistry: (only one)\n":
            Q6 = i
        if all_lines[i] == "Project out gradient direction [No]:\n":
            Q7 = i
    	if all_lines[i] == "Atom symbol - ground state atom number - transition state atom number - isotope\n":
	        Q8 = i
        if all_lines[i] == "Put at least one space between each input. DO NOT use any other characters to signify spaces!\n":
            Q9 = i

    GS_chkpt_file = []
    TS_chkpt_file = []
    hyperchem_files = []
    temp = []
    pressure = []
    scale_factor = []
    iso_masses = []
    gradient_direction = []
    isotope_changes = []
    #Pulls data from text file and enters them into respective list
    for j in range(0, len(all_lines)):
        if j > Q1 and j < Q2:
            GS_chkpt_file.append(all_lines[j])
            GS_chkpt_file = map(lambda s: s.strip(), GS_chkpt_file)
    	if j > Q2 and j < Q3:
            TS_chkpt_file.append(all_lines[j])
            TS_chkpt_file = map(lambda s: s.strip(), TS_chkpt_file)
        if j > Q3 and j < Q4:
            hyperchem_files.append(all_lines[j])
            hyperchem_files = map(lambda s: s.strip(), hyperchem_files)
        if j > Q4 and j < Q5:
            temp.append(all_lines[j])
            temp = map(lambda s: s.strip(), temp)
        if j > Q5 and j < Q6:
            pressure.append(all_lines[j])
            pressure = map(lambda s: s.strip(), pressure)
        if j > Q6 and j < Q7:
            scale_factor.append(all_lines[j])
            scale_factor = map(lambda s: s.strip(), scale_factor)
        if j > Q7 and j < Q8 and j != Q8:
            gradient_direction.append(all_lines[j])
            gradient_direction = map(lambda s: s.strip(), gradient_direction)
        if j > Q9:
            all_lines[j] = all_lines[j].split()
            isotope_changes.append(all_lines[j])

    merge = list(itertools.chain.from_iterable(isotope_changes)) #Combines all lists of lists into a single list
    isotope_changes = [isotope_changes[x:x+4] for x in xrange(0, len(isotope_changes), 4)] #'Chunks' the list into lists of 4 (atom symbol, GS atom num, TS atom num, isotope mass) - makes list of lists
    isotope_changes = list(itertools.chain(*isotope_changes)) #Make an iterator that returns elements from the first iterable until it is exhausted, then proceeds to the next iterable, until all of the iterables are exhausted
    isotope_changes = sorted(isotope_changes, key=itemgetter(1)) #Sorts the lists (within the bigger list) based on the GS atom num
    isotope_changes = [[try_int(x) for x in lst] for lst in isotope_changes] #Calls 'try_int' module to turn any values into integers possible

    #Find the number of atoms in the GS checkpoint file
    for j in GS_chkpt_file:
        IN_FILE_GS = open(GS_chkpt_file[0])
        number_atoms_GS = []
        for line in IN_FILE_GS:
            if "Number of atoms" in line:
                number_atoms_GS.append(line.split())
                number_atoms_GS = ''.join(number_atoms_GS[0]) #Joins all characters in the line together to form "NumberofatomsI#" where # is the number of atoms
                number_atoms_GS = ''.join(number_atoms_GS[14:]) #First 14 characters are "NumberofatomsI" (leaves just the number of atoms)
                number_atoms_GS = int(number_atoms_GS) #Makes the number an integer rather than a string

    #Find the number of atoms in the TS checkpoint file
    i = 0
    for j in TS_chkpt_file:
        IN_FILE_TS = open(TS_chkpt_file[0])
        number_atoms_TS = []
        for line in IN_FILE_TS:
            if "Number of atoms" in line:
                number_atoms_TS.append(line.split())
                number_atoms_TS = ''.join(number_atoms_TS[0]) #Joins all characters in the line together to form "NumberofatomsI#" where # is the number of atoms
                number_atoms_TS = ''.join(number_atoms_TS[14:]) #First 14 characters are "NumberofatomsI" (leaves just the number of atoms)
                number_atoms_TS = int(number_atoms_TS) #Makes the number an integer rather than a string

    IN_TEXT.close()
    IN_FILE_GS.close()
    IN_FILE_TS.close()
    return GS_chkpt_file, TS_chkpt_file, hyperchem_files, temp, pressure, scale_factor, gradient_direction, isotope_changes, number_atoms_GS, number_atoms_TS

def parse_txt_log(text_file):
    """
        Input: Text file (see def no_input for format)
        Output: Each input for the four separate freqchks
        Preparing temperary text files that will be fed into the 'freqchk' command
    """
    IN_TEXT = open(text_file)
    all_lines = IN_TEXT.readlines()
    
    lines = []
    for i in range(0, len(all_lines)):
        if all_lines[i] == "Unlabeled transition state file:\n":
            Q1 = i
        if all_lines[i] == "Labeled transition state file:\n":
            Q2 = i
        if all_lines[i] == "Unlabeled ground state file:\n":
            Q3 = i
        if all_lines[i] == "Labeled ground state file:\n":
            Q4 = i
        if all_lines[i] == "Temperature(in Kelvin):\n":
            Q5 = i
        if all_lines[i] == "Scaling factor:\n":
            Q6 = i

    TS_unlabeled = []
    TS_labeled = []
    GS_unlabeled = []
    GS_labeled = []
    scale_factor = []
    temp = []

    #Pulls data from text file and enters them into respective list
    for j in range(0, len(all_lines)):
        if j > Q1 and j < Q2:
            TS_unlabeled.append(all_lines[j])
            TS_unlabeled = map(lambda s: s.strip(), TS_unlabeled)
        if j > Q2 and j < Q3:
            TS_labeled.append(all_lines[j])
            TS_labeled = map(lambda s: s.strip(), TS_labeled)
        if j > Q3 and j < Q4:
            GS_unlabeled.append(all_lines[j])
            GS_unlabeled = map(lambda s: s.strip(), GS_unlabeled)
        if j > Q4 and j < Q5:
            GS_labeled.append(all_lines[j])
            GS_labeled = map(lambda s: s.strip(), GS_labeled)
        if j > Q5 and j < Q6:
            temp.append(all_lines[j])
            temp = map(lambda s: s.strip(), temp)
        if j > Q6:
            all_lines[j] = all_lines[j].split()
            scale_factor.append(all_lines[j])
    IN_TEXT.close()
    return TS_unlabeled, TS_labeled, GS_unlabeled, GS_labeled, temp, scale_factor

def run_freqchk_TS_no_marker(chkpt_file, hyperchem, temp, pressure, scale_factor, gradient):
    """
        Input: Checkpoint file name and typical answers to 'freqchk' prompts with the exception of changing the isotopes
        Output: file with freqchk output
        Running 'freqchk' using text file created below
    """
    iso_masses = "Yes" #forces freqchk to use natural isotopes
    
    #Turning a list into strings
    hyperchem = ''.join(hyperchem)
    temp = ''.join(temp)
    pressure = ''.join(pressure)
    scale_factor = ''.join(scale_factor)
    gradient = ''.join(gradient)
    
    #Create and write to a temporary text file that will act as the input to 'freqchk'
    TEMP_FREQCHK_FILE = open("temp_freqchk_file1.txt","w")
    TEMP_FREQCHK_FILE.write(chkpt_file + "\n" + hyperchem + "\n" + temp + "\n" + pressure + "\n" + scale_factor + "\n" + iso_masses + "\n" + gradient)
    TEMP_FREQCHK_FILE.close()
    
    #Run freqchk with temporary text file as input
    os.system("freqchk " + "< temp_freqchk_file1.txt" "> freq_TS_no_marker.txt")
    
    #Remove temporary text file
    os.system("rm temp_freqchk_file1.txt")

def run_freqchk_GS_no_marker(chkpt_file, hyperchem, temp, pressure, scale_factor, gradient):
    """
        Input: Checkpoint file name and typical answers to 'freqchk' prompts with the exception of changing the isotopes
        Output: file with freqchk output
        Running 'freqchk' using text file created below
    """
    iso_masses = "Yes" #forces freqchk to use natural isotopes

    #Turning a list into strings
    hyperchem = ''.join(hyperchem)
    temp = ''.join(temp)
    pressure = ''.join(pressure)
    scale_factor = ''.join(scale_factor)
    gradient = ''.join(gradient)

    #Create and write to a temporary text file that will act as the input to 'freqchk'
    TEMP_FREQCHK_FILE = open("temp_freqchk_file2.txt","w")
    TEMP_FREQCHK_FILE.write(chkpt_file + "\n" + hyperchem + "\n" + temp + "\n" + pressure + "\n" + scale_factor + "\n" + iso_masses + "\n" + gradient)
    TEMP_FREQCHK_FILE.close()

    #Run freqchk with temporary text file as input
    os.system("freqchk " + "< temp_freqchk_file2.txt" "> freq_GS_no_marker.txt")

    #Remove temporary text file
    os.system("rm temp_freqchk_file2.txt")

def run_freqchk_TS_marker(chkpt_file, hyperchem, temp, pressure, scale_factor, gradient, atom_num, iso_mass, num_atoms):
    """
        Input: Checkpoint file name and the typical answers to 'freqchk' prompts with the exception of changing the isotopes
        Output: file with freqchk output
        Creating the input for main body of code
    """
    iso_masses = "No" #forces freqchk to use changes to isotopes
    
    # Build the input for the changes to isotopes
    #Blanks entered for atoms that are not changing
    iso_input = []
    count = 1
    while count / (num_atoms+1) != 1:
        if count in atom_num:
            i = atom_num.index(count)
            count += 1
            iso_input.append(str(iso_mass[i]))
        else:
            count += 1
            iso_input.append("")
    iso_input = [iso_input[x:x+num_atoms] for x in xrange(0, len(iso_input), num_atoms)] #'Chunks' the list into a list the length of num_atoms

    #Turning a list into strings
    for each_iso_input in iso_input:
        hyperchem = ''.join(hyperchem)
        temp = ''.join(temp)
        pressure = ''.join(pressure)
        scale_factor = ''.join(scale_factor)
        iso_masses = ''.join(iso_masses)
        gradient = ''.join(gradient)

        #Create and write to a temporary text file that will act as the input to 'freqchk'
        TEMP_FREQCHK_FILE = open("temp_freqchk_file3.txt","w")
        TEMP_FREQCHK_FILE.write(chkpt_file + "\n" + hyperchem + "\n" + temp + "\n" + pressure + "\n" + scale_factor + "\n" + iso_masses + "\n" + "\n".join(each_iso_input) + "\n" + gradient)
        TEMP_FREQCHK_FILE.close()

        #Run freqchk with temporary text file as input
        os.system("freqchk " + "< temp_freqchk_file3.txt" "> freq_TS_marker.txt")

        #Remove temporary text file
        os.system("rm temp_freqchk_file3.txt")

def run_freqchk_GS_marker(chkpt_file, hyperchem, temp, pressure, scale_factor, gradient, atom_num, iso_mass, num_atoms):
    """
        Input: Checkpoint file name and the typical answers to 'freqchk' prompts with the exception of changing the isotopes
        Output: file with freqchk output
        Creating the input for main body of code
    """
    iso_masses = "No" #forces freqchk to use changes to isotopes

    # Build the input for the changes to isotopes
    #Blanks entered for atoms that are not changing
    iso_input = []
    count = 1
    while count / (num_atoms+1) != 1:
        if count in atom_num:
            i = atom_num.index(count)
            count += 1
            iso_input.append(str(iso_mass[i]))
        else:
            count += 1
            iso_input.append("")
    iso_input = [iso_input[x:x+num_atoms] for x in xrange(0, len(iso_input), num_atoms)] #'Chunks' the list into a list the length of num_atoms

    #Turning a list into strings
    for each_iso_input in iso_input:
        hyperchem = ''.join(hyperchem)
        temp = ''.join(temp)
        pressure = ''.join(pressure)
        scale_factor = ''.join(scale_factor)
        iso_masses = ''.join(iso_masses)
        gradient = ''.join(gradient)

        #Create and write to a temporary text file that will act as the input to 'freqchk'
        TEMP_FREQCHK_FILE = open("temp_freqchk_file4.txt","w")
        TEMP_FREQCHK_FILE.write(chkpt_file + "\n" + hyperchem + "\n" + temp + "\n" + pressure + "\n" + scale_factor + "\n" + iso_masses + "\n" + "\n".join(each_iso_input) + "\n" + gradient)
        TEMP_FREQCHK_FILE.close()
        
        #Run freqchk with temporary text file as input
        os.system("freqchk " + "< temp_freqchk_file4.txt" "> freq_GS_marker.txt")

        #Remove temporary text file
        os.system("rm temp_freqchk_file4.txt")

def get_frequencies_Bigeleisen(log_file_name):
    #Text
    frequency = []
    
    with open(log_file_name, 'r+') as f:
        lines = f.readlines()
        for i in range(0, len(lines)):
            line = lines[i]
            if "Frequencies" in line:
                frequency.extend(line.split()) #Looks like ['Frequencies', '--', '#', '#', '#']
                frequency.remove('Frequencies') #Leaves ['--', '#', '#', '#']
                frequency.remove('--') #Leaves ['#', '#', '#']
    return frequency

def get_frequencies(log_file_name, hart_part, kcal):
    """
        Input: Output of 'freqchk' (see run_freqchk modules above)
        Output: all frequencies in one list
        Populate frequencies
    """

    HF_all = []
    HF = []
    Rot_all = []
    Erot = []
    Srot = []
    Trans_all = []
    Etrans = []
    Strans = []
    frequency = []

    with open(log_file_name, 'r+') as f:
        lines = f.readlines()
        for i in range(0, len(lines)):
            line = lines[i]
            if "SCF Done:" in line:
                HF_all.append(line.split())
                HF = HF_all[0][4]
            if "E (Thermal)" in line:
                Rot_all.append(lines[i+5].split())
                Trans_all.append(lines[i+4].split())
                Erot = Rot_all[0][1]
                Erot = float(Erot)/(hart_part/kcal)
                Srot = Rot_all[0][3]
                Srot = float(Srot)/(1000*hart_part/kcal)
                Etrans = Trans_all[0][1]
                Etrans = float(Etrans)/(hart_part/kcal)
                Strans = Trans_all[0][3]
                Strans = float(Strans)/(1000*hart_part/kcal)
            if "Frequencies" in line:
                frequency.extend(line.split()) #Looks like ['Frequencies', '--', '#', '#', '#']
                frequency.remove('Frequencies') #Leaves ['--', '#', '#', '#']
                frequency.remove('--') #Leaves ['#', '#', '#']

    return frequency, HF, Erot, Srot, Etrans, Strans

def calc_u(frequency, temp, scale):
    """
    Input: array with frequencies (see get_frequencies module)
    Output: values of u
    Step 2 to calculate KIE
    """
    temp = int(temp) #Makes integer
    scale = ''.join(scale)
    scale = float(scale) #Makes integer
    h = 6.626e-34
    K = 1.381e-23
    hKT = h/(K*temp)*30000000000 #Equivalent to h/KT

    u = []
    for each_freq in frequency:
        u = [x * hKT for x in frequency] #frequency * h/KT * scale_factor = u
        u = [x * scale for x in u]
    return u

def exp_u_half(u):
    """
        Input: u array (see calc_u module)
        Output: values of exp(u/2)
        Step 3 to calculate KIE
    """
    e = 2.71828182845904

    exp = []
    for each_u in u:
        exp = [x/2 for x in u] #exp(u/2)
        exp = [e**x for x in exp ]
    return exp

def calc_one_minus_exp(u):
    """
        Input: u array
        Output: values of 1-exp(-u)
        Step 4 to calculate KIE
    """
    e = 2.71828182845904
    one_minus_exp = []
    for each_u in u:
        one_minus_exp = [(-1)*x for x in u] # 1-exp(-u)
        one_minus_exp = [e**x for x in one_minus_exp]
        one_minus_exp = [1-x for x in one_minus_exp]
    return one_minus_exp

def calc_prod(u, exp, one_minus_exp):
    """
        Input: u array, exp_u_half array, one_minus_exp array
        Output: values of product
        Step 5 to calcualte KIE
    """
    prod = []
    prod = np.multiply(exp, one_minus_exp) #prod = (exp(u/2) * 1-exp(-u)) / u
    prod = np.divide(prod, u)
    return prod

def calc_FTS_TS(prod, neg_u):
    """
        Input: prod array and negative frequency
        Output: FTS value for transition state
        Step 6 to calculate KIE
    """
    FTS = np.prod(np.array(prod)) #Multiply all prod together and divide by the negative u value
    FTS = FTS / neg_u
    return FTS

def calc_FTS(prod):
    """
        Input: prod array
        Output: FTS value for ground state
        Step 6.5 to calculate KIE
    """
    FTS = np.prod(np.array(prod)) #Multiply all prod together
    return FTS

def calc_qt(u):
    """
        Input: u array
        Output: qt value
        Step 7 to calculate KIE
    """
    qt = u/2 #(u/2)/sin(u/2) where u is the negative u value
    sine_qt = u/2
    sine_qt = mt.sin(sine_qt)
    qt = qt/sine_qt
    return qt

def calc_IE_Bigeleisen(GS_unlabeled, GS_labeled, TS_unlabeled, TS_labeled, max_temp, min_temp, scale_factor, step_size):
    min_temp = float(min_temp)
    max_temp = float(max_temp)
    step_size = float(step_size)

    output_BM = []
    y_BM = []
    x_BM = []

    #get frequencies from .txt file
    frequency_TS_natural = map(float, get_frequencies_Bigeleisen(TS_unlabeled))
    frequency_TS_isotope = map(float, get_frequencies_Bigeleisen(TS_labeled))
    frequency_GS_natural = map(float, get_frequencies_Bigeleisen(GS_unlabeled))
    frequency_GS_isotope = map(float, get_frequencies_Bigeleisen(GS_labeled))

    step = ((max_temp - min_temp) / step_size) + 1

    interval = np.linspace(min_temp, max_temp, num=step)

    for x in interval:

        # creat array with u values
        u_TS_natural = np.array(calc_u(frequency_TS_natural, x, scale_factor))
        u_TS_isotope = np.array(calc_u(frequency_TS_isotope, x, scale_factor))
        u_GS_natural = np.array(calc_u(frequency_GS_natural, x, scale_factor))
        u_GS_isotope = np.array(calc_u(frequency_GS_isotope, x, scale_factor))
        if u_TS_natural[0] < 0:
            u_neg_TS_natural= u_TS_natural[0] #negative u value
            u_TS_natural = u_TS_natural[1:] #allows calculation on all u values that are positive
            u_neg_TS_isotope = u_TS_isotope[0] #negative u value
            u_TS_isotope = u_TS_isotope[1:] #allows calculation on all u values that are positive
        else:
            u_neg_TS_natural = []
            u_neg_TS_isotope = []

        # create array with exp(u/2) values
        exp_TS_natural = np.array(exp_u_half(u_TS_natural))
        exp_TS_isotope = np.array(exp_u_half(u_TS_isotope))
        exp_GS_natural = np.array(exp_u_half(u_GS_natural))
        exp_GS_isotope = np.array(exp_u_half(u_GS_isotope))

        # create array with 1-exp(-u) values
        one_minus_exp_TS_natural = np.array(calc_one_minus_exp(u_TS_natural))
        one_minus_exp_TS_isotope = np.array(calc_one_minus_exp(u_TS_isotope))
        one_minus_exp_GS_natural = np.array(calc_one_minus_exp(u_GS_natural))
        one_minus_exp_GS_isotope = np.array(calc_one_minus_exp(u_GS_isotope))

        # create array with prod values
        prod_TS_natural = np.array(calc_prod(u_TS_natural, exp_TS_natural, one_minus_exp_TS_natural))
        prod_TS_isotope = np.array(calc_prod(u_TS_isotope, exp_TS_isotope, one_minus_exp_TS_isotope))
        prod_GS_natural = np.array(calc_prod(u_GS_natural, exp_GS_natural, one_minus_exp_GS_natural))
        prod_GS_isotope = np.array(calc_prod(u_GS_isotope, exp_GS_isotope, one_minus_exp_GS_isotope))

        # calculate FTS
        if u_neg_TS_natural:
            FTS_TS_natural = calc_FTS_TS(prod_TS_natural, u_neg_TS_natural)
            FTS_TS_isotope = calc_FTS_TS(prod_TS_isotope, u_neg_TS_isotope)
            FTS_GS_natural = calc_FTS(prod_GS_natural)
            FTS_GS_isotope = calc_FTS(prod_GS_isotope)
        else:
            FTS_TS_natural = calc_FTS(prod_TS_natural)
            FTS_TS_isotope = calc_FTS(prod_TS_isotope)
            FTS_GS_natural = calc_FTS(prod_GS_natural)
            FTS_GS_isotope = calc_FTS(prod_GS_isotope)

        # calculate qt for TS
        if u_neg_TS_natural:
            qt_TS_natural = calc_qt(u_neg_TS_natural)
            qt_TS_isotope = calc_qt(u_neg_TS_isotope)
        else:
            qt_TS_natural = calc_qt(u_TS_natural[0])
            qt_TS_isotope = calc_qt(u_TS_isotope[0])

        # calculate KIE
        KIE1 = (FTS_TS_isotope / FTS_TS_natural)
        KIE2 = (FTS_GS_natural / FTS_GS_isotope)
        KIE = KIE2 * KIE1
        KIE_tunneling = KIE * qt_TS_natural / qt_TS_isotope

        y_BM.append(x)
        x_BM.append(KIE)

        print "BM", "\t", x, "\t", KIE, "\t", KIE_tunneling

    return y_BM, x_BM

def nu_bar_calc(frequencies, scale_factor):
    #nu_bar calculation
    nu_bar = []
    scale_factor = ''.join(scale_factor)
    for each_freq in frequencies:
        nu_bar_calc = float(each_freq) * float(scale_factor) #Scale factor
        nu_bar.append(nu_bar_calc)
    return nu_bar

def mu_calc(nu_bar, h, c, k, temp):
    #mu calculation
    mu = []
    for each_nu in nu_bar:
        numerator = h * c * each_nu #h*c*nu_bar
        denom = k * temp #k*temp
        mu_calc = numerator / denom
        mu.append(mu_calc)
    return mu

def calc_exp_mu_neg_1(mu):
    #exp^mu-1 calculation
    exp_mu_neg_1 = []
    for each_mu in mu:
        calculation = mt.exp(each_mu)
        calculation = calculation - 1
        exp_mu_neg_1.append(calculation)
    return exp_mu_neg_1

def calc_Hvib(nu_bar, e_mu_neg_1, k, N, h, c, temp):
    Hvib = []
    for each_mu in e_mu_neg_1:
        i = e_mu_neg_1.index(each_mu)
        numbers = N * h * c #N * h * c
        fraction = nu_bar[i] / each_mu
        calculation = numbers * fraction
        Hvib.append(calculation)
    return Hvib

def calc_e_mu(mu):
    #Calculate 1-e^(-mu)
    e_mu = []
    for each_mu in mu:
        neg = each_mu * -1
        calculation = mt.exp(neg)
        one_minus = 1 - calculation
        e_mu.append(one_minus)
    return e_mu

def calc_Svib(mu, exp_mu_neg_1, e_mu):
    #Calculate Svib/R
    Svib = []
    for each_mu in mu:
        i = mu.index(each_mu)
        fraction = each_mu / exp_mu_neg_1[i]
        ln = np.log(e_mu[i])
        calculation = fraction - ln
        Svib.append(calculation)
    return Svib

def calc_ZPVE(nu_bar, N, h, c, hart_part):
    #Calculate ZPVE in kcal/mol
    total = sum(nu_bar)
    calculation = 0.5 * total * N * h * c # N * h * c
    calculation = calculation / hart_part
    return calculation

def calc_Hvib_final(Hvib, hart_part):
    #Calculate final Hvib in kcal/mol
    calculation = sum(Hvib)
    calculation = calculation / hart_part
    return calculation

def calc_Svib_final(Svib, N, k, hart_part):
    #Calculate final Svib in kcal/mol
    product = sum(Svib)
    calculation = N * k * product #N * k * sum(Svib)
    calculation = calculation / hart_part
    return calculation

def calc_Hcorr(ZPVE, Hvib, Erot, Etrans, N, k, temp, hart_part):
    multiplication = N * k * temp
    fraction = multiplication / hart_part
    add = ZPVE + Hvib + Erot + Etrans + fraction
    return add

def calc_enthalpy(Hcorr, HF):
    #Calculate enthalpy
    calculation = float(HF) + Hcorr
    return calculation

def calc_Gcorr(Hcorr, temp, Svib, Srot, Strans):
    add = Svib + Srot + Strans
    multiply = temp * add
    minus = Hcorr - multiply
    return minus

def calc_Gibbs(Gcorr, E):
    calculation = float(E) + Gcorr
    return calculation

def calc_DeltaG(TS, GS):
    calculation = (GS - TS) * 627.5095
    return calculation

def calc_DeltaDeltaG(unlabeled, labeled):
    calculation = unlabeled - labeled
    return calculation

def calc_isotope_effect(DDG, temp):
    multiplication = 0.00198588 * temp
    division = (-1 * DDG) / multiplication
    exponent = mt.exp(division)
    calculation = mt.pow(exponent, -1)
    return calculation

def calc_isotope_tun(IE, tun_un, tun_lab):
    #Calculate IE with tunneling corrections
    calculation = IE * tun_un/tun_lab
    return calculation

def neg_u_calc(neg_freq, scale, temp):
    #Calculate neg u
    scale = ''.join(scale)
    scale = float(scale) #Makes integer
    hkT = 6.626E-34/(1.381E-23*temp)*30000000000
    calculation = neg_freq * scale * hkT
    return calculation

def calc_IE_rigid_rotor_method(GS_unlabeled, GS_labeled, TS_unlabeled, TS_labeled, max_temp, min_temp, scale_factor, step_size):
    #Constants
    h = 6.626076e-34
    c = 2.99792e10
    k = 1.38066e-23
    N = 6.02214e23
    hartree = 4.3597e-18
    hart_part = hartree * N
    kcal = 4184
    min_temp = float(min_temp)
    max_temp = float(max_temp)
    step_size = float(step_size)

    thermal_output = []

    y_RR = []
    x_RR = []

    #Extract frequencies from log file
    frequencies_TS_unlabeled, HF_TS_un, Erot_TS_un, Srot_TS_un, Etrans_TS_un, Strans_TS_un = get_frequencies(TS_unlabeled, hart_part, kcal)
    frequencies_TS_labeled, HF_TS_lab, Erot_TS_lab, Srot_TS_lab, Etrans_TS_lab, Strans_TS_lab = get_frequencies(TS_labeled, hart_part, kcal)
    frequencies_GS_unlabeled, HF_GS_un, Erot_GS_un, Srot_GS_un, Etrans_GS_un, Strans_GS_un = get_frequencies(GS_unlabeled, hart_part, kcal)
    frequencies_GS_labeled, HF_GS_lab, Erot_GS_lab, Srot_GS_lab, Etrans_GS_lab, Strans_GS_lab = get_frequencies(GS_labeled, hart_part, kcal)

    neg_frequency_unlabeled = ''.join(frequencies_TS_unlabeled[0])
    neg_frequency_labeled = ''.join(frequencies_TS_unlabeled[0])

    #Removing neg frequency from calculations
    if float(neg_frequency_unlabeled) < 0:
        frequencies_TS_unlabeled = frequencies_TS_unlabeled[1:]
        frequencies_TS_labeled = frequencies_TS_labeled[1:]

    step = ((max_temp - min_temp) / step_size) + 1

    interval = np.linspace(min_temp, max_temp, num=step)

    for x in interval:
        #Calculate nu-bar
        nu_bar_TS_unlabeled = nu_bar_calc(frequencies_TS_unlabeled, scale_factor)
        nu_bar_TS_labeled = nu_bar_calc(frequencies_TS_labeled, scale_factor)
        nu_bar_GS_unlabeled = nu_bar_calc(frequencies_GS_unlabeled, scale_factor)
        nu_bar_GS_labeled = nu_bar_calc(frequencies_GS_labeled, scale_factor)

        #Calculate mu
        mu_TS_unlabeled = mu_calc(nu_bar_TS_unlabeled, h, c, k, x)
        mu_TS_labeled = mu_calc(nu_bar_TS_labeled, h, c, k, x)
        mu_GS_unlabeled = mu_calc(nu_bar_GS_unlabeled, h, c, k, x)
        mu_GS_labeled = mu_calc(nu_bar_GS_labeled, h, c, k, x)

        #Calculate e^mu-1
        exp_mu_neg_1_TS_unlabeled = calc_exp_mu_neg_1(mu_TS_unlabeled)
        exp_mu_neg_1_TS_labeled = calc_exp_mu_neg_1(mu_TS_labeled)
        exp_mu_neg_1_GS_unlabeled = calc_exp_mu_neg_1(mu_GS_unlabeled)
        exp_mu_neg_1_GS_labeled = calc_exp_mu_neg_1(mu_GS_labeled)

        #Calculate Hvib
        Hvib_TS_unlabeled = calc_Hvib(nu_bar_TS_unlabeled, exp_mu_neg_1_TS_unlabeled, k, N, h, c, x)
        Hvib_TS_labeled = calc_Hvib(nu_bar_TS_labeled, exp_mu_neg_1_TS_labeled, k, N, h, c, x)
        Hvib_GS_unlabeled = calc_Hvib(nu_bar_GS_unlabeled, exp_mu_neg_1_GS_unlabeled, k, N, h, c, x)
        Hvib_GS_labeled = calc_Hvib(nu_bar_GS_labeled, exp_mu_neg_1_GS_labeled, k, N, h, c, x)

        #Calculate 1-e^(-mu)
        e_mu_TS_unlabeled = calc_e_mu(mu_TS_unlabeled)
        e_mu_TS_labeled = calc_e_mu(mu_TS_labeled)
        e_mu_GS_unlabeled = calc_e_mu(mu_GS_unlabeled)
        e_mu_GS_labeled = calc_e_mu(mu_GS_labeled)

        #Calculate Svib/R
        Svib_TS_unlabeled = calc_Svib(mu_TS_unlabeled, exp_mu_neg_1_TS_unlabeled, e_mu_TS_unlabeled)
        Svib_TS_labeled = calc_Svib(mu_TS_labeled, exp_mu_neg_1_TS_labeled, e_mu_TS_labeled)
        Svib_GS_unlabeled = calc_Svib(mu_GS_unlabeled, exp_mu_neg_1_GS_unlabeled, e_mu_GS_unlabeled)
        Svib_GS_labeled = calc_Svib(mu_GS_labeled, exp_mu_neg_1_GS_labeled, e_mu_GS_labeled)

        #Calculate ZPVE
        ZPVE_TS_unlabeled = calc_ZPVE(nu_bar_TS_unlabeled, N, h, c, hart_part)
        ZPVE_TS_labeled = calc_ZPVE(nu_bar_TS_labeled, N, h, c, hart_part)
        ZPVE_GS_unlabeled = calc_ZPVE(nu_bar_GS_unlabeled, N, h, c, hart_part)
        ZPVE_GS_labeled = calc_ZPVE(nu_bar_GS_labeled, N, h, c, hart_part)

        #Calculate Hvib
        Hvib_final_TS_unlabeled = calc_Hvib_final(Hvib_TS_unlabeled, hart_part)
        Hvib_final_TS_labeled = calc_Hvib_final(Hvib_TS_labeled, hart_part)
        Hvib_final_GS_unlabeled = calc_Hvib_final(Hvib_GS_unlabeled, hart_part)
        Hvib_final_GS_labeled = calc_Hvib_final(Hvib_GS_labeled, hart_part)

        #Calculate Svib
        Svib_final_TS_unlabeled = calc_Svib_final(Svib_TS_unlabeled, N, k, hart_part)
        Svib_final_TS_labeled = calc_Svib_final(Svib_TS_labeled, N, k, hart_part)
        Svib_final_GS_unlabeled = calc_Svib_final(Svib_GS_unlabeled, N, k, hart_part)
        Svib_final_GS_labeled = calc_Svib_final(Svib_GS_labeled, N, k, hart_part)

        #Calculate Hcorr
        Hcorr_TS_unlabeled = calc_Hcorr(ZPVE_TS_unlabeled, Hvib_final_TS_unlabeled, Erot_TS_un, Etrans_TS_un, N, k, x, hart_part)
        Hcorr_TS_labeled = calc_Hcorr(ZPVE_TS_labeled, Hvib_final_TS_labeled, Erot_TS_lab, Etrans_TS_lab, N, k, x, hart_part)
        Hcorr_GS_unlabeled = calc_Hcorr(ZPVE_GS_unlabeled, Hvib_final_GS_unlabeled, Erot_GS_un, Etrans_GS_un, N, k, x, hart_part)
        Hcorr_GS_labeled = calc_Hcorr(ZPVE_GS_labeled, Hvib_final_GS_labeled, Erot_GS_lab, Etrans_GS_lab, N, k, x, hart_part)

        #Calculate Enthalpy
        Enthalpy_TS_unlabeled = calc_enthalpy(Hcorr_TS_unlabeled, HF_TS_un)
        Enthalpy_TS_labeled = calc_enthalpy(Hcorr_TS_labeled, HF_TS_lab)
        Enthalpy_GS_unlabeled = calc_enthalpy(Hcorr_GS_unlabeled, HF_GS_un)
        Enthalpy_GS_labeled = calc_enthalpy(Hcorr_GS_labeled, HF_GS_lab)

        #Calculate Gcorr
        Gcorr_TS_unlabeled = calc_Gcorr(Hcorr_TS_unlabeled, x, Svib_final_TS_unlabeled, Srot_TS_un, Strans_TS_un)
        Gcorr_TS_labeled = calc_Gcorr(Hcorr_TS_labeled, x, Svib_final_TS_labeled, Srot_TS_lab, Strans_TS_lab)
        Gcorr_GS_unlabeled = calc_Gcorr(Hcorr_GS_unlabeled, x, Svib_final_GS_unlabeled, Srot_GS_un, Strans_GS_un)
        Gcorr_GS_labeled = calc_Gcorr(Hcorr_GS_labeled, x, Svib_final_GS_labeled, Srot_GS_lab, Strans_GS_lab)

        #Calculate Gibbs Free Energy
        Gibbs_TS_unlabeled = calc_Gibbs(Gcorr_TS_unlabeled, HF_TS_un)
        Gibbs_TS_labeled = calc_Gibbs(Gcorr_TS_labeled, HF_TS_lab)
        Gibbs_GS_unlabeled = calc_Gibbs(Gcorr_GS_unlabeled, HF_GS_un)
        Gibbs_GS_labeled = calc_Gibbs(Gcorr_GS_labeled, HF_GS_lab)

        #Calculate DeltaG
        DeltaG_unlabeled = calc_DeltaG(Gibbs_TS_unlabeled, Gibbs_GS_unlabeled)
        DeltaG_labeled = calc_DeltaG(Gibbs_TS_labeled, Gibbs_GS_labeled)

        #Calculate DeltaDeltaG
        DeltaDeltaG = calc_DeltaDeltaG(DeltaG_unlabeled, DeltaG_labeled)

        #Prep for tunneling corrections
        neg_u_lab = neg_u_calc(float(neg_frequency_labeled), scale_factor, x)
        neg_u_un = neg_u_calc(float(neg_frequency_unlabeled), scale_factor, x)

        #Tunneling correction
        tun_corr_labeled = calc_qt(neg_u_lab)
        tun_corr_unlabeled = calc_qt(neg_u_un)

        #Calculate Equilibrium constant
        isotope_effect = calc_isotope_effect(DeltaDeltaG, x)
        isotope_effect_tun = calc_isotope_tun(isotope_effect, tun_corr_unlabeled, tun_corr_labeled)

        y_RR.append(x)
        x_RR.append(isotope_effect)

        print "RR", "\t", x, "\t", isotope_effect, "\t", isotope_effect_tun

        #prepare for output
        x = str(x)
        Hcorr_GS_unlabeled = str(Hcorr_GS_unlabeled)
        Hcorr_GS_labeled = str(Hcorr_GS_labeled)
        Hcorr_TS_unlabeled = str(Hcorr_TS_unlabeled)
        Hcorr_TS_labeled = str(Hcorr_TS_labeled)
        Gcorr_GS_unlabeled = str(Gcorr_GS_unlabeled)
        Gcorr_GS_labeled = str(Gcorr_GS_labeled)
        Gcorr_TS_unlabeled = str(Gcorr_TS_unlabeled)
        Gcorr_TS_labeled = str(Gcorr_TS_labeled)

        thermal_output.append("GS_unlabeled")
        thermal_output.append(x)
        thermal_output.append(Hcorr_GS_unlabeled)
        thermal_output.append(Gcorr_GS_unlabeled)
        thermal_output.append("GS_labeled")
        thermal_output.append(x)
        thermal_output.append(Hcorr_GS_labeled)
        thermal_output.append(Gcorr_GS_labeled)
        thermal_output.append("TS_unlabeled")
        thermal_output.append(x)
        thermal_output.append(Hcorr_TS_unlabeled)
        thermal_output.append(Gcorr_TS_unlabeled)
        thermal_output.append("TS_labeled")
        thermal_output.append(x)
        thermal_output.append(Hcorr_TS_labeled)
        thermal_output.append(Gcorr_TS_labeled)

    return y_RR, x_RR, thermal_output

def calculate_IE(GS_un, GS_lab, TS_un, TS_lab, max_temp, min_temp, scale_factor, step_size, var):
    if GS_un == '' or GS_lab == '' or TS_un == '' or TS_lab == '':
        print "Please select four files before proceeding."
        return

    try:
        if float(min_temp) < 0 or float(max_temp) < 0 or float(scale_factor) < 0 or float(step_size) < 0:
            print "Cannot have negative values for temperatures and/or scaling factor!"
            return
    except ValueError:
        print "Temperatures and/or scaling factor should be a number."
        return

    print "\n", "ISOTOPE EFFECTS"
    print "--------------------------------------------------"
    print "Method", "\t", "Temp", "\t", "Isotope effect", "\t", "IE with tunneling", "\n", "--------------------------------------------------"

    y_BM, x_BM = calc_IE_Bigeleisen(GS_un, GS_lab, TS_un, TS_lab, max_temp, min_temp, scale_factor, step_size)

    y_RR, x_RR, thermal_output = calc_IE_rigid_rotor_method(GS_un, GS_lab, TS_un, TS_lab, max_temp, min_temp, scale_factor, step_size)

    print "--------------------------------------------------"

    if var == 1:
        Plot = pyplot.plot(x_BM, y_BM, label='BM Method')
        Plot = pyplot.plot(x_RR, y_RR, label='RR Method')
        pyplot.legend(loc='upper right')
        pyplot.title('Effects of temperature\non isotope effects')
        pyplot.xlabel('Isotope Effect')
        pyplot.ylabel('Temperature (K)')
        pyplot.show(Plot)

    print "\n", "THERMAL CORRECTIONS"
    print "------------------------------------------------------------"
    print "File", "\t\t", "Temp", "\t", "to Enthalpy", "\t", "to Gibbs Free Energy", "\n", "------------------------------------------------------------"

    for i in range(len(thermal_output)/4):
        print "\t".join(thermal_output[i*4:(i+1)*4])

    print "------------------------------------------------------------"

def main():
    window = GUI()

    """
    master_output = []
    if len(sys.argv) > 1:
        print "\n", "Please wait. Your KIE's are being calculated.", "\n"
    else:
        print "\n", "Please fill out the text file that was just populated in the current directory (timestamped) and input it in the command line following the script.", "\n"
        ERROR_OUTPUT = no_input() #Module that populates an empty text file for the user to fill out
        sys.exit()
    
    heading = ["Method", "TS_unlabled", "TS_labeled", "GS_unlabeled", "GS_labeled", "Temperature", "Scale_factor", "IE", "IE_tunneling"]
    print '\t'.join(heading)
    
    TS_unlabeled, TS_labeled, GS_unlabeled, GS_labeled, temp, scale_factor = parse_txt_log(sys.argv[1])

    i = 0
    j = 0
    for each_TS_un in TS_unlabeled:
        i = TS_unlabeled.index(each_TS_un)
        each_TS_lab = TS_labeled[i]
        for each_GS_un in GS_unlabeled:
            j = GS_unlabeled.index(each_GS_un)
            each_GS_lab = GS_labeled[j]
            for each_temp in temp:

                each_TS_un = ''.join(each_TS_un)
                each_TS_lab = ''.join(each_TS_lab)
                each_GS_un = ''.join(each_GS_un)
                each_GS_lab = ''.join(each_GS_lab)
                scale_factor = list(itertools.chain.from_iterable(scale_factor))
                each_temp = ''.join(each_temp)

                IE_Bigeleisen = calc_IE_Bigeleisen(each_TS_un, each_TS_lab, each_GS_un, each_GS_lab, each_temp, scale_factor)

                IE_rigid_rotor_method = calc_IE_rigid_rotor_method(each_TS_un, each_TS_lab, each_GS_un, each_GS_lab, each_temp, scale_factor)

    print "\n", "BM = Bigeleisen-Mayor; RR = Rigid-Rotor Harmonic Oscillator"
    """

if __name__=="__main__":
    main()
