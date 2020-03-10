#!/usr/bin/python

import numpy as np

f = open("myhits_normal_charge_intein_seq.fa", "r").readlines()

def calculateNC(sequence):
    pos = sequence.count("R") + sequence.count("K") # arg + lys
    neg = sequence.count("D") + sequence.count("E") # asp + glu

    return pos - neg # charges

def calcMW(sequence):
    # Single to three letter aa code dictionary
    ot3 = {"I" : "ILE", "L" : "LEU", "V" : "VAL", "A" : "ALA", "G" : "GLY",
           "P" : "PRO", "F" : "PHE", "M" : "MET", "W" : "TRP", "Y" : "TYR",
           "H" : "HIS", "T" : "THR", "S" : "SER", "N" : "ASN", "Q" : "GLN",
           "D" : "ASP", "E" : "GLU", "K" : "LYS", "R" : "ARG", "C" : "CYS",
           "X" : "UNK", "U" : "SEC","-":"MIS"}

    # Residue weights in Da
    mass=dict([("ALA", 71.1), ("ARG", 156.2), ("ASN", 114.1),
               ("ASP", 115.1), ("CYS", 103.1), ("GLU", 129.1),
               ("GLN", 128.1), ("GLY", 57.1), ("HIS", 137.1),
               ("ILE", 113.2), ("LEU", 113.2), ("LYS", 128.2),
               ("MET", 131.2), ("PHE", 147.2), ("PRO", 97.1),
               ("SER", 87.1), ("THR", 101.1), ("TRP", 186.2),
               ("TYR", 163.2), ("VAL", 99.1), ("UNK", 110.0),
               ("SEC", 150.1),("MIS",0.001)])
    mw = 0.0
    for letter in sequence:
        mw = mass[ot3[letter]] + mw


    return mw

def calc_K_percent(sequence):
    K_num = sequence.count("K")
    resi_numb = len(sequence)
    return K_num / resi_numb

def calc_R_percent(sequence):
    R_num = sequence.count("R")
    resi_numb = len(sequence)
    return R_num / resi_numb

def calc_D_percent(sequence):
    D_num = sequence.count("D")
    resi_numb = len(sequence)
    return D_num / resi_numb

def calc_E_percent(sequence):
    E_num = sequence.count("E")
    resi_numb = len(sequence)
    return E_num / resi_numb
def calc_Charge_residue(sequence):
    pos = sequence.count("R") + sequence.count("K") # arg + lys
    neg = sequence.count("D") + sequence.count("E") # asp + glu

    return pos + neg # charges
def calcSASA(mw):
    return 6.3 * np.power(mw, 0.73) # square angstroms
name = ""
sequence = ""
#print("name\tNCD\tNCpercent\tPos_percent\tNega_percent\tKper\tRper\tDper\tEper")
for line in f:   
    if ">" in line and name != "" and sequence != "":
      print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(name,NCD,NCpercent,Pos_percent,Nega_percent,K_percent,R_percent,D_percent,E_percent,ABS_charge_percent))
      name = "" + line[1:-1]
      sequence = ""

    elif ">" in line :
         name = name + line[1:-1]
         pass

    elif ">" not in line:
        sequence = sequence + line.strip()
        net_charge = calculateNC(sequence)
        mw = calcMW(sequence)
        sasa = calcSASA(mw)
        NCD = net_charge / sasa
        NCpercent = net_charge / len(sequence)
        ABS_charge_percent = calc_Charge_residue(sequence)/len(sequence)
        K_percent = calc_K_percent(sequence)
        R_percent = calc_R_percent(sequence)
        D_percent = calc_D_percent(sequence)
        E_percent = calc_E_percent(sequence)
        Pos_percent = K_percent + R_percent
        Nega_percent = D_percent + E_percent
print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(name,NCD,NCpercent,Pos_percent,Nega_percent,K_percent,R_percent,D_percent,E_percent,ABS_charge_percent))    
