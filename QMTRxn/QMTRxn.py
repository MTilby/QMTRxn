import ast
import csv
import os
import re
import extractor
import functions
import numpy as np
from argparse import ArgumentParser

### set current directory ###

current_dir = os.getcwd()

### Function for command line arguments ###

def GetArgs():
    """Get command line arguments"""

    parser = ArgumentParser()

    parser.add_argument("-s", "--software", type=str, choices=["ORCA","G16"], default="ORCA", help="software used for the calculation")
    parser.add_argument("-csv", "--csv", type=str, help="name of .csv file with information on data to extract")
    parser.add_argument("-p", "--project", type=str, help="name of the directory to iterate")
    parser.add_argument("-dir", "--directory", action="store_true", help="searches in subdirectories for ORCA output files")

    parser.add_argument("-iso", "--isotope", type=str, help="element to change mass of")
    parser.add_argument("-m", "--mass", type=float, help="mass to change a give element set by --isotope/-iso to")
    parser.add_argument("-amap", "--atom_mapping", type=dict, help="atom mapping to change atom in the ith position to the provided mass")

    parser.add_argument("-T", "--temperature", type=float, help="temperature to recalcuate thermochemistry")
    parser.add_argument("-atm", "--atm", type=float, help="atmosphere to recalcuate thermochemistry")
    parser.add_argument("-conc", "--concentration", type=float, help="concentration to recalcuate thermochemistry")
    parser.add_argument("-symno", "--symmetry_number", type=float, help="provided symmetry number of the molecule in thermochemistry recalculation")

    parser.add_argument("-fs", "--freq_scale", type=float, help="frequency scale factor to apply in thermochemistry recalculation")
    parser.add_argument("-inv", "--freq_invert", choices=['All', 'nonTS', None], default=None, help="invert immaginary frequencies in thermochemistry recalculation, the default is false") 
    parser.add_argument("-fco", "--freq_cutoff", type=float, help="cutoff value to convert low frquency numbers to a given value in thermochemistry recalculation")
    parser.add_argument("-lf", "--low_freq", type=float, help="value to convert low frquency numbers to a given value in thermochemistry recalculation")
    parser.add_argument("-vmap", "--freq_mapping", type=str, help="mapping to change the frequencies")

    parser.add_argument("-qhH", "--quasi_harmonic_enthalpy", action="store_true", help="switch on quasi harmonic approx. for enthalpy contribution in thermochemistry recalculation")
    parser.add_argument("-emet", "--entropy_method", type=str, choices=["Truhlar", "RRHO", None], default=None, help="correction method to apply to the entropy in thermochemistry recalculation, automatically uses quasi-RRHO")

    parser.add_argument("-v0H", "--v0H", type=float, help="frequncy value in cm-1 to use for quasi-RRHO enthalpy correction in thermochemistry recalculation")
    parser.add_argument("-alphaH", "--alphaH", type=float, help="alpha value to use for quasi-RRHO enthalpy in thermochemistry recalculation")
    parser.add_argument("-v0S", "--v0S", type=float, help="frequncy value in cm-1 to use for quasi-RRHO and Truhlar entropy correction in thermochemistry recalculation")
    parser.add_argument("-alphaS", "--alphaS", type=float, help="alpha value to use for quasi-RRHO entropy correction in thermochemistry recalculation")

    args = parser.parse_args()

    return args

if __name__ == "__main__":
    args = GetArgs()

if args.atm and args.concentration:
    print('Both atm and conc are given will use atmosphere')

if args.atom_mapping:
    s = re.sub(r'([{,]\s*)([A-Za-z_]\w*)(\s*:)', r'\1"\2"\3', args.atom_mapping)
    amap = ast.literal_eval(s)

else:
    amap = None

if args.freq_mapping:
    s = re.sub(r'([{,]\s*)([A-Za-z_]\w*)(\s*:)', r'\1"\2"\3', args.freq_mapping)
    vmap = ast.literal_eval(s)

else:
    vmap = None

if bool(args.isotope) != bool(args.mass):
    print(f"Requires specification of atom and a new mass")

freq_trig = [
    "freq_scale",
    "freq_invert",
    "freq_cutoff",
    "low_freq"
]

thermo_trig = [
    "quasi_harmonic_enthalpy",
    "entropy_method",
    "temperature",
    "atm",
    "concentration",
    "symmetry_number",
    "v0H",
    "alphaH",
    "v0S",
    "alphaS",
]

args.do_freq = any(getattr(args, name) not in (None, False)for name in freq_trig)
args.do_thermo = any(getattr(args, name) not in (None, False)for name in thermo_trig)

### Function for file extraction from .csv - potentially move to other script for clarity - and allow passing of other data types in a jupyter notebook ###

def csvfile(csvfile):
    """Read the CSV file and return a dictionary of calculation files and parameters."""
    with open(csvfile, newline='', encoding='utf-8-sig') as csvfile:
        reader = csv.DictReader(csvfile)
        csv_data = {}

        for row in reader:
            row = {key.lower().replace("_", "" ""): value if isinstance(value, str) else value for key, value in row.items()}
            filepath = os.path.join(current_dir,row.get("calc", "").strip())
            
            if os.path.isfile(f"{filepath}.property.txt"):
                csv_data[row.get("calc", "").strip()] = f"{filepath}.property.txt"

            elif os.path.isfile(f"{filepath}.out"):
                csv_data[row.get("calc", "").strip()] = f"{filepath}.out"

    return csv_data


### .csv handling ###

if args.csv:
    csv_filename = f"{args.csv}.csv" if not args.csv.endswith(".csv") else args.csv
    csv_file = os.path.join(current_dir, csv_filename)
    calc_files = csvfile(csv_file)

### set directories to use ###

else:

    if args.project:
        base_dir = os.path.join(current_dir, args.project)
    else:
        base_dir = current_dir

    if args.directory:
        roots = [os.path.join(base_dir, d) for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]
    else:
        roots = [base_dir]

    ### obtain data from ORCA output files ###

    calc_files = {}

    for root in roots:
        for file in os.listdir(root):

            if file.endswith(".property.txt"):
                calc = file.removesuffix(".property.txt")
                calc_files[calc] = os.path.join(root, file)

            elif file.endswith(".out") and file.count(".") == 1:
                calc = file.removesuffix(".out")
                calc_files[calc] = os.path.join(root, file)

            elif file.endswith(".log") and args.software == "G16":
                calc = file.removesuffix(".log")
                calc_files[calc] = os.path.join(root, file)

properties = {}
for calc in calc_files:
    energies= {}
    vib = None

    if args.software == "ORCA":
        calc_data = extractor.orca(calc, calc_files[calc])

    elif args.software =='G16':
        calc_data = extractor.g16(calc, calc_files[calc])

    xyz = calc_data.xyz()
    
    if args.isotope or args.atom_mapping:
        no_neg = sum(1 for v in calc_data.vib().values() if v < 0)
        xyz = functions.XYZ(xyz).isotope(getattr(args, "isotope", None), getattr(args, "mass", None), getattr(args, "atom_mapping", None))
        hess = calc_data.hessian()
        hess_mw = functions.hessian(xyz).mass_weight(hess) 
        vib = functions.hessian(xyz).frequency(hess_mw, no_neg)

    if args.do_thermo or args.do_freq:
        if not vib:
            vib = calc_data.vib()

        if args.do_freq:                
            vib = functions.frequencies(vib).correction(getattr(args, "freq_scale", None), getattr(args, "freq_invert", None), getattr(args, "freq_cutoff", None), getattr(args, "low_freq", None), vmap)

    elif vib:
        pass 

    else:
        energies = calc_data.energies() 
        

    if not energies:
        print(calc)
        energy, spin = calc_data.e_spin()

        energies = functions.thermo(xyz, vib, getattr(args, "temperature", None)).thermo_data(energy, spin, getattr(args, "quasi_harmonic_enthalpy", None), getattr(args, "v0H", None), getattr(args, "alphaH", None), 
            getattr(args, "atm", None), getattr(args, "concentration", None), getattr(args, "symmetry_number", None), getattr(args, "entropy_method", None), getattr(args, "v0S", None), getattr(args, "alphaS", None))

    properties[calc] = energies

print(properties)





