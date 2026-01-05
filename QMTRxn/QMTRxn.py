import csv
import os
from argparse import ArgumentParser

### set current directory ###

current_dir = os.getcwd()

### Function for command line arguments ###

def GetArgs():
    """Get command line arguments"""

    parser = ArgumentParser()

    parser.add_argument("-csv", "--csv", type=str, help="name of .csv file with information on data to extract")
    parser.add_argument("-p", "--project", type=str, help="name of the directory to iterate")
    parser.add_argument("-dir", "--directory", action="store_true", help="searches in subdirectories for ORCA output files")

    args = parser.parse_args()

    return args

if __name__ == "__main__":
    args = GetArgs()

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
            path = os.path.join(root, file)


            if file.endswith(".property.txt"):
                calc = file.removesuffix(".property.txt")
                calc_files[calc] = path

            elif file.endswith(".out") and file.count(".") == 1:
                calc = file.removesuffix(".out")
                calc_files.setdefault(calc, path)

print(calc_files)






