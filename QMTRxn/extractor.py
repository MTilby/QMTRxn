import constants

class orca:
    """Class to extract data from ORCA output files."""
    def __init__(self, file):
        self.file = file
    
    def extractor(self):
        """Extract relevant data from the ORCA output file."""

        thermo_data = {}

        if self.file.endswith(".property.txt"):
            with open(self.file, 'r') as file:
                lines = file.readlines()

            if not any("NORMAL TERMINATION" in line.upper() for line in lines):
                return thermo_data
            
            keywords = {
                '&VDW': 'vdw',
                '&TRANSENERGY': 'Utrans',
                '&ROTENERGY': 'Urot',
                '&VIBENERGY': 'Uvib',
                '&QEL': 'Qel',
                '&QTRANS': 'Qtrans',
                '&QROT': 'Qrot',
                '&QVIB': 'Qvib',
                '&TEMPERATURE': 'T',
                '&PRESSURE': 'P',
                '&TOTALMASS': 'Mass',
                '&FINALEN': 'E',
                '&ZPE': 'ZPE',
                '&INNERENERGYU': 'U',
                '&ENTHALPYH': 'H',
                '&ENTROPYS': 'ST',
                '&FREEENERGYG': 'G',
                '&corr': 'corr'
                }

            
            thermo_data.update({key: None for key in keywords.values()})

            for i in range(len(lines) - 1, -1, -1):
                line = lines[i]
                upper_line = line.upper()
                for keyword, key in keywords.items():
                    if thermo_data[key] is None and keyword in upper_line:
                        thermo_data[key] = float(line.split()[3])

                if "&NUMOFFREQS" in upper_line:
                    vib_count = int(line.split()[3])

                    for j in range(i + 5, i + 5 + vib_count):
                        if float(lines[j].split()[1]) < 0:
                            thermo_data['im_freq'] = float(lines[j].split()[1])
                            break


            if thermo_data['vdw'] is not None and thermo_data['E'] is not None:
                thermo_data['E'] += thermo_data['vdw']

            if thermo_data['G'] is not None and thermo_data['E'] is not None:
                thermo_data['corr'] = thermo_data['G'] - thermo_data['E']
            else:
                thermo_data['corr'] = None

            return thermo_data


        elif self.file.endswith(".out"):
            with open(self.file, 'r') as file:
                lines = file.readlines()

            if not any("ORCA TERMINATED NORMALLY"in line.upper().replace("*", "") for line in lines):
                return thermo_data
                   
            if any("the optimization did not converge"in line.lower() for line in lines):
                return thermo_data
            
            thermo_data.update({'vdw':None})

            keywords = {
                'thermaltranslationalcorrection...': 'Utrans',
                'thermalrotationalcorrection...': 'Urot',
                'thermalvibrationalcorrection...': 'Uvib',
                'Electronicentropy...': 'Qel',
                'translationalentropy...': 'Qtrans',
                'rotationalentropy...': 'Qrot',
                'vibrationalentropy...': 'Qvib',
                'pressure...': 'P',
                'totalmass...': 'Mass',
                'temperature...': 'T',
                'electronicenergy...': 'E',
                'zeropointenergy...': 'ZPE',
                'totalthermalenergy...': 'U',
                'totalenthalpy...': 'H',
                'totalentropycorrection...': 'ST',
                'finalgibbsfreeenergy...': 'G',
                'g-e(el)...': 'corr'
                }
                
            thermo_data.update({key: None for key in keywords.values()})
            
            for i in range(len(lines) - 1, -1, -1):
                line = lines[i]
                lower_line = line.lower().replace(" ","")
                for keyword, key in keywords.items():
                    if thermo_data[key] is None and keyword in lower_line:

                        thermo_data[key] = float(line.split('...')[1].split()[0])

                if 'imaginary mode' in line.replace('*',''):
                    thermo_data['im_freq'] = line.split()[1]

                if thermo_data['E'] is None and line.startswith('FINAL SINGLE POINT ENERGY'):
                    thermo_data['E'] = float(line.split()[4])

                if thermo_data['vdw'] is None and line.startswith('Dispersion correction'):
                    thermo_data['vdw'] = float(line.split()[2])

            return thermo_data
        
class g16:
    """Class to extract data from G16 output files."""
    def __init__(self, file):
        self.file = file
    
    def extractor(self):
        """Extract relevant data from the G16 output file."""

        thermo_data = {}

        if self.file.endswith(".out") or self.file.endswith(".log"):
            with open(self.file, 'r') as file:
                lines = file.readlines()

            if not any("normal termination of gaussian 16" in line.lower() for line in lines):
                return thermo_data

            keywords = {
                '&TRANSENERGY': 'Utrans',
                '&ROTENERGY': 'Urot',
                '&VIBENERGY': 'Uvib',
                '&QEL': 'Qel',
                '&QTRANS': 'Qtrans',
                '&QROT': 'Qrot',
                '&QVIB': 'Qvib',
                '&TEMPERATURE': 'T',
                '&PRESSURE': 'P',
                '&MASS': 'Mass',
                '&FINALEN': 'E',
                'zero-point correction': 'ZPE',
                'sum of electronic and thermal energies': 'U',
                'sum of electronic and thermal enthalpies': 'H',
                '&ENTROPYS': 'ST',
                'sum of electronic and thermal free energies': 'G',
                'thermal correction to gibbs free energy': 'corr'
                }        

            thermo_data.update({key: None for key in keywords.values()})

            for i in range(len(lines) - 1, -1, -1):
                line = lines[i]
                lower_line = line.lower()
                for keyword, key in keywords.items():
                    if thermo_data[key] is None and keyword in lower_line:
                        thermo_data[key] = float(line.split("=")[1].split()[0])

                if 'temperature' in lower_line and 'pressure' in lower_line:
                    thermo_data['T'] = float(line.split()[1])
                    thermo_data['P'] = float(line.split()[4])

                if 'molecular mass' in lower_line:
                    thermo_data['Mass'] = float(line.split(":")[1].split()[0])

                if 'e (thermal)' in lower_line:
                    for j in range(i + 2, i+7):
                        if lines[j].split()[0] == "Total":
                            thermo_data['ST'] = float(lines[j].split()[3]) * 10e-3  / constants.constants.ha_kcal

                        elif lines[j].split()[0] == "Electronic":
                            thermo_data['Qel'] = float(lines[j].split()[3]) * 10e-3  / constants.constants.ha_kcal

                        elif lines[j].split()[0] == "Translational":
                            thermo_data['Utrans'] = float(lines[j].split()[1])  / constants.constants.ha_kcal
                            thermo_data['Qtrans'] = float(lines[j].split()[3]) * 10e-3  / constants.constants.ha_kcal 

                        elif lines[j].split()[0] == "Rotational":
                            thermo_data['Urot'] = float(lines[j].split()[1])  / constants.constants.ha_kcal  
                            thermo_data['Qrot'] = float(lines[j].split()[3]) * 10e-3  / constants.constants.ha_kcal      

                        elif lines[j].split()[0] == "Vibrational":
                            thermo_data['Uvib'] = float(lines[j].split()[1])  / constants.constants.ha_kcal  
                            thermo_data['Qvib'] = float(lines[j].split()[3]) * 10e-3  / constants.constants.ha_kcal 

                              
                if thermo_data['E'] is None and 'scf done' in lower_line:
                    thermo_data['E'] = float(line.split(":")[1].split()[2])

            
            if thermo_data['ST'] is not None:
                thermo_data['ST'] = thermo_data['ST'] * thermo_data['T']

            return thermo_data                
