import constants
import os
import numpy as np
from functions import hessian

class orca:
    """Class to extract data from ORCA output files."""
    def __init__(self, calc, file):
        self.calc = calc
        self.file = file

        with open(self.file, 'r') as text:
            self.lines = text.readlines()

    def e_spin(self):
        "extract just the single point energy and the spin"
        e_spin = {}

        if self.file.endswith(".property.txt"):

            if not any("NORMAL TERMINATION" in line.upper() for line in self.lines):
                return None, None

            keywords = {
                '&FINALENERGY': 'E',
                '&VDW': 'VDW',
                'GCP_ENERGY': 'gCP',
                "&ELENERGY": 'El_e',
                "&FINALEN": 'E_therm',
                '&MULT': 'Spin',
                }
        
            e_spin.update({key: 0 for key in keywords.values()})

            for i in range(len(self.lines) - 1, -1, -1):
                line = self.lines[i]
                upper_line = line.upper()
                for keyword, key in keywords.items():
                    if e_spin[key] is 0 and keyword in upper_line:
                        e_spin[key] = float(line.split()[3])   

            if e_spin['E'] == 0:
                if e_spin['E_therm'] != 0:
                    e_spin['E'] = e_spin['E_therm']
                else:
                    e_spin['E'] = e_spin['El_e'] + e_spin['VDW'] + e_spin['gCP'] 

        elif self.file.endswith(".out"):

            if not any("ORCA TERMINATED NORMALLY"in line.upper().replace("*", "") for line in self.lines):
                return None, None
                   
            if any("the optimization did not converge"in line.lower() for line in self.lines):
                return None, None

            e_spin = {'E':0,'Spin':0}

            for i in range(len(self.lines) - 1, -1, -1):
                line = self.lines[i]
                
                if e_spin['E'] == 0 and 'FINAL SINGLE POINT ENERGY' in line:
                    e_spin['E'] = float(line.split()[4])

                if e_spin['Spin'] == 0 and 'Multiplicity' in line:
                    e_spin['Spin'] = float(line.split(':')[1])       

        return e_spin['E'], e_spin['Spin'] 

    def energies(self):
        """Extract relevant thermochemical data from the ORCA output file."""

        thermo_data = {}
        vdw = None
        gcp = None

        if self.file.endswith(".property.txt"):

            if not any("NORMAL TERMINATION" in line.upper() for line in self.lines):
                return thermo_data
            
            keywords = {
                '&TOTALMASS': 'Mass',
                '&PRESSURE': 'P',
                '&TEMPERATURE': 'T',
                '&MULT': 'Spin',
                '&TRANSENERGY': 'Utrans',
                '&ROTENERGY': 'Urot',
                '&VIBENERGY': 'Uvib',
                '&QEL': 'Qel',
                '&QTRANS': 'Qtrans',
                '&QROT': 'Qrot',
                '&QVIB': 'Qvib',
                '&FINALENERGY': 'E',
                '&ZPE': 'ZPE',
                '&INNERENERGYU': 'U',
                '&ENTHALPYH': 'H',
                '&ENTROPYS': 'ST',
                '&FREEENERGYG': 'G',
                '&corr': 'corr'
                }
        
            thermo_data.update({key: None for key in keywords.values()})

            for i in range(len(self.lines) - 1, -1, -1):
                line = self.lines[i]
                upper_line = line.upper()
                for keyword, key in keywords.items():
                    if thermo_data[key] is None and keyword in upper_line:
                        thermo_data[key] = float(line.split()[3])

                if "&NUMOFFREQS" in upper_line:
                    vib_count = int(line.split()[3])

                    for j in range(i + 5, i + 5 + vib_count):
                        if float(self.lines[j].split()[1]) < 0:
                            thermo_data['im_freq'] = float(self.lines[j].split()[1])
                            break

                if "&VDW" in upper_line:
                    vdw = float(line.split()[3])

                if "&GCP_ENERGY" in upper_line:
                    gcp = float(line.split()[3])

            if thermo_data['E'] is None:
                for i in range(len(self.lines) - 1, -1, -1):
                    line = self.lines[i]
                    upper_line = line.upper()
                    if "&ELENERGY" in upper_line:
                        thermo_data['E'] = float(line.split()[3])
                        break
                    elif "&FINALEN" in upper_line:
                        therm_data['E'] = float(line.split()[3])
                        if vdw is not None:
                            thermo_data['E'] += vdw
                        if gcp is not None:
                            thermo_data['E'] += gcp             

            if thermo_data['G'] is not None and thermo_data['E'] is not None:
                thermo_data['corr'] = thermo_data['G'] - thermo_data['E']
            else:
                thermo_data['corr'] = None

            return thermo_data


        elif self.file.endswith(".out"):

            if not any("ORCA TERMINATED NORMALLY"in line.upper().replace("*", "") for line in self.lines):
                return thermo_data
                   
            if any("the optimization did not converge"in line.lower() for line in self.lines):
                return thermo_data
            
            thermo_data.update({'vdw':None})
            thermo_data.update({'Spin':None})

            keywords = {
                'totalmass...': 'Mass',
                'pressure...': 'P',
                'temperature...': 'T',
                '&MULT': 'Spin',
                'thermaltranslationalcorrection...': 'Utrans',
                'thermalrotationalcorrection...': 'Urot',
                'thermalvibrationalcorrection...': 'Uvib',
                'Electronicentropy...': 'Qel',
                'translationalentropy...': 'Qtrans',
                'rotationalentropy...': 'Qrot',
                'vibrationalentropy...': 'Qvib',
                'electronicenergy...': 'E',
                'zeropointenergy...': 'ZPE',
                'totalthermalenergy...': 'U',
                'totalenthalpy...': 'H',
                'totalentropycorrection...': 'ST',
                'finalgibbsfreeenergy...': 'G',
                'g-e(el)...': 'corr'
                }
                
            thermo_data.update({key: None for key in keywords.values()})
            
            for i in range(len(self.lines) - 1, -1, -1):
                line = self.lines[i]
                lower_line = line.lower().replace(" ","")
                for keyword, key in keywords.items():
                    if thermo_data[key] is None and keyword in lower_line:

                        thermo_data[key] = float(line.split('...')[1].split()[0])
                if thermo_data['Spin'] is None and 'Multiplicity' in line:
                    thermo_data['Spin'] = float(line.split(':')[1])

                if thermo_data['E'] is None and 'FINAL SINGLE POINT ENERGY' in line:
                    thermo_data['E'] = float(line.split()[4])

                if 'imaginary mode' in line.replace('*',''):
                    thermo_data['im_freq'] = float(line.split()[1])

            return thermo_data

    def xyz(self, data=False):
        """Extract relevant xyz data from the ORCA output file."""
        energy = None 
        vdw = None 
        xyz_dict = {}

        if self.file.endswith('.property.txt'): 
            e_calc = None
            energy = None
            vdw = None

            if not any("NORMAL TERMINATION" in line.upper() for line in self.lines):
                return None

            for i in range(len(self.lines) - 1, -1, -1):
                lower_line = self.lines[i].lower()

                if '&natoms' in lower_line:
                    num_atoms = int(lower_line.split()[3])

                if '&cartesiancoordinates' in lower_line and not xyz_dict:
                    units = lower_line.split()[5].replace("\"", "").replace("]", "")
                    
                    for j, k in enumerate(range(i + 1, i + 1 + num_atoms)):
                        line_split = self.lines[k].split()
                        atom = line_split[0]

                        if units == 'bohr':
                            xyz = (float(line_split[1]) * constants.constants.bohr_A, float(line_split[2]) * constants.constants.bohr_A, float(line_split[3]) * constants.constants.bohr_A) 

                        else:
                            xyz = (float(line_split[1]), float(line_split[2]), float(line_split[3]))
                        
                        xyz_dict[f"{atom}_{j}"] = xyz

        elif self.file.endswith('.out'):    

            if not any("ORCA TERMINATED NORMALLY"in line.upper().replace("*", "") for line in self.lines):
                return None
                
            if any("the optimization did not converge"in line.lower() for line in self.lines):
                return None

            for i in range(len(self.lines) - 1, -1, -1):
                upper_line = self.lines[i].upper()

                if 'CARTESIAN COORDINATES (ANGSTROEM)' in upper_line and not xyz_dict:
                    for j, k in enumerate(range(i + 2, len(self.lines) - 1)):
                        
                        line_split = self.lines[k].split()

                        if not line_split:
                            break

                        atom = line_split[0]

                        xyz = (float(line_split[1]), float(line_split[2]), float(line_split[3]))

                        xyz_dict[f"{atom}_{j}"] = xyz
                
        return xyz_dict 
        
    def vib(self): 
        """Extract relevant vib data from the ORCA output file."""
        vib_dict = {}
        freq_found = False

        if self.file.endswith(".property.txt"):

            if not any("NORMAL TERMINATION" in line.upper() for line in self.lines):
                return vib_dict

            for i in range(len(self.lines) - 1, -1, -1):

                if freq_found == True:
                    break

                upper_line = self.lines[i].upper()   

                if "&FREQ " in upper_line:
                    freq_found = True
                    freq_start_index = self.lines.index(self.lines[i]) + 3 
                    for freq_line in self.lines[freq_start_index:]:
                        freq_line_upper = freq_line.upper()
                        if freq_line_upper.strip().startswith("&ZPE"): 
                            break
                        freq = float(freq_line.split()[1])

                        index = len(vib_dict) 
                        vib_dict[f"v_{index}"] = freq
     
            return vib_dict 
 

        elif self.file.endswith(".out"):

            if not any("ORCA TERMINATED NORMALLY"in line.upper().replace("*", "") for line in self.lines):
                return vib_dict
                   
            if any("the optimization did not converge"in line.lower() for line in self.lines):
                return vib_dict
            
            for i in range(len(self.lines) - 1, -1, -1):

                if freq_found == True:
                    break

                upper_line = self.lines[i].upper()   

                if "VIBRATIONAL FREQUENCIES" in upper_line:
                    freq_found = True
                    freq_start_index = self.lines.index(self.lines[i]) + 5 
                    for freq_line in self.lines[freq_start_index:]:
                        freq_line_upper = freq_line.upper()
                        if not freq_line_upper.strip(): 
                            break
                        
                        freq = float(freq_line.split()[1])

                        index = len(vib_dict) 
                        vib_dict[f"v_{index}"] = freq
     
            return vib_dict                    
        
        elif self.file.endswith(".hess"):
            nm_m_weight = {}


            for i in range(len(self.lines) - 1, -1, -1):

                lower_line = self.lines[i].lower()   

                if "$normal_modes" in lower_line:
                    dim = int(self.lines[i+1].split()[0])

                    for n in range(dim):
                        nm_m_weight[str(n)] = []

                    for j in range(i + 2, len(self.lines) - 1, dim + 1):  
                        if '#' in self.lines[j] or not self.lines[j].strip():
                            break
                        
                        for k in range(j + 1, j + 1 + dim):
                            nm_all = self.lines[k].split()
                            
                            if nm_all:
                                idx = nm_all[0]
                                nm = nm_all[1:]
                                
                                nm_m_weight[str(idx)].extend(map(float, nm))

            nm_dict = hessian(xyz).mass_weight_norm(nm_m_weight)
                        
            return nm_dict

            

    def hessian(self):
        """Extract the raw hessian data from the ORCA output file."""         
        hess_dict = {}

        if self.file.endswith(".property.txt"):

            if not any("NORMAL TERMINATION" in line.upper() for line in self.lines):
                return hess_dict      

            for i in range(len(self.lines) - 1, -1, -1):

                upper_line = self.lines[i].upper()   

                if "&HESSIAN" in upper_line:
                    dim = int(self.lines[i].split()[4].replace('(','').split(',')[0])

                    for n in range(dim):
                        hess_dict[str(n)] = []

                    for j in range(i + 1, len(self.lines) - 1, dim + 2):  
                        if '&modes' in self.lines[j].lower():
                            break
                        
                        for k in range(j + 1, j + 2 + dim):
                            hess_all = self.lines[k].split()
                            
                            if hess_all:
                                idx = hess_all[0]
                                hess = hess_all[1:]
                                
                                hess_dict[str(idx)].extend(map(float, hess))
                        
            return hess_dict
                           
        elif self.file.endswith(".hess"):
            hess_dict = {}

            for i in range(len(self.lines) - 1, -1, -1):

                lower_line = self.lines[i].lower()   

                if "$hessian" in lower_line:
                    dim = int(self.lines[i+1].split()[0])

                    for n in range(dim):
                        hess_dict[str(n)] = []

                    for j in range(i + 2, len(self.lines) - 1, dim + 1):  
                        if '#' in self.lines[j] or not self.lines[j].strip():
                            break
                        
                        for k in range(j + 1, j + 1 + dim):
                            nm_all = self.lines[k].split()
                            
                            if nm_all:
                                idx = nm_all[0]
                                nm = nm_all[1:]
                                
                                hess_dict[str(idx)].extend(map(float, nm))
                        
            return hess_dict

class g16:
    """Class to extract data from G16 output files."""
    def __init__(self, calc, file):
        self.calc = calc
        self.file = file
        
        with open(self.file, 'r') as text:
            self.lines = text.readlines()

    def e_spin(self):
         if self.file.endswith(".out") or self.file.endswith(".log"):
            
            if not any("normal termination of gaussian 16" in line.lower() for line in self.lines):
                return None  

            e_spin = {'E':None,'Spin':None}

            for i in range(len(self.lines) - 1, -1, -1):
                line = self.lines[i]
                lower_line = line.lower()

                if e_spin['E'] is None and 'scf done' in lower_line:
                    e_spin['E'] = float(line.split(":")[1].split()[2])

                if e_spin['Spin'] is None and 'multiplicity' in lower_line:
                    e_spin['Spin'] = float(line.split()[5])

            return e_spin['E'], e_spin['Spin']


    def energies(self):
        """Extract relevant data from the G16 output file."""

        thermo_data = {}

        if self.file.endswith(".out") or self.file.endswith(".log"):

            if not any("normal termination of gaussian 16" in line.lower() for line in self.lines):
                return thermo_data

            keywords = {
                '&MASS': 'Mass',
                '&PRESSURE': 'P',
                '&TEMPERATURE': 'T',
                '&Mult': 'Spin',
                '&TRANSENERGY': 'Utrans',
                '&ROTENERGY': 'Urot',
                '&VIBENERGY': 'Uvib',
                '&QEL': 'Qel',
                '&QTRANS': 'Qtrans',
                '&QROT': 'Qrot',
                '&QVIB': 'Qvib',
                '&FINALEN': 'E',
                'zero-point correction': 'ZPE',
                'sum of electronic and thermal energies': 'U',
                'sum of electronic and thermal enthalpies': 'H',
                '&ENTROPYS': 'ST',
                'sum of electronic and thermal free energies': 'G',
                'thermal correction to gibbs free energy': 'corr'
                }        

            thermo_data.update({key: None for key in keywords.values()})

            for i in range(len(self.lines) - 1, -1, -1):
                line = self.lines[i]
                lower_line = line.lower()
                for keyword, key in keywords.items():
                    if thermo_data[key] is None and keyword in lower_line:
                        thermo_data[key] = float(line.split("=")[1].split()[0])

                if 'molecular mass' in lower_line:
                    thermo_data['Mass'] = float(line.split(":")[1].split()[0])

                if 'temperature' in lower_line and 'pressure' in lower_line:
                    thermo_data['T'] = float(line.split()[1])
                    thermo_data['P'] = float(line.split()[4])

                if 'multiplicity' in lower_line:
                    thermod_data['Spin'] = float(line.split()[5])

                if 'e (thermal)' in lower_line:
                    for j in range(i + 2, i+7):
                        if self.lines[j].split()[0] == "Total":
                            thermo_data['ST'] = float(self.lines[j].split()[3]) * 10e-3  / constants.constants.ha_kcal

                        elif self.lines[j].split()[0] == "Electronic":
                            thermo_data['Qel'] = float(self.lines[j].split()[3]) * 10e-3  / constants.constants.ha_kcal

                        elif self.lines[j].split()[0] == "Translational":
                            thermo_data['Utrans'] = float(self.lines[j].split()[1])  / constants.constants.ha_kcal
                            thermo_data['Qtrans'] = float(self.lines[j].split()[3]) * 10e-3  / constants.constants.ha_kcal 

                        elif self.lines[j].split()[0] == "Rotational":
                            thermo_data['Urot'] = float(self.lines[j].split()[1])  / constants.constants.ha_kcal  
                            thermo_data['Qrot'] = float(self.lines[j].split()[3]) * 10e-3  / constants.constants.ha_kcal      

                        elif self.lines[j].split()[0] == "Vibrational":
                            thermo_data['Uvib'] = float(self.lines[j].split()[1])  / constants.constants.ha_kcal  
                            thermo_data['Qvib'] = float(self.lines[j].split()[3]) * 10e-3  / constants.constants.ha_kcal 

                              
                if thermo_data['E'] is None and 'scf done' in lower_line:
                    thermo_data['E'] = float(line.split(":")[1].split()[2])

            if thermo_data['ST'] is not None:
                thermo_data['ST'] = thermo_data['ST'] * thermo_data['T']

            return thermo_data   

    def xyz(self):
        """Extract relevant xyz data from the G16 output file."""             
        energy = None  
        xyz_dict = {}

        if self.file.endswith(".out") or self.file.endswith(".log"):

            if not any("normal termination of gaussian 16" in line.lower() for line in self.lines):
                return energy, xyz_dict

            for i in range(len(self.lines) - 1, -1, -1):
                line = self.lines[i]
                lower_line = line.lower()

                if 'coordinates (angstroms)' in lower_line and not xyz_dict:
                    for j, k in enumerate(range(i + 3, len(self.lines) - 1)):
                        
                        line_split = self.lines[k].split()

                        if "---" in self.lines[k]:
                            break
                        
                        atom_idx = int(line_split[1])
                        atom = list(constants.constants.atom_map.keys())[atom_idx-1]

                        xyz = (float(line_split[3]), float(line_split[4]), float(line_split[5]))

                        xyz_dict[f"{atom}_{j}"] = xyz

                if energy is None and 'scf done' in lower_line:
                    energy = float(line.split(":")[1].split()[2])
                
            return xyz_dict

    def vib(self):
        vib_dict = {}

        if self.file.endswith(".out") or self.file.endswith(".log"):

            if not any("normal termination of gaussian 16" in line.lower() for line in self.lines):
                return vib_dict
        
            for i in range(len(self.lines) - 1):
                  
                if self.lines[i].startswith(' Frequencies'):
                        freq_1 = float(self.lines[i].split()[2])

                        if len(self.lines[i].split()) == 3:
                            freq2 = None

                        else: 
                            freq_2 = float(self.lines[i].split()[3])

                        if len(self.lines[i].split()) == 5:
                            freq_3 = float(self.lines[i].split()[4])
                        else:
                            freq3 = None

                        index = len(vib_dict) 
                        vib_dict[f"v_{index}"] = freq_1
                        vib_dict[f"v_{index+1}"] = freq_2

                        if freq_3:
                            vib_dict[f"v_{index+2}"] = freq_3  

            return vib_dict                  

    def hessian(self):
        hess_dict = {}
        
        if self.file.endswith(".out") or self.file.endswith(".log"):

            if not any("normal termination of gaussian 16" in line.lower() for line in self.lines):
                return hess_dict
            data = ''
            for i in range(len(self.lines) - 1):
                if 'NImag=' in self.lines[i]:
                    for j in range(i, len(self.lines) - 1):
                        if ' The archive entry for this job was punched.' in self.lines[j]:
                            break
                        
                        data += self.lines[j].strip()
            
            hess_diag = data.split('NImag=')[1].split('\\')[2].split(',')
            h_elements = [float(x) for x in hess_diag]
            n = int((-1 + np.sqrt(1 + 8 * len(h_elements))) / 2)
            
            for k in range(n):
                hess_dict[k] = [0.0] * n
            
            idx = 0
            for l in range(n):
                for m in range(l + 1):
                    val = h_elements[idx]
                    hess_dict[l][m] = val 
                    hess_dict[m][l] = val 
                    idx += 1

            return hess_dict


                



