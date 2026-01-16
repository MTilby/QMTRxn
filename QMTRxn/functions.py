import constants
import math
import re
import numpy as np
from itertools import combinations
from pointgroup import PointGroup

class frequencies:
    def __init__(self, vib):
        self.vib = vib

    def correction(self, scale, invert, cutoff, lowfreq, map): 

        new_vib = {}
        neg_freq_found = False

        for v, value in self.vib.items():

            if value < 0:
                neg_freq_found = True

            if not scale:
                scale = 1.0

            if invert == 'All':
                freq = abs(value * float(scale))

            elif invert == 'nonTS':
                if self.freq < 0 and neg_freq_found == False:
                    freq = value * float(scale)
                    neg_freq_found = True

                else:
                    freq = abs(value * float(scale))
            
            else:
                freq = value * float(scale)

            if cutoff and abs(value) <= cutoff and value != 0:
                if lowfreq:
                    freq = lowfreq
                else:
                    freq = cutoff
            elif lowfreq and abs(value) <= lowfreq and value != 0:
                freq = lowfreq

            if map:
                index = int(''.join(filter(str.isdigit, v)))

                for v_new, value_new in map.items():
                    if type(v_new) == int:
                        map_index = v_new - 1
                    else:
                        map_index = int(''.join(filter(str.isdigit, v_new))) - 1

                    if map_index == index:
                        freq =  value_new

            new_vib[v] = freq
            
        return new_vib

class XYZ:
    def __init__(self, xyz):
        self.xyz = xyz

    def Mw(self):
        mass_list = []
        for atom_id in self.xyz:
            if '=' in atom_id:
                mass = float(re.search(r'M\(([\d\.]+)\)', atom_id).group(1))
                mass_list.extend([mass])
            else:
                mass = float(constants.constants.atom_map.get(''.join(filter(str.isalpha, atom_id))))
                mass_list.extend([mass])                

        Mw = np.sum(mass_list)
        return mass_list, Mw

    def isotope(self, change_atom, mass, amap): 
        new_xyz = {}
   
        if amap:
            for atom, value in self.xyz.items():
                index = int(''.join(filter(str.isdigit, atom)))
                element = ''.join(filter(str.isalpha, atom))
                atom_label = atom

                for map_atom in amap:
                    map_index = int(''.join(filter(str.isdigit, map_atom))) - 1
                    map_element = ''.join(filter(str.isalpha, map_atom))


                    if map_index == index:
                        if map_element != element and map_element:
                            print(f"Index of element '{map_atom}' in mapping does not match the index of the element {element}{index+1}, continuing anyway")
                            if not map_element:
                                map_element = element

                        atom_label = f"{map_element}=M({amap[map_atom]}){atom[len(element):]}"
                        break
                        
                new_xyz[atom_label] = value
                
            self.xyz = new_xyz

            return new_xyz

        elif change_atom and mass:
            for atom, value in self.xyz.items():
                element = ''.join(filter(str.isalpha, atom))
                if element == change_atom:
                    atom_label = f"{element}=M({mass}){atom[len(element):]}"
                else:
                    atom_label = atom
                new_xyz[atom_label] = value


        else:
            new_xyz = self.xyz
        
        return new_xyz
            

    def bond_vector(self, atom1, atom2):
        bv = np.array(atom1, dtype=float) - np.array(atom2, dtype=float)
        norm_bv = bv / np.linalg.norm(bv)

        return norm_bv
    
    def linear(self):
        if len(self.xyz) < 3:
            return True
        
        else:
            atom_pairs = list(combinations(self.xyz, 2))
            first_pair = atom_pairs[0]
            bv_first = self.bond_vector(self.xyz[first_pair[0]], self.xyz[first_pair[1]])
            
            for atom in self.xyz:
                if atom != first_pair[0] and atom != first_pair[1]:
                    bv_other = self.bond_vector(self.xyz[first_pair[0]], self.xyz[atom])
                    if np.linalg.norm(np.cross(bv_first, bv_other)) < 1e-6:
                        return True
                    
            return False
        
    def inertia(self, use_matrix=False):
        mass, total_mass = self.Mw()
        coords = np.array([self.xyz[atom] for atom in self.xyz])

        cent = np.dot(mass, coords) / total_mass

        cent_coord = coords - cent
        
        inertia_matrix = sum(
            m * (np.identity(3) * np.dot(xyz, xyz) - np.outer(xyz, xyz))
            for m, xyz in zip(mass, cent_coord))
        
        inertia, eigenvector = np.linalg.eigh(inertia_matrix) 

        if use_matrix == True:
            return eigenvector, cent_coord

        else:
            return inertia
    
    def symmno(self):
        elements = [''.join(filter(str.isalpha, atom)) for atom in self.xyz]
        coords = np.array([self.xyz[atom] for atom in self.xyz]) 

        pg = PointGroup(coords, elements)
        
        symno = constants.constants.pg_mapping.get(pg.get_point_group())

        return symno

class hessian:
    def __init__(self, xyz):
        self.xyz = XYZ(xyz) 

    def mass_extract(self):
        mass_list = []
        for atom_id in self.xyz.xyz:
            if '=' in atom_id:
                mass = float(re.search(r'M\(([\d\.]+)\)', atom_id).group(1))
                mass_list.extend([mass]*3)


            else:
                mass = float(constants.constants.atom_map.get(''.join(filter(str.isalpha, atom_id))))
                mass_list.extend([mass]*3)

        return np.array(mass_list)

    def frequency(self, hessian_matrix, negs):
        freqs = {}
        hess_array = np.array([v for v in hessian_matrix.values()])

        eigenvalues = np.linalg.eigvalsh(hess_array) * constants.constants.ha_j / (constants.constants.amu_kg * constants.constants.bohr_m **2)
        freq_squared = eigenvalues / (4 * (constants.constants.PI ** 2) * (constants.constants.c ** 2)) 
        frequencies = np.sign(freq_squared) * np.sqrt(np.abs(freq_squared)) 
 
        ### future work would project out these modes ###

        if self.xyz.linear() == True:
            im_freq = frequencies[:negs]
            low_freq = frequencies[negs:5+negs]
            low_freq = np.concatenate((low_freq, im_freq))
            frequencies[:5+negs] = low_freq
            frequencies[:5] = 0
        else:
            im_freq = frequencies[:negs]
            low_freq = frequencies[negs:6+negs]
            low_freq = np.concatenate((low_freq, im_freq))
            frequencies[:6+negs] = low_freq
            frequencies[:6] = 0

        for i, freq in enumerate(frequencies):
            freqs[f"v_{i}"] = float(freq)

        return freqs

    def mass_weight(self, matrix):
        modes = {}
        matrix_array = np.array([v for v in matrix.values()])

        mass_array = self.mass_extract()

        mass_matrix = np.diag(1/np.sqrt(mass_array))
           
        matrix_mw = mass_matrix @ matrix_array @ mass_matrix
         
        for i in matrix:
            modes[i] = matrix_mw[int(i)].tolist()

        return modes

class thermo:
    def __init__(self, xyz, vib_data, temperature):
        self.xyz = XYZ(xyz)
        vib_list = np.array(list(vib_data.values()))
        self.vib_list = vib_list[vib_list > 0]
        if temperature:
            self.temperature = float(temperature)
        else:
            self.temperature = 298.15

    def dampening(self, freq, alpha):
        damper = [1 / (1 + (freq / x) ** alpha) for x in self.vib_list]
        return damper

    def thermal_trans(self):
        utrans = 1.5 * (constants.constants.kb / constants.constants.ha_j) * self.temperature 
        return utrans    

    def thermal_rot(self):
        
        if self.xyz.linear() == True:
            urot = (constants.constants.R / constants.constants.ha_j / constants.constants.na) * self.temperature 
        else:
            urot = 1.5 * (constants.constants.R / constants.constants.ha_j / constants.constants.na) * self.temperature 
        
        return urot

    def thermal_vib(self):
        vib_unit = [constants.constants.h * constants.constants.c / constants.constants.ha_j * x for x in self.vib_list]

        ZPE = 0.5 * np.sum(vib_unit)

        thermal_vib = [x / (np.exp(x * constants.constants.ha_j / (constants.constants.kb * self.temperature)) - 1) for x in vib_unit]
        uvib = np.sum(thermal_vib).item()

        return ZPE, uvib

    def quasi_u(self, v0H, alphaH):
        vib_unit = [constants.constants.h * constants.constants.c / constants.constants.ha_j * x for x in self.vib_list]
        thermal_vib = []
        damper = self.dampening(v0H, alphaH)
        for i, dampening in enumerate(damper):
            thermal_vib.append(dampening * vib_unit[i] / (math.exp(vib_unit[i] * constants.constants.ha_j / (constants.constants.kb * self.temperature)) - 1) + (1- dampening) * 0.5 * self.temperature * constants.constants.kb / constants.constants.ha_j)
        return np.sum(thermal_vib)

    def Q_el(self, mult):
        sel = np.log(mult) * (constants.constants.kb / constants.constants.ha_j) * self.temperature 
        return sel

    def Q_trans(self, atm, conc):
        mass, total_mass = self.xyz.Mw()
        if atm:
            vol_com = np.log(constants.constants.kb * self.temperature / (atm * constants.constants.atm_pa))

        else:
            vol_com = np.log(1 / (constants.constants.na * conc * 1e3))

        loggit = 1.5 * np.log(2 * constants.constants.PI * total_mass * constants.constants.amu_kg * constants.constants.kb * self.temperature / constants.constants.h ** 2)

        strans = constants.constants.kb * (2.5 + loggit + vol_com) / constants.constants.ha_j * self.temperature

        return strans

    def Q_rot(self, symno):
        linear = self.xyz.linear()
        inertia = self.xyz.inertia()  
        if linear == False:
            rotemp = [constants.constants.h ** 2 / (8 * constants.constants.PI ** 2 * constants.constants.kb * i * constants.constants.ang_m ** 2 * constants.constants.amu_kg) for i in inertia]
            qrot = (constants.constants.PI * self.temperature ** 3 / np.prod(rotemp)) ** 0.5 / symno 
            srot = constants.constants.kb * (np.log(qrot) + 1.5) / constants.constants.ha_j * self.temperature
        else:
            inerita = next(i for i in inertia if i > 0)
            rotemp = constants.constants.h ** 2 / (8 * math.pi ** 2 * constants.constants.kb * inerita * constants.constants.ang_m ** 2 * constants.constants.amu_kg) 
            qrot = self.temperature / (symno * rotemp)  
            srot = constants.constants.kb * (np.log(qrot) + 1) / constants.constants.ha_j * self.temperature

        return srot

    def quasi_s(self, RRHO, freqS, alphaS):
        inertia = self.xyz.inertia()
        inert_unit = [i * constants.constants.ang_m ** 2 * constants.constants.amu_kg for i in inertia]
        avg_inertia = np.sum(inert_unit) / 3

        mu = [constants.constants.h / (8 * constants.constants.PI **2 * x * constants.constants.c) for x in self.vib_list]
        muI = [avg_inertia * x / (avg_inertia + x) for x in mu]

        q_fr = [0.5 + 0.5 * np.log(8 * constants.constants.PI ** 3 * x * constants.constants.kb * self.temperature / constants.constants.h **2) for x in muI]

        damper = self.dampening(freqS, alphaS)

        q_vib = []
        for i, dampening in enumerate(damper):
            q_vib.append(dampening * RRHO[i] + (1 - dampening) * q_fr[i])

        return q_vib

    def Q_vib(self, method, freqS, alphaS):
        if method == "Truhlar":
            vib_list = [freqS if v < freqS else v for v in self.vib_list]
        else:
            vib_list = self.vib_list
        
        theta = [constants.constants.h * constants.constants.c * x / (constants.constants.kb * self.temperature) for x in vib_list]

        qvib = [x / (np.exp(x) - 1) - np.log(1 - np.exp(-x)) for x in theta]
        
        if method == "Quasi":
            qvib = self.quasi_s(qvib, freqS, alphaS)
        
        svib = constants.constants.kb * self.temperature * np.sum(qvib) / constants.constants.ha_j

        return svib

    def thermo_data(self, energy, spin, quasiH, v0H, alphaH, atm, conc, symno, method, v0S, alphaS):
        thermo_data = {'Mass':self.xyz.Mw()[1].item()}

        if atm is None and conc is None:
            thermo_data['P'] = 1
            atm = 1
        elif atm:
            thermo_data['P'] = atm

        if atm is None and conc:
            thermo_data['conc'] = conc

        thermo_data['T'] = self.temperature
        thermo_data['Spin'] = spin


        thermo_data['Utrans'] = self.thermal_trans()
        thermo_data['Urot'] = self.thermal_rot()
        ZPE, thermo_data['Uvib'] = self.thermal_vib()

        if quasiH or v0H or alphaH:
            if not v0H:
                v0H = 100

            if not alphaH:
                alphaH = 4
            thermo_data['Uvib'] = self.quasi_u(v0H, alphaH)

            
        if not symno:
            symno = self.xyz.symmno()

        if method or v0S or alphaS:
            if not method:
                method = "Quasi"

            if not v0S:
                v0S = 100

            if not alphaS:
                alphaS = 4 
        
        thermo_data['Qel'] = self.Q_el(spin).item()
        thermo_data['Qtrans'] = self.Q_trans(atm, conc).item()
        thermo_data['Qrot'] = self.Q_rot(symno).item()
        thermo_data['Qvib'] = self.Q_vib(method, v0S, alphaS).item()

        thermo_data['E'] = energy
        thermo_data['ZPE'] = ZPE.item()
        thermo_data['U'] = thermo_data['E'] + thermo_data['Utrans'] + thermo_data['Urot'] + ZPE.item() +  thermo_data['Uvib']
        thermo_data['H'] = thermo_data['U'] + constants.constants.kb * self.temperature / constants.constants.ha_j
        thermo_data['ST'] = thermo_data['Qel'] + thermo_data['Qtrans'] + thermo_data['Qrot'] + thermo_data['Qvib']
        thermo_data['G'] = thermo_data['H'] - thermo_data['ST']
        thermo_data['corr'] = thermo_data['G'] - energy

        return thermo_data







class QMT:
    def __init__(self, Vr, Vp, mu):
        self.Vr = Vr
        self.Vp = Vp
        self.mu = mu

### check naming ###

    def D(self, V_max, F_s):
        Beta = (V_max ** 0.5 + (V_max - (self.Vp - self.Vr)) ** 0.5) ** 2
        Alpha = (Beta * F_s / (2 * V_max * (V_max - (self.Vp - self.Vr)))) ** 0.5
        D = 2 * constants.PI * abs((2 * self.mu * Beta - (Alpha/2)**2))**0.5 / Alpha
        return Alpha, D

    def A(self,E,ALPHA):
        return  2 * constants.PI * (2 * self.mu * E)**0.5 / ALPHA

    def B(self, E,V_p,V_r,ALPHA):
        return  2 * constants.PI * (2 * self.mu * (E - (V_p - V_r)))**0.5 / ALPHA

#Calculation of Transmission Probabilty of Eckart Potential
    def T(a,b,d):
        T = (math.cosh(a+b) - math.cosh(a-b))/(math.cosh(a+b) + math.cosh(d))
        return T