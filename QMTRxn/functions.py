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

    def isotope(self, change_atom, mass, map): 
        new_xyz = {}
   
        if map:
            for atom, value in self.xyz.items():
                index = int(''.join(filter(str.isdigit, atom)))
                element = ''.join(filter(str.isalpha, atom))
                atom_label = atom

                for map_atom in map:
                    map_index = int(''.join(filter(str.isdigit, map_atom))) - 1
                    map_element = ''.join(filter(str.isalpha, map_atom))


                    if map_index == index:
                        if map_element != element and map_element:
                            print(f"Index of element '{map_atom}' in mapping does not match the index of the element {element}{index+1}, continuing anyway")
                            if not map_element:
                                map_element = element

                        atom_label = f"{map_element}=M({map[map_atom]}){atom[len(element):]}"
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
        mass = [constants.constants.atom_map.get(''.join(filter(str.isalpha, atom))) for atom in self.xyz]
        coords = np.array([self.xyz[atom] for atom in self.xyz])

        total_mass = np.sum(mass)
        cent = np.dot(mass, coords) / total_mass

        cent_coord = coords - cent
        inertia_matrix = sum(
            m * (np.identity(3) * np.dot(xyz, xyz) - np.outer(xyz, xyz))
            for m, xyz in zip(mass, cent_coord))
        
        inertia, eigenvector = np.linalg.eigh(inertia_matrix) 

        if use_matrix == True:
            return inertia_matrix, eigenvector, cent_coord

        else:
            return inertia
    
    def symmno(self):
        elements = [''.join(filter(str.isalpha, atom)) for atom in self.xyz]
        coords = np.array([self.XYZ[atom] for atom in self.xyz]) 

        pg = PointGroup(coords, elements)
        
        symno = constants.pg_mapping.get(pg.get_point_group())

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

    def frequency(self, hessian_matrix, normal_modes):
        freqs = {}
        hess_array = np.array([v for v in hessian_matrix.values()])
        normal_array = np.array([v for v in normal_modes.values()])
        hess_nm =  normal_array.T @ hess_array @  normal_array

        eigenvalues = np.diag(hess_nm) * constants.constants.ha_j / (constants.constants.amu_kg * constants.constants.bohr_m **2)

        frequencies = np.sqrt((eigenvalues) / (4 * (constants.constants.PI ** 2) * (constants.constants.c ** 2)))  

        for i, freq in enumerate(frequencies):
            freqs[f"v_{i}"] = float(freq)

        return freqs

    def mass_weight(self, matrix):
        normal_modes = {}
        mass_array = self.mass_extract()
        matrix_array = np.array([v for v in matrix.values()]) 
        matrix_mw = matrix_array / np.sqrt(mass_array)
         
        norms = np.linalg.norm(matrix_mw, axis=0)
        norms[norms < 1e-10] = 1.0

        matrix_mw = matrix_mw / norms

        for i in matrix:
            normal_modes[i] = matrix_mw[int(i)].tolist()

        return normal_modes

class thermo:
    def __init__(self, vib_data, inertia, **kwargs):
        self.vib_data = vib_data
        self.inertia = inertia

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