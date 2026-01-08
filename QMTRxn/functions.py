import constants
import math
import numpy as np
from itertools import combinations
from pointgroup import PointGroup

class frequencies:
    def __init__(self, freq):
        self.freq = freq

    def correction(self, scale, invert, cutoff, lowfreq, neg_freq_found): 
        if not scale:
            scale = 1.0

        if invert == 'All':
            freq = abs(self.freq * float(scale))

        elif invert == 'nonTS':
            if self.freq < 0 and neg_freq_found == False:
                freq = self.freq * float(scale)
                neg_freq_found = True

            else:
                freq = abs(self.freq * float(scale))
        
        else:
            freq = self.freq * float(scale)

        if cutoff and abs(self.freq) <= cutoff and self.freq != 0:
            if lowfreq:
                freq = lowfreq
            else:
                freq = cutoff
        elif lowfreq and abs(self.freq) <= lowfreq and self.freq != 0:
            freq = lowfreq

        return freq, neg_freq_found

class XYZ:
    def __init__(self, xyz):
        self.xyz = xyz

    def isotope(self, change_atom, mass, map): ###recalc hessian
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

            self.xyz = new_xyz
        
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
            print(first_pair)
            bv_first = self.bond_vector(self.xyz[first_pair[0]], self.xyz[first_pair[1]])
            
            for atom in self.xyz:
                if atom != first_pair[0] and atom != first_pair[1]:
                    bv_other = self.bond_vector(self.xyz[first_pair[0]], self.xyz[atom])
                    cross_product = np.dot(bv_first, bv_other)
                    if cross_product == float(1):
                        return True
                    
            return False
        
    def inertia(self):
        mass = [constants.constants.atom_map.get(''.join(filter(str.isalpha, atom))) for atom in self.xyz]
        coords = np.array([self.xyz[atom] for atom in self.xyz])

        total_mass = np.sum(mass)
        print(total_mass)
        cent = np.dot(mass, coords) / total_mass

        cent_coord = coords - cent
        inertia_matrix = sum(
            m * (np.identity(3) * np.dot(xyz, xyz) - np.outer(xyz, xyz))
            for m, xyz in zip(mass, cent_coord)
        )
        inertia, _ = np.linalg.eigh(inertia_matrix)

        return inertia
    
    def symmno(self):
        elements = [''.join(filter(str.isalpha, atom)) for atom in self.xyz]
        coords = np.array([self.XYZ[atom] for atom in self.xyz]) 

        pg = PointGroup(coords, elements)
        
        symno = constants.pg_mapping.get(pg.get_point_group())

        return symno

class hessian:
    def __init__(self, xyz):
        self.xyz = xyz

    def mass_hess(self, matrix):
        print('gets mass w hessian')

    def hess(self, matrix):
        for atom_id in self.xyz:
            if '=' in atom_id:
                mass = int(''.join(filter(str.isdigit, atom_id.split('=')[1].split('_')[0])))
                print(mass)


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