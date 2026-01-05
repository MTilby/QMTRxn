class thermal:
    def __init__(self, vib_data, inertia, **kwargs):
        self.vib_data = vib_data
        self.inertia = inertia
        self.vib_list = list(self.vib_data.values())
        self.vib_list = [v for v in self.vib_list if v > 0]
        self.temperature = float(kwargs.get('temperature', 298.15))
        self.atm = kwargs.get('atm')
        self.conc = kwargs.get('conc')
        self.method = kwargs.get('entropy_method', 'Quasi')
        self.freqH = kwargs.get('v0H', 100)
        self.alphaH = kwargs.get('alphaH', 4)
        self.freqS = kwargs.get('v0S', 100) 
        self.alphaS = kwargs.get('alphaS', 4) 

    def dampening(self, freq, alpha):
        damper = [1 / (1 + (freq / x) ** alpha) for x in self.vib_list]
        return damper
                    
    def thermal_trans(self):
        utrans = 1.5 * (constants.kb / constants.ha_j) * self.temperature 
        return utrans
    
    def thermal_rot(self, linear):
        if linear == True:
            urot = (constants.R / constants.ha_j / constants.na) * self.temperature 
        else:
            urot = 1.5 * (constants.R / constants.ha_j / constants.na) * self.temperature 
        
        return urot
    
    def thermal_vib(self):
        vib_unit = [constants.h * constants.c / constants.ha_j * x for x in self.vib_list]
    
        ZPE = 0.5 * np.sum(vib_unit)

        thermal_vib = [x / (math.exp(x * constants.ha_j / (constants.kb * self.temperature)) - 1) for x in vib_unit]
        uvib = np.sum(thermal_vib)

        return ZPE, uvib

    def quasi_u(self):
        vib_unit = [constants.h * constants.c / constants.ha_j * x for x in self.vib_list]
        thermal_vib = []
        damper = self.dampening(self.freqH, self.alphaH)
        for i, dampening in enumerate(damper):
            thermal_vib.append(dampening * vib_unit[i] / (math.exp(vib_unit[i] * constants.ha_j / (constants.kb * self.temperature)) - 1) + (1- dampening) * 0.5 * self.temperature * constants.kb / constants.ha_j)
        return np.sum(thermal_vib)

    def U_calc(self, E, linear):
        utrans = thermal_outputs.thermal_trans(self)
        urot = thermal_outputs.thermal_rot(self, linear)
        ZPE, uvib = thermal_outputs.thermal_vib(self)
        U = E + utrans + urot + ZPE + uvib
        return ZPE, U
    
    def enthalpy(self,U):
        H = U + constants.kb * self.temperature / constants.ha_j
        return H

    def Q_trans(self, mass):
        if self.atm is not None:
            vol_com = math.log(constants.kb * self.temperature / (float(self.atm) * constants.atm_pa))

        else:
            vol_com = math.log(1 / (constants.na * self.conc * 1e3))

        loggit = 1.5 * math.log(2 * math.pi * mass * constants.amu_kg * constants.kb * self.temperature / constants.h ** 2)

        strans = constants.kb * (2.5 + loggit + vol_com) / constants.ha_j * self.temperature

        return strans
    
    def Q_rot(self, symno, linear=False):  
        if linear == False:
            rotemp = [constants.h ** 2 / (8 * math.pi ** 2 * constants.kb * i * constants.ang_m ** 2 * constants.amu_kg) for i in self.inertia]
            qrot = (math.pi * self.temperature ** 3 / np.prod(rotemp)) ** 0.5 / symno 
            srot = constants.kb * (np.log(qrot) + 1.5) / constants.ha_j * self.temperature
        else:
            inerita = next(i for i in self.inertia if i > 0)
            rotemp = constants.h ** 2 / (8 * math.pi ** 2 * constants.kb * inerita * constants.ang_m ** 2 * constants.amu_kg) 
            qrot = self.temperature / (symno * rotemp)  
            srot = constants.kb * (np.log(qrot) + 1) / constants.ha_j * self.temperature

        return srot
            
    def Q_vib(self):
        if self.method == "Truhlar":
            vib_list = [self.freqS if v < self.freqS else v for v in self.vib_list]
        else:
            vib_list = self.vib_list
        
        theta = [constants.h * constants.c * x / (constants.kb * self.temperature) for x in vib_list]

        qvib = [x / (math.exp(x) - 1) - math.log(1 - math.exp(-x)) for x in theta]
        
        if self.method == "Quasi":
            qvib = self.quasi_s(qvib)
        
        svib = constants.kb * self.temperature * np.sum(qvib) / constants.ha_j

        return svib

    def quasi_s(self, RRHO):
        inert_unit = [i * constants.ang_m ** 2 * constants.amu_kg for i in self.inertia]
        avg_inertia = np.sum(inert_unit) / 3

        mu = [constants.h / (8 * math.pi **2 * x * constants.c) for x in self.vib_list]
        muI = [avg_inertia * x / (avg_inertia + x) for x in mu]

        q_fr = [0.5 + 0.5 * math.log(8 * math.pi ** 3 * x * constants.kb * self.temperature / constants.h **2) for x in muI]

        damper = self.dampening(self.freqS, self.alphaS)

        q_vib = []
        for i, dampening in enumerate(damper):
            q_vib.append(dampening * RRHO[i] + (1 - dampening) * q_fr[i])

        return q_vib

class xyz_functions:
    def __init__(self, XYZ):
        self.XYZ = XYZ

    def bond_vector(self, atom1, atom2):
        bv = np.array(atom1, dtype=float) - np.array(atom2, dtype=float)
        norm_bv = bv / np.linalg.norm(bv)
        
        return norm_bv

    def linear(self):
        if len(self.XYZ) < 3:
            return True
        
        else:
            atom_pairs = list(combinations(self.XYZ, 2))
            first_pair = atom_pairs[0]
            bv_first = self.bond_vector(self.XYZ[first_pair[0]], self.XYZ[first_pair[1]])
            
            for atom in self.XYZ:
                if atom != first_pair[0] and atom != first_pair[1]:
                    bv_other = self.bond_vector(self.XYZ[first_pair[0]], self.XYZ[atom])
                    cross_product = np.dot(bv_first, bv_other)
                    if cross_product == float(1):
                        return True
                    
            return False

 ### might need this if I want to introduce isotopes late stage ###

    def mass_calc(self):
        for atom in self.XYZ:
            print(atom)


    def inertia(self):
        mass = [constants.atom_mapping.get(''.join(filter(str.isalpha, atom))) for atom in self.XYZ]
        coords = np.array([self.XYZ[atom] for atom in self.XYZ])

        total_mass = np.sum(mass)
        cent = np.dot(mass, coords) / total_mass

        cent_coord = coords - cent
        inertia_matrix = sum(
            m * (np.identity(3) * np.dot(xyz, xyz) - np.outer(xyz, xyz))
            for m, xyz in zip(mass, cent_coord)
        )
        inertia, _ = np.linalg.eigh(inertia_matrix)

        return inertia

    def symm_no_calc(self):
        elements = [''.join(filter(str.isalpha, atom)) for atom in self.XYZ]
        coords = np.array([self.XYZ[atom] for atom in self.XYZ]) 

        pg = PointGroup(coords, elements)
        
        symno = constants.pg_mapping.get(pg.get_point_group())

        return symno