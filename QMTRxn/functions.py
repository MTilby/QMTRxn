import constants
import math

class Eckart:
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
    