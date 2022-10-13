import math
M_PI = math.pi
M_2PI = M_PI * 2

class sector:
    def __init__(self, cy_index, pos, r, t, theta, h, mtl, T_init = (273+25)):
        self.cy_index = cy_index
        self.pos = pos
        self.r = r
        self.t = t
        self.theta = theta
        self.h = h

        self.r_out = r + t
        
        self.A_top = M_PI * (self.r_out**2 - self.r**2) * (self.theta/(M_2PI))
        self.A_in = r * theta * h
        self.A_out = self.r_out * theta * h
        self.A_side = t * h
        
        self.V = self.A_top * h

        self.mtl = mtl
        self.T = T_init

        self.fixed_T = False
        self.q_gen = 0

    def fix_T(self):
        self.fixed_T = True

    def unfix_T(self):
        self.fixed_T = False

    def get_heat_capacity(self):
        return self.V * self.mtl.get_density() * self.mtl.get_specific_heat(self.T)

    def remove_heat(self, Q):
        if not self.fixed_T:
            self.T -= Q/self.get_heat_capacity()

    def add_heat(self, Q):
        if not self.fixed_T:
            self.T += Q/self.get_heat_capacity()

class mesh:
    def __init__(self, n_l, n_r, n_c, sectors=[]):
        self.n_l = n_l
        self.n_r = n_r
        self.n_c = n_c
        self.sectors = sectors
