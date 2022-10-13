import math
M_PI = math.pi
M_2PI = M_PI * 2

import random
import matplotlib.pyplot as plt

from material import *
from meshing import *

SS = SS304L()
CCZ = CuCrZr()

# length, inner radius, outer radius, material, n_length, n_radial, n_circumferential
def generate_pipe_uniform_initial_conditions(L, r_in, r_out, mtl, n_l, n_r, n_c, T_init):
    new_sectors = []

    t = r_out - r_in

    step_l = (L/n_l)
    step_r = (t/n_r)
    step_c = (M_2PI/n_c)

    for cursor_l in range(n_l):
        current_l = cursor_l * step_l
        for cursor_r in range(n_r):
            current_r = r_in + cursor_r * step_r
            for cursor_c in range(n_c):
                current_theta = cursor_c * step_c
                
                new_index = [cursor_l, cursor_r, cursor_c]
                
                new_pos = [current_l,
                           current_r * math.cos(current_theta),
                           current_r * math.sin(current_theta)]

                new_r = current_r
                new_t = step_r
                new_theta = step_c
                new_h = step_l
                
                new_sector = sector(new_index, new_pos, new_r, new_t, new_theta, new_h, mtl, T_init)
                new_sectors.append(new_sector)

    new_mesh = mesh(n_l, n_r, n_c, new_sectors)
    return new_mesh

def set_T_at_end(workmesh, T, fix=False):
    sectors = workmesh.sectors
    for sector in sectors:
        if sector.cy_index[0] == 0:
            sector.T = T
            if fix:
                sector.fix_T()

def set_T_at_inner_corner(workmesh, T, fix=False):
    sectors = workmesh.sectors
    for sector in sectors:
        if sector.pos[1] <= 0 and sector.pos[2] <= 0 and sector.cy_index[1] == 0 and sector.cy_index[0] < workmesh.n_l/2:
            sector.T = T
            if fix:
                sector.fix_T()

def conduct_heat(workmesh, time, dt):
    n_t = int(time/dt)

    sectors = workmesh.sectors
    n_l = workmesh.n_l
    n_r = workmesh.n_r
    n_c = workmesh.n_c

    for current_step in range(n_t):
        step_transfers = []
        for current_sector in sectors:
            sector_transfers = []
            transfer_up = False
            transfer_in = False
            transfer_CW = False

            if not current_sector.cy_index[0] == 0:
                transfer_up = True

            if not current_sector.cy_index[1] == 0:
                transfer_in = True

            if not current_sector.cy_index[2] == 0:
                transfer_CW = True

            if transfer_up:
                top_sector_index = [current_sector.cy_index[0]-1, current_sector.cy_index[1], current_sector.cy_index[2]]
                top_sector = sectors[top_sector_index[0] * n_r * n_c + top_sector_index[1] * n_c + top_sector_index[2]]
                
                k = current_sector.mtl.get_thermal_conductivity(current_sector.T)
                A = current_sector.A_top
                dT = current_sector.T - top_sector.T
                dL = current_sector.pos[0] - top_sector.pos[0]
                
                if not dL == 0: # failsafe
                    Q = (k * A * dT / dL) * dt
                else:
                    Q = 0
                
                sector_transfers.append( [Q, current_sector.cy_index, top_sector_index] )

            if transfer_in:
                in_sector_index = [current_sector.cy_index[0], current_sector.cy_index[1]-1, current_sector.cy_index[2]]
                in_sector = sectors[in_sector_index[0] * n_r * n_c + in_sector_index[1] * n_c + in_sector_index[2]]

                k = current_sector.mtl.get_thermal_conductivity(current_sector.T)
                A = current_sector.A_in
                dT = current_sector.T - in_sector.T
                dL = ((current_sector.pos[1] - in_sector.pos[1])**2 + (current_sector.pos[2] - in_sector.pos[2])**2)**(0.5)
                
                if not dL == 0: # failsafe
                    Q = (k * A * dT / dL) * dt
                else:
                    Q = 0
                
                sector_transfers.append( [Q, current_sector.cy_index, in_sector_index] )

            if transfer_CW:
                CW_sector_index = [current_sector.cy_index[0], current_sector.cy_index[1], current_sector.cy_index[2]-1]
                CW_sector = sectors[CW_sector_index[0] * n_r * n_c + CW_sector_index[1] * n_c + CW_sector_index[2]]

                k = current_sector.mtl.get_thermal_conductivity(current_sector.T)
                A = current_sector.A_in
                dT = current_sector.T - CW_sector.T
                dL = ((current_sector.pos[1] - CW_sector.pos[1])**2 + (current_sector.pos[2] - CW_sector.pos[2])**2)**(0.5)
                
                if not dL == 0: # failsafe
                    Q = (k * A * dT / dL) * dt
                else:
                    Q = 0
                
                sector_transfers.append( [Q, current_sector.cy_index, CW_sector_index] )

            step_transfers.append(sector_transfers)

        for sec_transfers in step_transfers:
            for dir_transfer in sec_transfers:
                current_Q = dir_transfer[0]
                current_sector = sectors[dir_transfer[1][0] * n_r * n_c + dir_transfer[1][1] * n_c + dir_transfer[1][2]]
                other_sector = sectors[dir_transfer[2][0] * n_r * n_c + dir_transfer[2][1] * n_c + dir_transfer[2][2]]

                current_sector.remove_heat(current_Q)
                other_sector.add_heat(current_Q)

        if current_step % 100 == 0:
            print(str(round(current_step/n_t * 100, 2)) + "%")

    print("Done!")

def plot_temperatures(workmesh):
    T_max = 1000
    T_min = 200
    
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    sectors = workmesh.sectors

    xs = []
    ys = []
    zs = []
    colors = []
    for sector in sectors:
        xs.append(sector.pos[0]*1000) # convert to mm
        ys.append(sector.pos[1]*1000) # convert to mm
        zs.append(sector.pos[2]*1000) # convert to mm

        red = max(0, min(1, (sector.T - T_min)/(T_max - T_min)))
        blue = 1 - red
        colors.append((red, 0, blue))

    ax.scatter(xs, ys, zs, color=colors)
    ax.set_xlabel('X (mm)')
    ax.set_ylabel('Y (mm)')
    ax.set_zlabel('Z (mm)')
    plt.show()
