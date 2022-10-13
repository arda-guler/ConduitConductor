# import analysis module
from analysis import *

# set up a pipe mesh with uniform temperature distribution
pipe_mesh = generate_pipe_uniform_initial_conditions(30e-3, 15e-3, 25e-3, CCZ, 20, 3, 32, (273 + 25))

# set up boundary conditions
#set_T_at_end(pipe_mesh, 1000, True)
set_T_at_inner_corner(pipe_mesh, 1000, True)

# plot initial state
plot_temperatures(pipe_mesh)

# run the heat conduction simulation
conduct_heat(pipe_mesh, 1, 1e-03)

# plot final state
plot_temperatures(pipe_mesh)
