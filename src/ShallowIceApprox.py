#!/usr/bin/env python
# coding: utf-8

# # Shallow Ice Doodle Zone
# Area to experiment on the ice flow.

# In[ ]:

#get_ipython().magic(u'qtconsole')


# In[39]:

import dolfin as dfn
import numpy as np

# Constants Zone
n = 3.0             # Flow-law exponent
A = 10.0e-16        # 1 / (Pa^3 a)
spy = 31556926.0    # s / a
g = 9.81 * spy            # m / s**2

rho_i = 911.0       # kg / m**3
n = 3.0             # - dimensionless
a_dot = 0.3         # phi is a calculated value on page 4 of the paper.
k = 2.1             # W / (m^1 K^1)
c_p = 2009.0        # J/(kg * K)
T_0 = 273.15        # K
beta = 8.7e-4       # 1/(Km)
G = 42.0            # mW / M^2


T = 50000.0         # Maximum number of years to run the simulation
dt = 100.0         # timestep in years. 

# T is end time
# t is current time


# ## Diffusion Accumulation
# Implement the diffusion of the accumulation of the ice. 

# In[40]:

def dH_boundary(h, on_boundary):
    return on_boundary


mesh = dfn.RectangleMesh(0, 0, 1500000, 1500000, 31, 31)

FS = dfn.FunctionSpace(mesh, 'CG', 1)

H_n = dfn.Function(FS)
H_o = dfn.Function(FS)
dH = dfn.TrialFunction(FS)
w = dfn.TestFunction(FS) 
t = 0.0

def D(H):
    return ((2.0 * A * (rho_i * g)**n) / (n + 2.0)) * H**(n + 2.0) \
        * dfn.dot(dfn.grad(H), dfn.grad(H))**((n - 1.0) / 2.0)

F = (((H_n - H_o) / dt) * w) + (D(H_n) * \
    dfn.dot(dfn.grad(H_n), dfn.grad(w))) - (a_dot * w)
F = F * dfn.dx

F = dfn.action(F, H_n)

f = dfn.Constant(0)
bc = dfn.DirichletBC(FS, f, dH_boundary)


# Calculate the Jacobian Matrix
J = dfn.derivative(F, H_n, dH)
problem = dfn.NonlinearVariationalProblem(F, H_n, bc, J)
solver  = dfn.NonlinearVariationalSolver(problem)

slvr_prms = solver.parameters;
solver.parameters['nonlinear_solver'] = 'snes'
solver.parameters['snes_solver']['method']='virs'
solver.parameters['snes_solver']['maximum_iterations']=50
solver.parameters['snes_solver']['linear_solver']='mumps'
solver.parameters['snes_solver']['preconditioner'] = 'ilu'


icesim_file = dfn.File("icesim_out.pvd", "compressed")

# Bounds need to be functions
lb = dfn.interpolate(dfn.Expression('0'),FS)   # lower
ub = dfn.interpolate(dfn.Expression('0'),FS) # upper

t = dt
while t <= T:
    print "Time: %0.2f" % t
    t += dt
    solver.solve(lb, ub)
    H_o.vector()[:] = H_n.vector()
    icesim_file << (H_o.split()[0], t)


dfn.plot(H_n) 
dfn.interactive()


# In[10]:



