# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# # Shallow Ice Doodle Zone
# Area to experiment on the ice flow.

# <codecell>

# %qtconsole

# <markdowncell>

# Setup the variables as defined in _The EISMINT benchmarks for testing ice-sheet models_.

# <codecell>

import dolfin as dfn
import numpy as np

# Constants Zone
n = 3.0             # Flow-law exponent
A = 1.0e-16         # 1 / (Pa^3 a)
spy = 31556926.0    # s / a
g = 9.81            # m / s**2

rho_i = 911.0       # kg / m**3
n = 3.0             # - dimensionless
a_dot = 0.3         # phi is a calculated value on page 4 of the paper.
k = 2.1             # W / (m^1 K^1)
c_p = 2009.0        # J/(kg * K)
T_0 = 273.15        # K
beta = 8.7e-4       # 1/(Km)
G = 42.0            # mW / M^2


T = 250000     # Maximum number of years to run the simulation
dt = 5000.0        # timestep in years. 

# <codecell>

def D(h):
    a = 2 * A * (rho_i * g)**n / (n + 2)
    b = h**(n+2)
    c = dfn.dot(dfn.grad(h), dfn.grad(h))
    d = (n - 1) / 2
    return a*b*(c**d)

# <codecell>

mesh = dfn.RectangleMesh(0, 0, 1500000, 1500000, 30, 30)

# <codecell>

V = dfn.FunctionSpace(mesh, 'CG', 1)

# <codecell>

# Define the trial and test functions
dH = dfn.TrialFunction(V)
w = dfn.TestFunction(V)

# <codecell>

# Now define the initial conditions
H_o = dfn.Function(V)
H_n = dfn.Function(V)

# <markdowncell>

# Residual:
# \begin{align*}
#    0 &= \frac{\partial H}{\partial t} + \nabla \cdot D \nabla H + M
# \end{align*}
# 
# Conversion into weak form. 
# \begin{align*}
#   0 &= \int_\Omega\left[\frac{\partial H}{\partial t} * w + \cdot D \nabla H \cdot \nabla{w} - M * w\right]
# \end{align*}
# 
# Where $M = \dot a$ in this project. The first component assumes that $\dot a$ is a constant.

# <codecell>

a = (H_n - H_o) / dt * w
b = D(H_n) * dfn.dot(dfn.grad(H_n), dfn.grad(w))
c = a_dot * w
F = (a + b - c) * dfn.dx

# <codecell>

# F = dfn.action(F, H_coef)

# <markdowncell>

# Compute the Jacobian matrix. 

# <codecell>

J = dfn.derivative(F, H_n, dH)

# <markdowncell>

# Setup the boundaries to contain a constant value (we do not want any Von Neuman boundaries)

# <codecell>

def boundary(h, on_boundary):
    return on_boundary

bc = dfn.DirichletBC(V, dfn.Constant(0.0), boundary)

# <markdowncell>

# Setup the the problem so we can solve it. 

# <codecell>

problem = dfn.NonlinearVariationalProblem(F, H_n, bc, J)

# <markdowncell>

# Construct the solver

# <codecell>

solver = dfn.NonlinearVariationalSolver(problem)

# <markdowncell>

# Set the solver parameters

# <codecell>

solver.parameters['nonlinear_solver'] = 'snes'
solver.parameters['snes_solver']['method']='virs'
solver.parameters['snes_solver']['maximum_iterations']=50
solver.parameters['snes_solver']['linear_solver']='mumps'
solver.parameters['snes_solver']['preconditioner'] = 'ilu'

prm = solver.parameters
prm['newton_solver']['absolute_tolerance'] = 1E-8
prm['newton_solver']['relative_tolerance'] = 1E-7
prm['newton_solver']['maximum_iterations'] = 25
prm['newton_solver']['relaxation_parameter'] = 1.0

# <codecell>

t = 0
while t <= T:
    solver.solve()
    H_o.vector()[:] = H_n.vector()
    print t
    dfn.plot(H_n)
    t += dt

# <codecell>

dfn.plot(H_n)

# <codecell>

dfn.interactive()

# <codecell>


