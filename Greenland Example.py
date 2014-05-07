# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# # Greenland Example
# Run through Jesse's Greenland Example using FEniCS and Paraview

# <codecell>

import dolfin as dfn
import numpy as np

# <codecell>

# Load the data files
expth = "./example/Greenland_10km/"
mesh = dfn.Mesh(expth + "GIS_mesh_10km.xml")

FS = dfn.FunctionSpace(mesh, "CG", 1)

adot = dfn.Function(FS)
S_o  = dfn.Function(FS)
H_o  = dfn.Function(FS)
u_o  = dfn.Function(FS)

dfn.File(expth + "adot.xml") >> adot
dfn.File(expth + "S.xml") >> S_o
dfn.File(expth + "H.xml") >> H_o
dfn.File(expth + "vmag.xml") >> u_o

# <codecell>

#FEniCS freezes after closing the plot. 
dfn.plot(mesh)
dfn.plot(S_o)
dfn.plot(H_o)
dfn.plot(u_o)
dfn.interactive()

# <codecell>



