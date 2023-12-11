from mesh import Mesh
import numpy as np
import scipy.sparse
from scipy.sparse.linalg import lsmr
import importlib

m = Mesh("shell/slice.obj")
nboundary = sum(m.boundary)

model = "shell"
attr = importlib.import_module(model + ".attributes")

#for c in range(m.ncorners): # lift all vertices of all horizons
#    if attr.horizon_id[c]>=0:
#        height = (1+attr.horizon_id[c])/37.76 # arbitrarily chosen coeff to get a visually nice result
#        m.V[m.org(c)][2] = m.V[m.dst(c)][2] = height

"""average_horrizon = {}
for i in range(1,max(attr.horizon_id)+1):
    average, count = 0,0
    for elm in range(len(attr.horizon_id)):
        if elm == i:
        count +=1
        average += """

for dim in range(2): # solve for x first, then for y
    A = scipy.sparse.lil_matrix((nboundary + m.ncorners, m.nverts))
    b = [0] * A.shape[0]
    
    for row in range(m.ncorners):
        i = m.org(row)
        j = m.dst(row)
        A[row, j] =  1
        A[row, i] = -1

    for (i,v) in enumerate(m.V):
        if m.on_border(i):
            A[row, i] = 1 *100 # quadratic penalty to lock boundary vertices
            b[row] = v[dim] *100
            row += 1
           
    for c in range(m.ncorners): # per-half-edge discretization of the derivative
        A[m.nverts+c, m.org(c)] = -1
        A[m.nverts+c, m.dst(c)] = 1

    A = A.tocsr() # convert to compressed sparse row format for faster matrix-vector muliplications
    x = lsmr(A, b)[0] # call the least squares solver
    for i in range(m.nverts): # apply the computed flattening
        m.V[i][dim] = x[i]

m.write_vtk("output.vtk")
#print(m) # output the deformed mesh

