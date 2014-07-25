__author__ = 'Kamil Koziara & Taiyeb Zahir'

import cProfile

import numpy
from utils import generate_pop, HexGrid, draw_hex_grid
from stops_ import Stops2

secretion = numpy.array([0])
reception = numpy.array([2])
receptors = numpy.array([-1])
bound=numpy.array([1,1,1,1])

base1=numpy.array([0,1,0,0])


trans_mat = numpy.array([[0,0,0,0], #Delta_ligand
                         [0.005,0,0,-1], #Delta
                         [-1,-1,0,0], #notch receptor
                         [0,0,0,0] #dummy
                         ])

init_pop = generate_pop([(2, base1)])
grid = HexGrid(1, 2, 1)

cell1, cell2, err, not_complete = 0, 0, 0, 0

m=1000 # number of simulations

def run():
    global cell1, cell2, err, not_complete
    x = Stops2(trans_mat, init_pop, grid.adj_mat, bound, secretion, reception, receptors, secr_amount=1, leak=0, max_con=1, max_dist=1.5, opencl=False)
    for i in range(1300):
        x.step()
        if (i+1)%1300 == 0:
            if (x.pop[0,2]==1 and x.pop[1,2]==0) :
                cell1+=1
            elif (x.pop[1,2]==1 and x.pop[0,2]==0) :
                cell2+=1
            elif (x.pop[0,2]==1 and x.pop[1,2]==1):
                err+=1
            else:
                not_complete+=1
            print x.pop
        
for i in range(m):
    print i
    cProfile.run("run()")

print "right cell=", cell1, "left cell =", cell2, "error=", err, "not complete=", not_complete
