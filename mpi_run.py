from mpi4py import MPI
import numpy as np
from Init.init_super_simulation import *

def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn

def main():
    
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nprocs = comm.Get_size()

    if rank == 0:
        data = [super_simulation,
                super_simulation,
                super_simulation,
                super_simulation]

        # determine the size of each sub-task
        ave, res = divmod(len(data), nprocs)
        counts = [ave + 1 if p < res else ave for p in range(nprocs)]

        # determine the starting and ending indices of each sub-task
        starts = [sum(counts[:p]) for p in range(nprocs)]
        ends = [sum(counts[:p+1]) for p in range(nprocs)]

        # converts data into a list of arrays 
        data = [data[starts[p]:ends[p]] for p in range(nprocs)]
    else:
        data = None

    data = comm.scatter(data, root=0)

    print('Process {} has data:'.format(rank), data)

    data[0](EMT(),'Cu')
    
if __name__ == "__main__":
    main()
