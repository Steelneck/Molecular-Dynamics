from mpi4py import MPI
from Init.init_super_simulation import *
import os
import numpy as np

os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'
os.environ['NUMEXPR_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'

def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

if rank == 0:
    jobs = np.array(['Cu', 'Ag', 'Pd', 'Ni'])
    job_array = np.array_split(jobs, size)
    print("We have", size, "processes.")
    for i in range(0, size):
        comm.isend(job_array[i], dest=i, tag=i)

my_jobs = comm.recv(source=0, tag=rank)
result = [super_simulation(EMT(), my_job) for my_job in my_jobs]
comm.isend(result, dest=0, tag=0)

if rank == 0:
    for i in range(0, size):
        result = comm.recv(source=i, tag=0)



# def test_calc(size):
#     A = np.ones((size*10, size*10))
#     A = A.dot(A) + np.diag([1]*size*10)
#     A = np.linalg.det(A)
#     print("Calc of", size, "gave", A)
#     return A

# def main():

#     comm = MPI.COMM_WORLD
#     rank = comm.Get_rank()
#     size = comm.Get_size()

#     if rank == 0:
#         jobs = np.arange(1,10)
#         job_array = np.array_split(jobs,size)
#         print("We have", size, "processes.")
#         for i in range(0, size):
#             comm.isend(job_array[i], dest=i, tag=i)

#     my_jobs = comm.recv(source=0, tag=rank)
#     result = [test_calc(my_job) for my_job in my_jobs]
#     comm.isend(result, dest=0, tag=0)

#     if rank == 0:
#         for i in range(0, size):
#             result = comm.recv(source=i, tag=0)

# if __name__ == "__main__":
#     main()

# def main():

#     comm = MPI.COMM_WORLD
#     rank = comm.Get_rank()
#     size = comm.Get_size()

#     if rank == 0:
#         print("We have", size, "processes.")
#         for i in range(0, size):
#             comm.isend("Hello "+str(i)+"!", dest=i, tag=i)

#     data = comm.recv(source=0, tag=rank)
#     msg = "Hello from "+str(rank)+"."
#     comm.isend(msg, dest=0, tag=0)

#     if rank == 0:
#         for i in range(0, size):
#             msg = comm.recv(source=i, tag=0)
#             print(msg)

# if __name__ == "__main__":
#     main()

# def main():

#     comm = MPI.COMM_WORLD
#     size = comm.Get_size()
#     rank = comm.Get_rank()

#     numDataPerRank = 2
#     data = None
#     if rank == 0:
#         data = np.linspace(1,size*numDataPerRank, numDataPerRank*size)

#     recvbuf = np.empty(numDataPerRank, dtype = 'd')
#     comm.Scatter(data, recvbuf, root=0)

#     print('Rank: ', rank, ', recvbuf received: ', recvbuf)

#     print(data)
#     print(recvbuf)

# if __name__ == "__main__":
#     main()

# def main():
    
#     comm = MPI.COMM_WORLD
#     rank = comm.Get_rank()
#     nprocs = comm.Get_size()

#     if rank == 0:
#         data = [super_simulation,
#                 super_simulation,
#                 super_simulation,
#                 super_simulation]

#         # determine the size of each sub-task
#         ave, res = divmod(len(data), nprocs)
#         counts = [ave + 1 if p < res else ave for p in range(nprocs)]

#         # determine the starting and ending indices of each sub-task
#         starts = [sum(counts[:p]) for p in range(nprocs)]
#         ends = [sum(counts[:p+1]) for p in range(nprocs)]

#         # converts data into a list of arrays 
#         data = [data[starts[p]:ends[p]] for p in range(nprocs)]
#     else:
#         data = None

#     data = comm.scatter(data, root=0)

#     print('Process {} has data:'.format(rank), data)

#     data[0](EMT(),'Cu')
    
# if __name__ == "__main__":
#     main()
