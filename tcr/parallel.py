"""
Functions for parallel processing in PyTCR
"""


import sys
import numpy as np


def domain_decompose(comm, rank, size, storm_list, verbose=False):
    """
    Decomposes the domain of storms across MPI ranks for parallel processing.

    Parameters:
    ----------
    - comm: MPI communicator object
    - rank: The rank of the current MPI process
    - size: The total number of MPI processes
    - storm: The list of storms to be distributed
    - verbose: A boolean flag to enable or disable verbose output

    Returns:
    --------
    - loc_storm_list: A list of storms assigned to the current MPI rank
    - loc_num_storms: The number of storms assigned to the current MPI rank
    """

    if rank == 0:
        print(f"Total number of storms: {len(storm_list)}")
        numpairs = np.shape(storm_list)[0]
        counts = np.arange(size, dtype=np.int32)
        displs = np.arange(size, dtype=np.int32)
        ave = int(numpairs / size)
        extra = numpairs % size
        offset = 0

        for i in range(0, size):
            loc_num_storms = ave if i < size - extra else ave + 1
            counts[i] = loc_num_storms

            if i == 0:
                loc_num_storms0 = loc_num_storms
                offset += loc_num_storms
                displs[i] = 0
            else:
                comm.send(offset, dest=i)
                comm.send(loc_num_storms, dest=i)
                offset += loc_num_storms
                displs[i] = displs[i - 1] + counts[i - 1]

            if verbose:     
                for j in range(offset - loc_num_storms, offset):
                    print(f"Rank: {i} - {storm_list[j]}")
                    sys.stdout.flush()

        offset = 0
        loc_num_storms = loc_num_storms0

    comm.Barrier()

    if rank != 0:  # workers
        offset = comm.recv(source=0)
        loc_num_storms = comm.recv(source=0)

    comm.Barrier()
    loc_storm_list = storm_list[offset: offset + loc_num_storms]
    return loc_storm_list, loc_num_storms
