"""
Functions for parallel processing in PyTCR
"""


import sys
import numpy as np


def domain_decompose(comm, rank, size, data, verbose=False):
    """
    Decomposes the domain of storms across MPI ranks for parallel processing.

    Parameters:
    ----------
    - comm: MPI communicator object
    - rank: The rank of the current MPI process
    - size: The total number of MPI processes
    - data: The list of data to be distributed
    - verbose: A boolean flag to enable or disable verbose output

    Returns:
    --------
    - loc_data: A list of elements assigned to the current MPI rank
    - loc_num_elements: The number of elements assigned to the current MPI rank
    """

    if rank == 0:
        print(f"Total number of data decomposed: {len(data)}")
        numpairs = np.shape(data)[0]
        counts = np.arange(size, dtype=np.int32)
        displs = np.arange(size, dtype=np.int32)
        ave = int(numpairs / size)
        extra = numpairs % size
        offset = 0

        for i in range(0, size):
            loc_num_element = ave if i < size - extra else ave + 1
            counts[i] = loc_num_element

            if i == 0:
                loc_num_element0 = loc_num_element
                offset += loc_num_element
                displs[i] = 0
            else:
                comm.send(offset, dest=i)
                comm.send(loc_num_element, dest=i)
                offset += loc_num_element
                displs[i] = displs[i - 1] + counts[i - 1]

            if verbose:
                for j in range(offset - loc_num_element, offset):
                    print(f"Rank: {i} - {data[j]}")
                    sys.stdout.flush()

        offset = 0
        loc_num_element = loc_num_element0

    comm.Barrier()

    if rank != 0:  # workers
        offset = comm.recv(source=0)
        loc_num_element = comm.recv(source=0)

    comm.Barrier()
    loc_data = data[offset: offset + loc_num_element]
    return loc_data, loc_num_element
