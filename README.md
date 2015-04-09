# Minimal Aggregated SHared memory Messaging (MASHM) - Message Passing Layer

Many high-performance distributed memory applications rely on
point-to-point messaging using the Message Passing Interface (MPI). Due
to latency of the network, and other costs, this communication can limit
the scalability of an application when run on large node-counts of
distributed memory supercomputers. Communication costs are further
increased on modern multi- and many-core architectures, when using more
than one MPI process per node, as each process sends and receives
messages independently, inducing multiple latencies and contention for
resources. 

We use shared memory constructs, available in
the MPI 3.0 standard, to implement an aggregated communication method to
minimize the number of inter-node messages while eliminating intra-node
messages altogether to reduce these costs. 
This is called the Minimal Aggregated SHared Memory (MASHM) communication method.

The MASHM library facilitates the use of the MASHM communication method in user codes. One only needs to set the number of messages, the destination and size of each message, and the library will set up the MASHM communication method and return pointers by which the original point-to-point data can be written and read. Through a straightforward sequence of calls, the user can then initiate minimal aggregated inter-node communication as well as exchange intra-node data. As MASHM is a layer on top of MPI, all MPI communication calls are handled internally to the library. The library provides consistent C and Fortran bindings. We hope that this library will encourage others to use the MASHM communication method.


# Using the MASHM library in applications

The usage of the MASHM in user codes assumes the following.

1. That a point-to-point MPI communication exchange currently exists
2. That the point-to-point messaging information (source, destination, and size of the messages) is explicitly available.

With the above information, one can use the API provided in this library to specify necessary information and then use the shared memory communication methods available.

Fortran bindings for this library are provided, although the Fortran API differs from the C API slightly to handle multi-dimensional pointers.

Examples are provided in the test/ directory. 

# Design decisions

The design of the API was chosen to balance simplifying the setup of shared memory communication schemes with the flexibility to allow for the overlap of computation and communication. The API calls handle many of the gory details of the setup and implementation of shared memory communication schemes, however the user needs to be aware of the order of API calls. In particular, since the intranodal communication can be separated from the internodal communication, once set up, a full communication exchange has the following form.

1. MashmInterNodeCommBegin(myMashm);
2. MashmIntraNodeCommBegin(myMashm);
3. MashmIntraNodeCommEnd(myMashm);
4. MashmInterNodeCommEnd(myMashm);

Although this requires four API calls, it provides the user the maximal opportunities to perform computation overlapped with communication.

# Dependencies

By default, this package requires the use of MPI and specifically requires a version which supports the following MPI 3.0 features.

1. MPI_Comm_split_type - MPI_COMM_TYPE_SHARED
2. MPI_Win_allocate_shared

These are required by the library to implement the shared memory messaging. These features have been found to be supported in the following implementations and versions

1. OpenMPI >= 1.7.5
2. MVAPICH >= 2.0
3. MPICH >= 6.0.2
4. IMPI >= 5.0.1

To use the Fortran bindings of this library one must have a Fortran compiler which supports the 2003 standard. In particular, the Fortran implementation must support the following.

1. iso_c_binding, c_ptr
2. c_f_pointer

# How to build and install the package

The library uses CMake to build and install the library, as well as building and running several tests and examples. CMake supports (recommended) out of place builds.

    jamroz@yslogin2:mashm-opt> cat config.sh
    #!/bin/bash

    rm -rf CMakeFiles CMakeCache.txt

    cmake \
      -DCMAKE_C_COMPILER="mpicc" \
      -DCMAKE_Fortran_COMPILER="mpif90" \
      -DCMAKE_C_FLAGS="-O3" \
      -DCMAKE_Fortran_FLAGS="-O3" \
      /path/to/source

    jamroz@yslogin2:mashm-opt> ./config.sh
    jamroz@yslogin2:mashm-opt> make -j 8 

This will produce executables under the following directory.

    /path/to/build/test

# How to build MASHM with GPTL timers

MASHM optionally can use the General Purpose Timing Library (GPTL), available at http://jmrosinski.github.io/GPTL/ , to provide timings of the communication routines. To enable these timers build and install GPTL (to say /path/to/gptl-install) and set the following configure time variable.

    -DGPTL_DIR=/path/to/gptl-install \

Ensure that the output of the configure step indicates that the GPTL library was found and then build the code as above (by typing make).

# TODO:

1. Extend the documentation to internal classes 
2. Compiler checks in CMake (MPI, F2003)
3. Remove unnecessary data from MashmPrivate
