# Minimal Aggregated SHared memory Messaging - Message Passing Layer (MASHM - MPL)

# Design decisions

The following library is designed to facilitate the use of MPI 3.0 shared memory features in user codes. Here, we specifically target nearest neighbor communication patterns which are common in high-performance scientific applications. 

The design of the API was chosen to balance simplifying the setup of shared memory communication schemes with the flexibility to allow for the overlap of computation and communication. The API calls handle many of the gory details of the setup and implementation of shared memory communication schemes, however the user needs to be aware of the order of API calls. In particular, since the intranodal communication can be separated from the internodal communication, once set up, a full communication exchange has the following form.

1. MashmInterNodeCommBegin(myMashm);
2. MashmIntraNodeCommBegin(myMashm);
3. MashmIntraNodeCommEnd(myMashm);
4. MashmInterNodeCommEnd(myMashm);

Although this requires four API calls, it provides the user the maximal opportunities to perform computation overlapped with communication.

# Usage

The assumptions of user codes are the following:

1. An existing point to point MPI communication exchange
2. Connectivity information explicitly available

With the above information, one can use the API provided in this library to specify necessary information and then use the shared memory communication methods available.

Fortran bindings for this library are provided, although the Fortran API differs from the C API slightly to handle multi-dimensional pointers.

Examples are provided in the test/ directory. Here, a two-dimensional domain is set up and a domain decomposition is performed. The resulting communication pattern is fed into the library which then sets upt he resulting communication reducing methods.

# Dependencies

By default, this package requires the use of MPI and specifically requires a version which supports the following MPI 3.0 features.

1. MPI_Comm_split_type - MPI_COMM_TYPE_SHARED
2. MPI_Win_allocate_shared

These are required by the library to implement the shared memory messaging. These features have been found to be supported in the following implementations and versions

1. OpenMPI >= 1.7.5
2. MVAPICH >= 2.0
3. MPICH >= 6.0.2
4. IMPI ? 5.0.0

To use the Fortran bindings of this library one must have a Fortran compiler which supports the 2003 standard. In particular, the Fortran implementation must support the following.

1. iso_c_binding, c_ptr
2. c_f_pointer

# How to build and install the package

The library uses CMake to build and install the library, as well as building and running several tests and examples. CMake supports (recommended) out of place builds.

cmake /path/to/source

# TODO:

1. Documentation - Doxygen
2. Compiler checks in CMake (MPI, F2003)
3. More robust examples
4. Refactor to eliminate dynamic message memory (require user to allocate number of MPI messages)
