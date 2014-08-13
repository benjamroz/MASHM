# Minimal Aggregated SHared memory Messaging - Message Passing Layer (MASHM-MPL)

# Design decisions

The following library is designed to facilitate the use of MPI 3.0 shared memory features in user codes. Here, we specifically target nearest neighbor communication patterns which are common in high-performance scientific applications. 

# Usage

The assumptions of user codes are the following:

1. An existing point to point MPI communication exchange
2. Connectivity information explicitly available

With the above information, one can use the API provided in this library to specify necessary information and then use the shared memory communication methods available.

# Dependencies

By default, this package requires the use of MPI and specifically requires a version which supports the following MPI 3.0 features.

1. MPI_Comm_split_type - MPI_COMM_TYPE_SHARED
2. MPI_Win_allocate_shared

These are required by the library to implement the shared memory messaging. These features have been found to be supported in the following implementations and versions

1. OpenMPI >= 1.7.5
2. MVAPICH >= 2.0
3. MPICH >= 6.0.2
4. IMPI ? 5.0.0

# How to build and install the package

The library uses CMake to build and install the library, as well as building and running several tests and examples. CMake supports (recommended) out of place builds.

cmake /path/to/source

