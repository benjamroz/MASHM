===============================
The MASHM Library User's Manual
===============================

Minimally Aggregated SHared memory Messaging (MASHM) Library
============================================================

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

Requirements
------------

The MASHM library requires the following software.

- CMake (>= 2.8)
- Modern C compiler
- Optional modern Fortran compiler
- MPI implementations supporting the MPI 3.0 standard

How can I get it?
=================

The best way to obtain the MASHM code is to check it out from the
publicly available git repo at github.

::

    $ git clone https://github.com/benjamroz/MASHM.git mashm

This is the most recent stable version of the source code. The trunk is
also available for download, if you choose to have the most up-to-date
version.

How do I set it up?
===================

Installation
------------

In this section, we describe how to install the PyReshaper package on a
unix-like system. The procedure is similar for a Mac, but we have not
tested the package on Windows.

As described in the previous section, first check out the source code
from the subversion repository. On unix-like systems, the command is
shown below.

::

    $ git clone https://github.com/benjamroz/MASHM.git mashm

The build system supports out of place builds. Note that this is the best option for doing development. Create a "build directory" in which to compile the code.

::

    $ mkdir mashm-build

and enter the directory.

::

    $ cd mashm-build

The best way to configure the mashm build is to create a config.sh script containing options to pass to CMake.

::

    $ cat config.sh
    #!/bin/bash

    rm -rf CMakeFiles CMakeCache.txt

    cmake \
      -DCMAKE_C_COMPILER=mpiicc \
      -DCMAKE_Fortran_COMPILER=mpiifort \
      -DCMAKE_C_FLAGS="-O3" \
      -DCMAKE_Fortran_FLAGS="-O3" \
      -DCMAKE_INSTALL_PREFIX=/path/to/install/directory/mashm \
      /path/to/source/mashm

Here, the parallel C and Fortran MPI compiler wrappers are specified, the compiler flags for each of these compilers are specified, the installation directory is specified by the ``-DCMAKE_INSTALL_PREFIX`` argument, and the last argument ``/path/to/source/mashm`` specifies the location of the source directory. Note that the MPI compilers must support the MPI 3.0 standard.

Next, to build the code simply type make.

::

    $ make -j 8 

This will compile the MASHM library as well as the example executables. Finally, to install the library and header files to the location specified above simply type the following.

::

    $ make install

How to build MASHM with GPTL timers
-----------------------------------

MASHM optionally can use the General Purpose Timing Library (GPTL), available at [http://jmrosinski.github.io/GPTL/](http://jmrosinski.github.io/GPTL/), to provide timings of the communication routines. To enable these timers build and install GPTL (to say /path/to/gptl-install) and set the following configure time variable.

    -DGPTL_DIR=/path/to/gptl-install \

Ensure that the output of the configure step indicates that the GPTL library was found and then build the code as above (by typing make).


Generating the User Documentation
---------------------------------

The ``README.rst`` file and this User Manual should be consulted for help
on installing and using the software. Both documents are included with
the source. The ``README.rst`` file is included with the top-level
MASHMdirectory, and the User Manual is contained in the
``UserManual.rst`` file. Both files are reStructuredText formatted
files, meaning they are simple text files that can be read with any text
viewer.

An HTML version of the User Manual will be created (provided you have 
sphinx installed) by running

::

    $ make manual

and will be located in the build directory under
``doc/manual/html/UserManual.html`` or after installation in the install directory under 
the same relative path.

If you are a developer, you may find the Doxygen-generated API
documentation helpful in understanding the design and functionality of
the MASHM code. To generate this documentation, you must have
Doxygen available and installed. If you do, the API documentation can be
easily generated with the following.

::

    $ make doc

The Doxygen API documentation will be placed in ``doc/html/index.html``.

Running the example exectuables
-------------------------------

Several examples are provided to demonstrate the usage of MASHM. These source files for these examples are located in the test subdirectory of the source code : nodalComm.c, nodalCommFtn.F90, and poisson3d.F90. 

The nodalComm.c driver performs a domain decomposition of a two-dimensional rectangular domain and assigns MPI processes to individual elements. The connectivity information between elements, including the MPI rank of neighbors and the number of points shared, is given for each process.  Then a standard non-blocking point-to-point MPI communication is set up and run to exchange data between processes. Next, MASHM is used to exchange the same information. Finally, reductions are performed to ensure that the two methods achieve the same result.

The driver nodalCommFtn.F90 is the same as the C program just described except that it uses Fortran wrappers to the C library. This is useful to see the translation of the C call statements into Fortran.

The poisson3d.F90 driver performs a relaxation of a three-dimensional anisotropic Laplace's equation using standard non-blocking point-to-point MPI communication as well as the MASHM communication method. Here a three-dimensional domain is decomposed across MPI processes, the MPI process connectivity information is given and used to set up the standard communication scheme as well as the MASHM library.


Using the MASHM Library in Applications
---------------------------------------

The usage of the MASHM in user codes assumes the following.

1. That a point-to-point MPI communication exchange currently exists
2. That the point-to-point messaging information (source, destination, and size of the messages) is explicitly available.

With the above information, one can use the API provided in this library to specify necessary information and then use the shared memory communication methods available.

Fortran bindings for this library are provided, although the Fortran API differs from the C API slightly to handle multi-dimensional pointers.

The design of the API was chosen to balance simplifying the setup of shared memory communication schemes with the flexibility to allow for the overlap of computation and communication. The API calls handle many of the gory details of the setup and implementation of shared memory communication schemes, however the user needs to be aware of the order of API calls. In particular, since the intranodal communication can be separated from the internodal communication, once set up, a full communication exchange has the following form.

1. MashmInterNodeCommBegin(myMashm);
2. MashmIntraNodeCommBegin(myMashm);
3. MashmIntraNodeCommEnd(myMashm);
4. MashmInterNodeCommEnd(myMashm);

Although this requires four API calls, it provides the user the maximal opportunities to perform computation overlapped with communication.

The usage of MASHM in an application has the following form. First one declares and initializes the MASHM object.

.. code-block:: c

    /* Declare the MASHM object */
    Mashm myMashm;

    /* Initialize the MASHM object */
    MashmInit(&myMashm, MPI_COMM_WORLD);

Note that the initialization requires an MPI communicator. One can print the MPI process decomposition onto nodes by the following call.

.. code-block:: c

    /* Print nodal MPI process decomposition info */
    MashmPrintInfo(myMashm);

The number of messages each MPI process sends, in a standard point-to-point communication method, is set with

.. code-block:: c

    /* Set the number of messages for each processes */
    MashmSetNumComms(myMashm, numNeighbors);

and the destination and size of each message is set with the following.

.. code-block:: c

    /* Add the destination and size of each message */
    for (i = 0; i < numNeighbors; i++) {
      MashmSetComm(myMashm, i, msgDest[i], msgSizes[i]);
    }

After all of the messages sizes and destinations have been specified the MASHM precalculation stage can occur. This allocates shared memory for nodal messages as well as intra-node messages. This is accomplished with a call to MashmCommFinish.

.. code-block:: c

    /* Perform MASHM precalculation */
    MashmCommFinish(myMashm);

One can print out all of the messages with a call to MashmPrintCommCollection.

.. code-block:: c

    /* Print the communication collection */
    MashmPrintCommCollection(myMashm);

Next, we need to retrieve the pointers to the MASHM shared memory regions. Thus we need an array of pointers for each the send and receive buffers.

.. code-block:: c

    /* Retrieve pointers for buffers */
    mashmSendBufferPtrs = (double**) malloc(sizeof(double*)*numNeighbors);
    mashmRecvBufferPtrs = (double**) malloc(sizeof(double*)*numNeighbors);

We can then retrieve the pointers for the shared memory send and receive regions with calls to MashmGetBufferPointer specifying the send or receive buffer as follows.

.. code-block:: c

    /* Retrieve the pointer to access the MASHM's storage for each message */
    for (i = 0; i < numNeighbors; i++) {
      mashmSendBufferPtrs[i] = MashmGetBufferPointer(myMashm, i, MASHM_SEND);
      mashmRecvBufferPtrs[i] = MashmGetBufferPointer(myMashm, i, MASHM_RECEIVE);
    }

Now that we have the pointers set up we can begin filling the MASHM memory regions with the data that we want messages. We can first set only the inter-node data.
   
.. code-block:: c

    /* Fill internode buffers */
    for (i = 0; i < numNeighbors; i++) {
      if (! MashmIsMsgIntraNodal(myMashm, i)) {
        for (j = 0; j < msgSizes[i]; j++) {
          mashmSendBufferPtrs[i][j] = rank*msgSizes[i]+j;
        }
      }
    }

and then begin the sending of the nodal messages with the following non-blocking call.

.. code-block:: c

    /* Begin sending internode messages */
    MashmInterNodeCommBegin(myMashm);

Next, we can set the data for the intra-node messages

.. code-block:: c

    /* Messages sent and receives posted 
     * Can asynchronously do work on intra-node data 
     */
    for (i = 0; i < numNeighbors; i++) {
      if (MashmIsMsgIntraNodal(myMashm, i)) {
        for (j = 0; j < msgSizes[i]; j++) {
          mashmSendBufferPtrs[i][j] = rank*msgSizes[i]+j;
        }
      }
    }

and begin sending the intra-node data.

.. code-block:: c

    /* Send intranode messages */
    MashmIntraNodeCommBegin(myMashm);

At this point, all of the communication has been posted in a non-blocking manner. One can overlap the communication with some computation etc. Here we can allocate an array.

.. code-block:: c

    /* Asynchronously do some computation */
    mashmData = (double*) malloc(sizeof(double)*sumMsgSizes);

We can then wait for the intra-node data to be successfully completed

.. code-block:: c

    MashmIntraNodeCommEnd(myMashm);

and then retrieve this data from MASHM into our application.

.. code-block:: c

    /* Asynchronously do work on nodal data */
    for (i = 0; i < numNeighbors; i++) {
      if (MashmIsMsgIntraNodal(myMashm, i)) {
        /* Unpack individual buffer */
        offset = msgOffsets[i];
        for (j = 0; j < msgSizes[i]; j++) {
          mashmData[offset+j] = mashmRecvBufferPtrs[i][j];
        }
      }
    }

Finally, we can wait on the nodal messages

.. code-block:: c

    /* Now wait on nodal messages */
    MashmInterNodeCommEnd(myMashm);

and retrieve this data as well.

.. code-block:: c

    for (i = 0; i < numNeighbors; i++) {
      if (! MashmIsMsgIntraNodal(myMashm, i)) {
        /* Unpack individual buffer */
        offset = msgOffsets[i];
        for (j = 0; j < msgSizes[i]; j++) {
          mashmData[offset+j] = mashmRecvBufferPtrs[i][j];
        }
      }
    }

To clean up, we can retire the pointers to the MASHM data

.. code-block:: c

    /* Retire the Mashm buffer pointers */
    for (i = 0; i < numNeighbors; i++) {
      MashmRetireBufferPointer(myMashm, &(mashmSendBufferPtrs[i]));
      MashmRetireBufferPointer(myMashm, &(mashmRecvBufferPtrs[i]));
    }

and destroy the MASHM object.

.. code-block:: c

    /* Destroy the Mashm object */
    if (rank == 0) printf("Calling Mashm Destroy.\n");
    MashmDestroy(&myMashm);



