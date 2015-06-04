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

Generating the User Documentation
---------------------------------

The ``README.rst`` file and this User Manual should be consulted for help
on installing and using the software. Both documents are included with
the source. The ``README.rst`` file is included with the top-level
MASHMdirectory, and the User Manual is contained in the
``UserManual.rst`` file. Both files are reStructuredText formatted
files, meaning they are simple text files that can be read with any text
viewer.

An HTML version of the User Manual will automatically be created by
Doxygen, as described in the previous section. A link will be created
from the "Related Pages" tab on the Doxygen-generated page.

::

    $ doxygen Doxyfile

If you are a developer, you may find the Doxygen-generated API
documentation helpful in understanding the design and functionality of
the PyReshaper code. To generate this documentation, you must have
Doxygen available and installed. If you do, the API documentation can be
easily generated with the following command from the top-level
PyReshaper directory.

::

    $ doxygen Doxyfile

The API documentation will be placed in the ``docs/html/`` directory.

Running the example exectuables
-------------------------------

Several examples are provided to demonstrate the usage of MASHM. These source files for these examples are located in the test subdirectory of the source code : nodalComm.c, nodalCommFtn.F90, and poisson3d.F90. 

The nodalComm.c driver performs a domain decomposition of a two-dimensional rectangular domain and assigns MPI processes to individual elements. The connectivity information between elements, including the MPI rank of neighbors and the number of points shared, is given for each process.  Then a standard non-blocking point-to-point MPI communication is set up and run to exchange data between processes. Next, MASHM is used to exchange the same information. Finally, reductions are performed to ensure that the two methods achieve the same result.

The driver nodalCommFtn.F90 is the same as the C program just described except that it uses Fortran wrappers to the C library. This is useful to see the translation of the C call statements into Fortran.

The poisson3d.F90 driver performs a relaxation of a three-dimensional anisotropic Laplace's equation using standard non-blocking point-to-point MPI communication as well as the MASHM communication method. Here a three-dimensional domain is decomposed across MPI processes, the MPI process connectivity information is given and used to set up the standard communication scheme as well as the MASHM library.

Example of Usage
----------------

::

    $ doxygen Doxyfile

If you are a developer, you may find the Doxygen-generated API
documentation helpful in understanding the design and functionality of
the PyReshaper code. To generate this documentation, you must have
Doxygen available and installed. If you do, the API documentation can be
easily generated with the following command from the top-level
PyReshaper directory.

::

    $ doxygen Doxyfile

The API documentation will be placed in the ``docs/html/`` directory.

Generating the User Documentation
---------------------------------

The ``README.md`` file and this User Manual should be consulted for help
on installing and using the software. Both documents are included with
the source. The ``README.md`` file is included with the top-level
PyReshaper directory, and the User Manual is contained in the
``docs/user/UserManual.md`` file. Both files are Markdown formatted
files, meaning they are simple text files that can be read with any text
viewer.

An HTML version of the User Manual will automatically be created by
Doxygen, as described in the previous section. A link will be created
from the "Related Pages" tab on the Doxygen-generated page.

Before Using the PyReshaper
---------------------------

After the Pyreshaper package has been installed using the procedure
above, you must add the installation site-packages directory to your
``PYTHONPATH``. If you installed with the ``--user`` option, this means
adding the ``$HOME/.local/lib/python2.X/site-packages`` directory to
your ``PYTHONPATH``. If you specified a different ``--prefix`` option,
then you must point to that prefix directory. For bash users, this is
done with the following command.

::

    $ export PYTHONPATH=$PYTHONPATH:$PREFIX/lib/python2.X/site-packages

where the ``$PREFIX`` is the root installation directory used when
installing the PyReshaper package (``$HOME/.local/`` if using the
``--user`` option), and the value of ``X`` will correspond to the
version of Python used to install the PyReshaper package.

If you want to use the command-line interface to the PyReshaper, you
must also add the PyReshaper executables directory to your ``PATH``.
Like for the ``PYTHONPATH``, this can be done with the following
command.

::

    $ export PATH=$PATH:$PREFIX/bin

How do I use it?
================

Some General Concepts
---------------------

Before we describe the various ways you can use the PyReshaper, we must
describe more about what, precisely, the PyReshaper is designed to do.

As we've already mentioned, the PyReshaper is designed to convert a set
of NetCDF files from time-slice (i.e., multiple time-dependent variables
with one time-value per file) format to time-series (one time-dependent
variable with multiple time-values per file) format. This statement
contains a number of assumptions that pertain to the time-slice (input)
data, which we list below.

1. Each time-slice NetCDF file has multiple time-dependent variables
   inside it, but can have many time-independent variables inside it, as
   well.
2. Each time-slice NetCDF file contains data for times that do not
   overlap with each other. (That is, each time-slice NetCDF file can
   contain data spanning a number of simulation time steps. However, the
   span of time contained in one time slice cannot overlap the span of
   time in another time-slice.)
3. Every time-slice NetCDF file contains the same time-dependent
   variables, just at differing times.

Similarly, there are a number of assumptions made about the time-series
data produced by the PyReshaper conversion process.

1. By default, every time-dependent variable will be written to its own
   time-series NetCDF file.
2. Any time-dependent variables that should be included in every
   time-series file (e.g., such as ``time`` itself), instead of getting
   their own time-series file, must be specified by name.
3. Every time-independent variable that appears in the time-slice files
   will be written to every time-series file.
4. Every time-series file written by the PyReshaper will span the total
   range of time spanned by all time-slice files specified.
5. Every time-series file will be named with the same prefix and suffix,
   according to the rule:

   time\_series\_filename = prefix + variable\_name + suffix

where the variable\_name is the name of the time-dependent variable
associated with that time-series file.

It is important to understand the implications of the last assumption on
the list above. Namely, it is important to note what this assumption
means in terms of NetCDF file-naming conventions. It is common for the
file-name to contain information that pertains to the time-sampling
frequency of the data in the file, or the range of time spanned by the
time-series file, or any number of other things. To conform to such
naming conventions, it may be required that the total set of time-slice
files that the user which to convert to time-series be given to the
PyReshaper in multiple subsets, or chunks. Throughout this manual, we
will refer to such "chunks" as streams. As such, every single PyReshaper
operation is designed to act on a single stream.

Using the PyReshaper from within Python
---------------------------------------

Obviously, one of the advantages of writing the PyReshaper in Python is
that it is easy to import features (modules) of the PyReshaper into your
own Python code, as you might link your own software tools to an
external third-party library. The library API for the PyReshaper is
designed to be simple and light-weight, making it easy to use in your
own Python tools or scripts.

Single-Stream Usage
~~~~~~~~~~~~~~~~~~~

Below, we show an example of how to use the PyReshaper from within
Python to convert a single stream from time-slice format to time-series
format.

.. code:: py

    from pyreshaper import specification, reshaper

    # Create a Specifier object (that defined a single stream to be converted
    specifier = specification.create_specifier()

    # Specify the input needed to perform the PyReshaper conversion
    specifier.input_file_list = [ "/path/to/infile1.nc", "/path/to/infile2.nc", ...]
    specifier.netcdf_format = "netcdf4c"
    specifier.output_file_prefix = "/path/to/outfile_prefix."
    specifier.output_file_suffix = ".000101-001012.nc"
    specifier.time_variant_metadata = ["time", "time_bounds", ...]

    # Create the Reshaper object
    rshpr = reshaper.create_reshaper(specifier, serial=False, verbosity=1)

    # Run the conversion (slice-to-series) process
    rshpr.convert()

    # Print timing diagnostics
    rshpr.print_diagnostics()

In the above example, it is important to understand the input given to
the PyReshaper. Namely, all of the input for this single stream is
contained by a single instantiation of a Specifier object (the code for
which is defined in the specification module). We will describe each
attribute of the Specifier object below.

Specifier Object Attributes
^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  ``input_file_list``: This specifies a list of input (time-slice) file
   paths that all conform to the input file assumptions (described
   above). The list of input files need not be time-ordered, as the
   PyReshaper will order them appropriately. (This means that this list
   can easily be generated by using filename globs.)

In the example above, each file path is full and absolute, for safety's
sake.

-  ``netcdf_format``: This is a string specifying what NetCDF format
   will be used to write the output (time-series) files.

In the above example, NetCDF4 with level-1 compression is requested.

Acceptable Options are:

-  ``"netcdf"``: NetCDF3
-  ``"netcdf4"``: NetCDF4 uncompressed
-  ``"netcdf4c"``: NetCDF4 compressed (level 1)

-  ``output_file_prefix``: This is a string specifying the common output
   (time-series) filename prefix. It is assumed that each time-series
   file will be named according to the rule:

   filename = prefix + variable\_name + suffix

It is important to understand, as in the example above, that the prefix
can include the full, absolute path information for the output
(time-series) files.

-  ``output_file_suffix``: This is a string specifying the common output
   (time-series) filename suffix. It is assumed that each time-series
   file will be named according to the above rule.

-  ``time_variant_metadata``: This specifies a list of variable names
   corresponding to variables that should be written to every output
   (time-series) NetCDF file.

Even though the PyReshaper is designed to work on a single stream at a
time, multiple streams can be defined as input to the PyReshaper. When
running the PyReshaper with multiple stream, multiple Specifier objects
must be created, one for each stream.

Multiple Stream Usage
~~~~~~~~~~~~~~~~~~~~~

In the example below, we show one way to define a multiple stream
PyReshaper run.

.. code:: py

    from pyreshaper import specification, reshaper

    # Assuming all data defining each stream is contained 
    # in a list called "streams"
    specifiers = {}
    for stream in streams:
        specifier = specification.create_specifier()

        # Define the Pyreshaper input for this stream
        specifier.input_file_list = stream.input_file_list
        specifier.netcdf_format = stream.netcdf_format
        specifier.output_file_prefix = stream.output_file_prefix
        specifier.output_file_suffix = stream.output_file_suffix
        specifier.time_variant_metadata = stream.time_variant_metadata

        # Append this Specifier to the dictionary of specifiers
        specifiers[stream.name] = specifier

    # Create the Reshaper object
    rshpr = reshaper.create_reshaper(specifiers, serial=False, verbosity=1)

    # Run the conversion (slice-to-series) process
    rshpr.convert()

    # Print timing diagnostics
    rshpr.print_diagnostics()

In the above example, we assume the properly formatted data (like the
data shown in the single-stream example above) is contained in the list
called *streams*. In addition to the data needed by each Specifier
(i.e., the data defining each stream), this example assumes that a name
has been given to each stream, contained in the attribute "stream.name".
Each Specifier is then contained in a dictionary with keys corresponding
to the stream name and values corresponding to the stream Specifier.
This name will be used when printing diagnostic information during the
``convert()`` and ``print_diagnostics()`` operations of the PyReshaper.

Alternatively, the specifiers object (in the above example) can be a
Python list, instead of a Python dictionary. If this is the case, the
list of Specifier objects will be converted to a dictionary, with the
keys of the dictionary corresponding to the list index (i.e., an
integer).

Arguments to the ``create_reshaper()`` Function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In both examples above, the Reshaper object (rshpr) is created by
passing the single Specifier object, list of Specifier objects, or
dictionary of named Specifier objects, to the function
``create_reshaper()``. This function returns a Reshaper object that has
the functions ``convert()`` and ``print_diagnostics()`` that perform the
time-slice to time-series conversion step and print useful timing
diagnostics, respectively.

Additionally, the ``create_reshaper()`` function takes the parameter
``serial``, which can be ``True`` or ``False``, indicating whether the
Reshaper ``convert()`` step should be done in serial (``True``) or
parallel (``False``). By default, parallel operation is assumed if this
parameter is not specified.

The ``create_reshaper()`` function also takes the parameter
``verbosity``, which specified what level of output (to ``stdout``) will
be produced during the ``convert()`` step. Currently, there are only
three (3) verbosity levels:

1. ``verbosity = 0``: This means that no output will be produced unless
   specifically requested (i.e., by calling the ``print_diagnostics()``
   function).
2. ``verbosity = 1``: This means that only output that would be produced
   by the head rank of a parallel process will be generated.
3. ``verbosity = 2``: This means that all output from all processors
   will be generated, but any output that is the same on all processors
   will only be generated once.

By setting the ``verbosity`` parameter in the ``create_reshaper()``
function to a value of 2 or above will result in the greatest amount of
output.

Arguments to the ``convert()`` Function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

While not shown in the above examples, there is an argument to the
``convert()`` function of the Reshaper object called ``output_limit``.
This argument sets an integer limit on the number of time-series files
generated during the ``convert()`` operation (per processor). This can
be useful for debugging purposes, as it can greatly reduce the length of
time consumed in the ``convert()`` function. (A value of ``0`` indicates
no limit, or all output files will be generated.)

Using the PyReshaper from the Unix Command-Line
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

While the most flexible way of using the PyReshaper is from within
Python, as described above, it is also possible to run the PyReshaper
from the command-line. In this section, we describe how to use the
Python script ``slice2series``, which provides a command-line interface
(CLI) to the PyReshaper. (This script will be installed in the
``$PREFIX/bin`` directory, where ``PREFIX`` is the installation root
directory.)

Below is an example of how to use the PyReshaper CLI, ``slice2series``,
for a serial run.

::

    $ slice2series --serial \
      --netcdf_format="netcdf4c" \
      --output_prefix="/path/to/outfile_prefix." \
      --output_suffix="000101-001012.nc" \
      -m "time" -m "time_bounds" \
      /path/to/infiles/*.nc

In this example, you will note that we have specified each
time-dependent metadata variable name with its own ``-m`` option. (In
this case, there are only 2, ``time`` and ``time_bounds``.) We have also
specified the list of input (time-slice) files using a wildcard, which
the Unix shell fills in with a list of all filenames that match this
pattern. (In this case, it is all files with the ``.nc`` file extension
in the directory ``/path/to/infiles``.) These command-line options and
arguments specify all of the same input passed to the Specifier objects
in the examples of the previous section.

For parallel operation, one must launch the ``slice2series`` script from
the appropriate MPI launcher. On the Yellowstone system
(``yellowstone.ucar.edu``), this is done with the following command.

::

    $ mpirun.lsf slice2series \
      --netcdf_format="netcdf4c" \
      --output_prefix="/path/to/outfile_prefix." \
      --output_suffix="000101-001012.nc" \
      -m "time" -m "time_bounds" \
      /path/to/infiles/*.nc

In the above example, this will launch the ``slice2series`` script into
the MPI environment already created by either a request for an
interactive session or from an LSF submission script.

Additional Arguments to the ``slice2series`` Script
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

While the basic options shown in the previous two (2) examples above are
sufficient for most purposes, two additional options are available. The
``--verbosity`` option can be used to set the verbosity level, just like
the ``verbosity`` argument to the ``create_reshaper()`` function
described in the previous sections. Additionally, the ``--limit``
command-line option can be used to set the ``output_limit`` argument of
the Reshaper ``convert()`` function, also described in the previous
sections.

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

MASHM optionally can use the General Purpose Timing Library (GPTL), available at [http://jmrosinski.github.io/GPTL/](http://jmrosinski.github.io/GPTL/), to provide timings of the communication routines. To enable these timers build and install GPTL (to say /path/to/gptl-install) and set the following configure time variable.

    -DGPTL_DIR=/path/to/gptl-install \

Ensure that the output of the configure step indicates that the GPTL library was found and then build the code as above (by typing make).

# TODO:

1. Extend the documentation to internal classes 
2. Compiler checks in CMake (MPI, F2003)
3. Remove unnecessary data from MashmPrivate
