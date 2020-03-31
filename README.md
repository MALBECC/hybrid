Hybrid
============

Hybrid is QM/MM code designed primarily to obtain potential energy profiles in complex systems.

REQUIREMENTS
------------

* LAPACK or INTEL MKL.
* GNU or INTEL Fortran Compiler.
* GNU Make.
* QM external software, LIO or Orca.

COMPILATION
------------

The program can be compiled using the make command. The following options can be used to modify compilation. For example, the following compiles intel compilers:

```
make intel=1
```

* _intel_ : Enables the use of intel compilers.
    0 = no intel compilers; use gnu fortran compiler and libraries.
    1 = use the intel compilers instead of the GNU ones.
    2 = use the intel compilers and the MKL library for lapack.

* _analytics_ : Enables profiling or debugging information and flags.
    0 = No extra information is provided.
    1 = Profiling information is enabled.
    2 = Basic debugging information is enabled.
    3 = Detailed debugging information is enabled.

* _qm_lio_ : Enable Lio as external program to compute QM & QM/MM energy and forces
    0 = disable
    1 = enable

* _qm_orca_ : Enable Orca as external program to compute QM & QM/MM energy and forces
    0 = disable
    1 = enable



INSTALLATION
------------

Compilation will produce an executable file in /bin.


INSTALLATION WITH LIO
-------------------------

Define LIOHOME=lio_path/lio and compile with _qm_lio=1_ option.


INSTALLATION WITH ORCA
-------------------------

Compile with _qm_orca=1_ option, the Orca path should be defined in .fdf file.


TESTS
-----

The test suite can be ran from the tests directory, each subfolder contains a "run.scr" script which performs the test.


CONTRIBUTING
------------
Any contribution will be appreciated, at this moment this is a work in progress project.

PUBLICATIONS
------------

1. 
