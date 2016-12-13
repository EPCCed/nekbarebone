nekbarebone is a cutdown version of nekbone, featuring just the original example2 nekbone test case.
nekbarebone is based on version 3.1 of nekbone taken from https://cesar.mcs.anl.gov/content/software/thermal_hydraulics.

All comms specific functionality (i.e., MPI calls) has been moved to the comm files in the jl sub-directory.
Furthermore, the compilation scripts have been altered such that the jl source is built has a static library,
which is then linked to the nekbarebone executable.

By default the jl source is compiled to exags_mpi.a, however, the nekbarebone executable can instead be made
to link to external exags libraries (such as one based on upc for example) by making suitable edits to the
./src/makefile.template file: first, set the IFEXAGSMPI variable to false and then set EXAGSLIBNAME and EXAGSLIBPATH
to the name and path of the external exags library. 
