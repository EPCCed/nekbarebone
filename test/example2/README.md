This example will run one nekbone test with:

50  elements per process
10  polynomial order

It will run without the multigrid preconditioner.

The following commands will build nekbarebone for
the example2 test case. Please note, you must first set
the SOURCE_ROOT environment variable in the makenek script.

makenek clean
makenek ex2

On ARCHER, nekbarebone can be compiled with the Cray, Intel
and gnu programming environments.

The "qsub submit.pbs" command will submit the nekbarebone
job to the standard queue. Please note, you must first set
the PBS account variable in the submit.pbs submission script.
