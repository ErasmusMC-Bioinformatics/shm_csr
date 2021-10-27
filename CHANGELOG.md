version 1.2.0
-------------
+ Add a container in which the tool can execute.
+ Remove redundant files after execution.

version 1.1.0
-------------
+ Added changeo as a dependency. Porting to python3 was necessary to achieve 
  this. This will make sure the shm_csr package can be installed on all 
  galaxies.
+ Make sure the wrapper script runs with `set -e -o pipefail` and fails on 
  error.
+ Updated all python scripts to work on python3
