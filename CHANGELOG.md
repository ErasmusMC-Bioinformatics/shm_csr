version 1.3.3
-------------
+ Fix a bug where tandem lengths were incorrectly calculated.

version 1.3.2
-------------
+ Fix a bug where the file removal process caused errors.

version 1.3.1
-------------
+ Fix issues with container discovery

version 1.3.0
-------------
+ Add missing dependencies to the requirements section.

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
