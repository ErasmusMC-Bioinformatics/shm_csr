version 1.8.1
-----------------
+ Fix a bug where input files with spaces could not be used.
+ Fix a bug were the data format was not properly set for filtered IMGT 
  archives.

version 1.8.0
-----------------
+ Add options "Everything is IGA" and "Everything is IGG". 
+ Make sure IGM naive and naive memory are part of the galaxy output.

version 1.7.0
-----------------
+ Use the name of the input file to generate the name of the output IMGT 
  archives.
+ Add same duplicate filters as immune repertoire pipeline.
+ Add a new "Everything is IGM" class filter for captured IGM sequences.
+ Fix bug where empty tables would cause crashes when generating plots.
+ Fix bug where R script errors where not written to stderr, causing galaxy to 
  mistake the jobs as being successful. 

version 1.6.0
-------------
+ Faster runtime due to faster gene identification, sequence overview creation 
  and IMGT TAR archive creation.
+ Two extra IMGT files are provided. One for IGM with less than 2% mutations
  (naive) and one for IGM with more than 2% mutations (naive memory).
+ All IMGT files per sequence class are always provided regardless of how the 
  ``Fast`` option is set. Previously this had to be set to ``no``.
+ Fix a bug in sequence overview where links to unmatched sequences where not
  working properly in the by_id.html file.

version 1.5.0
-------------
+ Add an option to download all output files in a zip file.

version 1.4.0
-------------
+ Fix a bug where synonymous mutations where incorrectly parsed.
+ Use a container from biocontainers.

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
