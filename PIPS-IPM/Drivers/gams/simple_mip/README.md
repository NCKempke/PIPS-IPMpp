#Running simple-mip examples

* Make sure PIPS is installed and compiled in the build folder
  ```
  <path_to_pips>/build
  ```
  using the installation guide. There should be two executables ```gmspips``` and ```gmschk```.
* Make sure to have an updated GAMS version (GAMS version 33 or newer).
* Set the environment variable ```$GAMSSYSDIR``` to the path under which GAMS is installed:
    ```
    export GAMSSYSDIR=<path_to_gams>
    ```
  You can also adapt your ```.bashrc``` for this and run ```source ~/.bashrc```
* Running the provided script
  ```
  ./run_simple_mip_models.sh <number_mpi_procs>
  ```
  will run the given examples using ```<number_mpi_procs>``` processes (this cannot be more than you CPU allows for). The script can actually be run from any directory.
  
***NOTE:*** the script assumes ```mpirun``` as mpi-run command. Should your environment use some other command (like e.g. ```srun```) you will have to adapt the script manually 

***NOTE:*** you have to make sure that the corresponding script is indeed executable (possibly by running ```chmod +x run_simple_models.sh```)
