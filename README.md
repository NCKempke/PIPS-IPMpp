# What is PIPS-IPM++?

PIPS-IPM++ is a (MPI + OpenMP) parallel interior-point method for doubly bordered block diagonal Linear Programms (QPs with linear constraints are currently not supported - we could do that if wished for). For more information on the current algorithm implemented in PIPS and on how to use it see [here](https://opus4.kobv.de/opus4-zib/files/7432/ip4energy.pdf). There exist two general purpose interfaces to PIPS-IPM++, one via [GAMS](https://www.gams.com/) and thus a GAMS license is required, and one esablished with the [METIS-project](https://www.fz-juelich.de/ias/jsc/EN/Research/Projects/_projects/metis.html) (on branch PyomoToPIPS). We encourage and support writing additional interfaces, too.

PIPS-IPM++ is a (highly modified) derivative of the [PIPS](https://github.com/Argonne-National-Laboratory/PIPS) solver originally developed at Argonne National Laboratory. You can find our most recent publications [here](https://www.sciencedirect.com/science/article/abs/pii/S0377221721005944) and [here](https://www.springerprofessional.de/en/first-experiments-with-structure-aware-presolving-for-a-parallel/18413092).

# INSTALLATION Instructions

Note that at this point we only support installation on Linux systems.

## General installation instructions
1. Install the packages **wget**, **cmake**, **mpich2**, and **boost**.
You can get them via the following command (xxxx stands for the name of the package):
In Linux(Ubuntu):\
   ```
   (sudo) apt install xxxx
   ```
   Note, that PIPS-IPM++ (and its third party libraries) heavily rely on a proper installation of **lapack** and **blas**. The correct libraries can be provided to PIPS' cmake in 3 ways (priority as follows):
 * The user can pass the variable MATH_LIBS to cmake passing all necessary linker information for custom defined lapack and blas to be used.
 * The user can define the environment variable MKLROOT to make PIPS use the Intel MKL libraries (suitable for Intel processing units). This environment variable can and should be set correctly by one of the scrips provided by Intel with each of its installations of MKL (for more information we recommend checking out the Intel MKL installation instructions).
 * The user can not pass/define anything and cmake will try to automatically find blas and lapack on the system (if installed)
In either way the cmake output shows the finally chosen lapack/blas routines.

2. Check out the current PIPS-IPM++ "master" **pipsipm** (or **PyomoToPIPS**, for using the Pyomo-Interface):
   ```
   git checkout linking-zib-pardiso-hierarchical
   ```
   It is vital that you check out the corresponding branch first since some of the installation scripts further down might have been modified in-between branches. If you have run a script before checking out linking-zib-pardiso-hierarchical and get errors from cmake later on, please delete the content of said folders and rerun the correct scripts. In addition it might help to delete the file CMakeCache.txt in you build folder (further down).
   

3. Go to the following folders and run the script wgetXXX.sh: **ThirdPartyLibs/METIS**.\
    For an example, use command ```sh wgetMETIS.sh``` in the folder ThirdPartyLibs/METIS.  


4. **OPTIONAL:** Download **MA27** and/or **MA57** from HSL and put the .tar.gz in the correct folder ThirdPartyLibs/MAXX and run the respective install script (see ThirdPartyLibs/MA27/README.txt and ThirdPartyLibs/MA57/README.txt for more details.)


5. **OPTIONAL:** Obtain a **PARDISO** (best >= 7.2) license and download PARDISO from [here](http://www.pardiso-project.org/).
   
    Copy the correct PARDISO library (the one compiled with GNU - this step might depend on your system and which version of gclib is used there - check out [this stack overflow answer](https://stackoverflow.com/questions/57155333/undefined-reference-to-log2fglibc-2-27-undefined-reference-to-logfglibc-2-2/61729571#61729571) for more info)) as **libpardiso.so** into the folder ThirdPartyLibs/PARDISO/src.
    ```
    cp pardiso_lib_name <pips_root>/ThirdPartyLibs/PARDISO/src
    ```
    Either copy the lincense (in a file named pardiso.lic) into your home directory or alternatively into the directory pips will be run from or, and this is the recommended way, set the environment variable PARDISO_LIC_PATH to the folder where the pardiso.lic file can be found (see also the [pardiso user guide](https://pardiso-project.org/manual/manual.pdf)).
    
    Even though this step is theoretically optional we **highly recommend** using PARDISO for proper performance of PIPS-IPM++. Currently MA27 and MA57 are not well integrated and tested (creating the folder will become deprecated in the future).

**NOTE:** At least one of the solvers MA27, MA57 or PARDISO must be available in order to proceed! So you must got through at least step 4 or 5!    

6. We should now be able to **compile PIPS-IPM++**. Assuming we are trying to install PIPS-IPM++ in the folder <pips_root>/build, where
   <pips_root> is the root installation folder, use the following commands to configure and install PIPS:
    ```
    cd <pips_root>
    mkdir build
    cd build
    cmake .. -DCMAKE_BUILD_TYPE=RELEASE
    make
    ```
    Alternatively:
    ```
    -DCMAKE_BUILD_TYPE=DEBUG
    ```
    Alternatively 
    ```
    make -j
    ``` 
    to use all your PCs available resources to compile PIPS-IPM++ or 
    ```
    make -jXX
    ```
    where XX is some max number of threads allowed eg. make -j10).
    
    Note: the ```cmake ..``` command should display a list of all solvers you installed. If your solver is not printed there it cannot be found and you might want to reconcile with steps 4 or 5.
   

7. Now the build directory should contain an executable called **pipsipmCallbackExample**. Run this to run PIPS for the first time and check whether everything is working.

## Running PIPS-IPM++

The actual PIPS-IPM++ executable is **gmspips** (or **pipsipmPyomo**). Consult the [best practice guide](https://gitlab.com/beam-me/bpg) Chapter 4.3 in order to learn on how to annotate models in GAMS and how to split them into multiple files so that PIPS-IPM++ can handle them.

1.  Let's assume, you already have a split model in SOMEFOLDER named model0.gdx, model1.gdx, ..., modelN.gdx where N is the number of blocks in your problem.
    Then, a typical call to PIPS-IPM++ running on **n processes** to solve the model looks like 
    ```
    mpirun -np n PIPSMAINPATH/build/gmspips N+1 SOMEFOLDER/model GAMSFOLDER scaleGeo stepLp presolve
    ```
    Note, that the actual MPI-Command **mpirun** might differ depending on your system and MPI installation (it could be **srun** etc.). You can always also run PIPS sequentially (leaving out the ```mpirun -np n``` part) but that would somehow defeat its purpose.
    Here GAMSFOLDER is the path your GAMS installations and the ```gams``` executable lie in. The last three arguments are optional but should be provided for best performance.\
    \
    When running **pipsipmPyomo** nothing really changes except for
    
 * The argument **GAMSFOLDER** is dropped.
 * Currently, the files have to be named Block0.txt, Block1.txt, ..., BlockN.txt
 * Instead of **N+1** the value **N** (so the actual amount of blocks must be specified)
   
    A typical call to **pipsipmPyomo** thus has the form
    ```
    mpirun -np n PIPSMAINPATH/build/gmspips N SOMEFOLDER/Block scaleGeo stepLp presolve
    ```
2. The **PIPSIPMpp.opt** file is PIPS-IPM++'s options file.\
There is not yet a detailed description of all adjustable parameters for the options file. Upon start PIPS-IPM++ will look for a file called PIPSIPMpp.opt in the current working directory (the directory where it is run from) and read in its options. Each line has the following structure:
    ```
    PARAMETERNAME VALUE TYPE
    ```
    where ```PARAMETERNAME``` is the name of the parameter to set, ```VALUE``` the value to set it to and ```TYPE``` one of ```bool```, ```int``` or ```double``` (depending on whether we are setting a bool, int or double parameter).
 
    The interested user can find a load of parameters and some short descriptions in the source files **PIPSIPMppOptions.C** and **Options.C** but expert knowledge of PIPS-IPM++ is currently required here and in question one should contact the developers (until a concise guide is available).
 

3. The three optional command line arguments turn on certain features.\
   The argument **presolve** activates PIPS-IPM++'s presolving, **scaleGeo** activates PIPS-IPM++'s geometric scaling and **stepLp** activates the use of primal and dual step lengths in the interior point method. We recommend to activate all of them.


4. The **number of threads** used by each MPI process in PIPS-IPM++ can be controlled by setting the environment variable ```OMP_NUM_THREADS```. PIPS-IPM++ will complain if this variable is not set. Best performance (from our point of view) can currently be achieved by setting
    ```
    export OMP_NUM_THREADS=2
    ```
    but a value of 1 is also fine and generally performance will depend on your instances.

### Profiling and timing for HPC 
PIPS-IPM++ + PARDISO offers best-in-class HPC performance. PIPS-IPM++ has built-in parallel performance profiling (mostly in the form of detailed timing and extended convergence reporting). To enable this feature, build PIPS with the -DWITH_TIMING option, for example, a typical build command would be
```{r, engine='bash', withtiming}
cmake -DCMAKE_TOOLCHAIN_FILE=../Toolchain.cmake -DWITH_TIMING=ON .. 
```
This is an expert option and requires in-depth knowledge of PIPS-IPM++'s algorithm. It is in refactoring currently.

# CONTRIBUTIONS
(in alphabetical order)

## PIPS-IPM++
Developed by:
  * Nils-Christian Kempke - Zuse-Institut Berlin/Technische Universität Berlin
  * Daniel Rehfeld - Zuse-Institut Berlin/Technische Universität Berlin
  * Charlie Vanaret - Technische Universität Berlin

Contributions from:
  * Svenja Uslu 

We thank for their work on interfaces/testing/bug reoprting and general support
  * Thomas Breuer (FZ Juelich)
  * Michael Bussieck (GAMS)
  * Henrik Büsing (FZ Juelich)
  * Frederick Fiand (GAMS)
  * Theresa Groß (FZ Juelich)
  * Maximilian Hoffmann (FZ Juelich)
  * Manuel Wetzel (DLR)

## PIPS-IPM
Originally developed by:
  * Cosmin G. Petra - Lawrence Livermore National Laboratory

Contributions from:
  * Miles Lubin - while at Argonne National Laboratory

# LICENSE

PIPS-IPM++ is derivative work of PIPS-IPM. While the original code retains it old licensing all PIPS-IPM++ additions are distirbuted under the [LPGLv3](https://www.gnu.org/licenses/lgpl-3.0.en.html).

See LICENSE file.

PIPS-IPM is derivative work of OOQP (http://pages.cs.wisc.edu/~swright/ooqp/) by E. Michael Gertz and Stephen. Wright

Currently PIPS is not public / open souce but this is about to change.
# ACKNOWLEDGMENTS

PIPS-IPM++ has been developed under the financial support of:
- BEAM-ME (BMWi project, ID: 03ET4023A-F)
- UNSEEN (BMWi project, ID: 03EI1004-C)
- Gauss Centre for Supercomputing e.V. (www.gauss-centre.eu) and the John von Neumann Institute for Computing (NIC) on the GCS Supercomputer JUWELS at Jülich Supercomputing Centre (JSC)
- Research Campus MODAL (BMBF grant numbers 05M14ZAM, 05M20ZBM)

PIPS has been developed under the financial support of: 
- Department of Energy, Office of Advanced Scientific Computing Research
- Department of Energy, Early Career Program 
- Department of Energy, Office of Electricity Delivery and Energy Reliability
