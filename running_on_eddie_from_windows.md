- Copy the input files from local machine to Eddie using WinSCP, stfp, scp, or rsync.  Copy to a directory in ``/exports/eddie/scratch/<uun>/``.
	- set winSCP to "Transfer Settings: Text" to avoid invisible windows carriage return ^M characters
- Write/edit a submit script ``submit.sh`` such as this one (MPI) and copy it to eddie (to ``/home/<uun>/`` or project scratch directory).
	- ```bash
	  #!/bin/sh
	  # Grid Engine options (lines prefixed with #$)
	  #$ -N ogs5mpi 
	  #$ -wd  /exports/eddie/scratch/<uun>/3D     #<-- change directory of input files 
	  #$ -pe scatter 2                            #<--  change number to match domain decomp
	  #$ -l h_rt=00:03:00 
	  #$ -l h_vmem=8G
	  
	  . /etc/profile.d/modules.sh
	  
	  module load openmpi/1.10.1
	  module load geos/apps/ogs/5.8-UoE
	  
	  mpirun -np 2 ogs 3D                     #<--  change number processes and input (.pcs) file name
	  ```
	- choose a running time and a virtual memory requirement by looking at the output of ``qacct -j <jobid>`` of a similar simulation.
 	- more options: https://www.wiki.ed.ac.uk/display/ResearchServices/Job+Submission
- Open a Powershell command line and log into eddie: 
  ``ssh <uun>@eddie.ecdf.ed.ac.uk``
- On eddie, navigate to directory of submit script and submit the job:  ``qsub submit.sh `` .
  
- Check progress with ``qstat`` :  ``qw`` for waiting, ``r`` for running, ``E`` for error/exitting.  Delete jobs with ``qdel <jobID>`` .
  	- If the job is successful it runs for the required time, then output data files and log files appear in the project directory.
	- Troubleshooting:  Look in ``/exports/eddie/scratch/<uun>/<project>`` for files ``<jobname>.e<jobID>`` (error and warnings log), ``<jobname>.o<jobID>`` (output log) .  Also check output of command ``qacct -j <jobID>`` .
		- Stuck at ``Eqw`` with no output likely means a problem in the resource request ``#$`` lines of the script.  Check the path at ``-wd`` exists, check for invisble characters (run ``dos2unix`` on the submit script file), and check fail reason in ``qacct -j <jobID>``.
    
- Copy output data back using WinSCP, stfp, scp, or rsync .

  ## Builds
  All build include the latest changes from https://github.com/ufz/ogs5 as well as University of Edinburgh - specific additions.
  ### MPI
  The version at ``geos/apps/ogs/5.8-UoE`` is an MPI build, compiled with setting ``DOGS-CONFIG=MPI`` and built with ``openmpi/1.10.1``.  To run all but non-MPI jobs, The module ``openmpi/1.10.1`` must also be loaded.
  #### Input
  In addition to the usual files, a domain decomposition file ``.ddc`` matching the number of MPI cores must be present in the directory.  No further changes are needed to run in parallel.
  #### Solvers, .num file options
ILU precoditioner is not available in the MPI build.  The message "``Linear solver BiCGSTab with ILU not available. Use Jacobi:``" when ILU solver (key ``100``) is requested is hard-coded into the latest software version.  Use the MPI build with Jacobi preconditioner (key ``1``).  
Solvers ``2`` "BiCGStab", ``3`` "BiCGStab", ``5`` "CG", and ``7`` "CGS" are available with MPI builds.
 #### Outputs
 With option ``TECPLOT`` is ``.out`` file, a single TECPLOT file ``_domain_tet_0.tec`` gathering all spatial output is produced.  This is why the MPI build is prefered to PETSC builds.
  Options ``BINARY`` and ``VTK`` produce files by time step.
  ### PETSC-MPI
  The alternative PETSC build is at ??.  It must also be loaded with ``openmpi/1.10.1``.
  #### Input
   In addition to the usual files, three ``partition.msh`` or ``partition.bin`` must be present in the input files directory.  They can be generated with ``partmesh`` at https://github.com/wenqing/mesh_partition .  
   The ``-q`` flag must be added where the simulation includes mesh deformations and the ``-n`` flag is appropriate for this version of ogs5.  Flag ``-asci`` to generate human-readbale mesh formats is optional; the simulation works from either ``.bin`` or ``.msh`` files.

   The solver line(s) in ``.num`` files must be changed to a petsc solver, i.e.
```
petsc cg jacobi 1.e-16  2000 1.0
```
#### Output
Each parallel domain writes a separate output file ``_domain_tet_<n>.tec`` .
Alternatively binary files can be output and postprocessed to spatially joined ``.vtk`` files, one for each time step, with https://github.com/wenqing/post_process_ogs5.  
  ### Serial run
  Either MPI or Petsc builds can be used to run in serial, a serial script would look like 
```bash
	  #!/bin/sh
	  # Grid Engine options (lines prefixed with #$)
	  #$ -N ogs5mpi 
	  #$ -wd  /exports/eddie/scratch/<uun>/3D     #<-- change directory of input files 
	  #$ -l h_rt=00:10:00 
	  #$ -l h_vmem=8G
	  
	  . /etc/profile.d/modules.sh
	
	  module load geos/apps/ogs/5.8-UoE
	  
	  ogs 3D                     #<--  change input (.pcs) file name
```

### Desktop Build
Defaults: ``OGS_CONFIG`` was ``FEM`` and ``OGS_LSOLVER`` was ``RF``.  All preconditioners and solvers are available.  ``spBICGSTAB`` is used when preconditioner, solver options ``100``, ``2`` are requested in the num file.
