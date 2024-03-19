- Copy the input files from local machine to Eddie using WinSCP, stfp, scp, or rsync.  Copy to a directory in ``/exports/eddie/scratch/<uun>``
- Open a Powershell command line and log into eddie: 
  ``ssh <uun>@eddie.ecdf.ed.ac.uk``
- Write/edit a submit script ``submit.sh`` such as this one (MPI) and copy it to eddie (to ``\home\<uun>\`` or project scratch directory).
	- ```bash
	  #!/bin/sh
	  # Grid Engine options (lines prefixed with #$)
	  #$ -N ogs5mpi 
	  #$ -wd  /exports/eddie/scratch/<uun>/3D     #<-- change directory of input files 
	  #$ -pe scatter 2                            #<--  change number to match domain decomp
	  #$ -l h_rt=00:03:00 
	  #$ -l h_vmem=8G
	  
	  . /etc/profile.d/modules.sh
	  
	  module load openmpi/3.1.6
	  module load geos/apps/ogs/5.8-UoE
	  
	  mpirun -np 2 ogs 3D                     #<--  change number processes and input (.pcs) file name
	  ```
	- choose a running time and a virtual memory requirement by looking at the output of ``qacct -j <jobid>`` of a similar simulation.
 	- more options: https://www.wiki.ed.ac.uk/display/ResearchServices/Job+Submission
- On eddie, navigate to directory of submit script and submit the job:  ``qsub submit.sh `` .
  
- Check progress with ``qstat`` :  ``qw`` for waiting, ``r`` for running, ``E`` for error/exitting.  Delete jobs with ``qdel <jobID>`` .
  	- If the job is successful it runs for the required time, then output data files and log files appear in the project directory.
	- Troubleshooting:  Look in ``/exports/eddie/scratch/<uun>/<project>`` for files ``<jobname>.e<jobID>`` (error and warnings log), ``<jobname>.o<jobID>`` (output log) .  Also check output of command ``qacct -j <jobID>`` .
- Copy output data files back using WinSCP, stfp, scp, or rsync .
