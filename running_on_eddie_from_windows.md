- Copy the input files from local machine to Eddie using WinSCP, stfp, scp, or rsync.  Copy to a directory in ``/exports/eddie/scratch/<uun>``
  logseq.order-list-type:: number
- Open a Powershell command line and log into eddie: 
  logseq.order-list-type:: number
  ``ssh <uun>@eddie.ecdf.ed.ac.uk``
- Write/edit a submit script ``submit.sh`` such as this one (MPI) and copy it to eddie (to ``\home\<uun>\`` or project scratch directory).
  logseq.order-list-type:: number
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
	  
	  mpirun -np 2 ogs 3D                     #<--  change number and input (.pcs) file name
	  ```
	- more options: https://www.wiki.ed.ac.uk/display/ResearchServices/Job+Submission
- Navigate to directory of submit script and submit the job:  ``qsub submit.sh `` .
  logseq.order-list-type:: number
  Check progress with ``qstat`` :  ``qw`` for waiting, ``r`` for running, ``E`` for error/exitting.  Delete jobs with ``qdel <jobID>`` .
	- Troubleshooting:  Look in ``/exports/eddie/scratch/<uun>/<project>`` for files ``<jobname>.e<jobID>`` (error and warnings log), ``<jobname>.o<jobID>`` (output log) .  Also check output of command ``qacct -j <jobID>`` .
- Copy output data files back using WinSCP, stfp, scp, or rsync .
  logseq.order-list-type:: number