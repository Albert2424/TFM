run:
	sbatch runwe.slurm

clean:
	rm -f *.h5 *.log *.err *.out
	rm -r traj_segs seg_logs istates
clust:
	@squeue -o "%.10i %.9P %.8j %.8T %.6M %.10l %.6D"
cancel:
	scancel -n WESTPA
	
help:
	@echo -e '################# HELP ##################'
	@echo ''
	@echo -e 'run: initializes and runs the WESTPA simulation.'
	@echo -e 'clust: provides information of the slurm queue.'
	@echo -e 'clean: erases all files produced on the previous simulation.'
	@echo -e 'cancel: cancels the job.'
	@echo ''
