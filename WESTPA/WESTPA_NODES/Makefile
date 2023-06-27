run:
	./set_input.sh
	sbatch -J $(shell basename "`pwd`") runwe.slurm
clean:
	rm -f *.h5 *.log *.err *.out *.pdf
	rm -r traj_segs seg_logs istates
new:
	clean
	run
clust:
	@squeue -o "%.10i %.9P %.8j %.8T %.6M %.10l %.6D %16S"
cancel:
	scancel -n $(shell basename "`pwd`")
analyse:
	wedap -h5 west.h5 -dt evolution --xlabel CS --bins 36 -o hist_evo.pdf
	wedap -h5 west.h5 -dt average -pm line --xlabel CS -o hist_av.pdf
	
help:
	@echo -e '################# HELP ##################'
	@echo ''
	@echo -e '    run : (initializes) runs the WESTPA simulation.'
	@echo -e '  clust : provides information of the slurm queue.'
	@echo -e '  clean : erases all files produced on the previous simulation.'
	@echo -e ' cancel : cancels the job.'
	@echo -e 'analyse : provides 1D energy plot and average plot.'
	@echo -e '    new : initializes and runs the WESTPA simulation.'
	@echo ''
