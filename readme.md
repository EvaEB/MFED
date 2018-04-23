# ABC-SMC fit of the MFED of HIV

## Files in this archive
* `analyze_SMC.py` - plots and calculates model probabilities and parameter estimates
* `getSummaryStats.py` - calculates summary statistics on the output of a simulation
* `initial_smc.py` - generates random parameter sets according to prior distribution
* `next_par.py` - generates new parameter sets based on accepted parameter sets from the previous iteration
* `patients` - text file containing the patient characteristics to match in the simulations. These were created from sequence sets (.fasta) available from Keele et al (2008) and Li et al (2010).
* `simulate_all.py` - created simulated dataset for a large number of parameter sets
* `simulate.py` - simulate a single single patient
* `simulation_functions.py` - all functions needed for the simulations


## Running ABC-SMC fit of MFED in early HIV infection
1. generate initial parameters: `python initial_smc.py [parameter_filename]`
this will generate 700 parameter sets. Call the script multiple times for more
parameter sets (that can be run in parallel)

2. perform simulation for all parameter sets:
`python simulate_all.py [parameter_filename] [simulation_directory] [outfile]`

 In order for the next script to work, outfile should be in the format
 `stats_set_[iter]/stats_[nr].npy`

 NOTE: these simulations are expensive. The simulation for a single parameter
 set takes ~25 seconds on a regular computer, so a parameter file as generated
 in step one will take ~5 hours to complete. It is however safe to stop the
 script while it is running, since it will save to [outfile] after every
 simulated parameter set and restart with the right parameter set.

3. generate new parameter sets based on the results of the simulations:
`python next_par.py [max_dist] [stats_dir] [new_par_dir] [name]`

  We did 7 iterations, with max_dist equal to 5, 2.2, 1.3, 0.8, 0.7, 0.6, 0.5

4. repeat steps 2 & 3, lowering [max_dist] as required.
5. analyze the final iteration: `analyze_SMC.py`
  this will generate a plot with the model probabilities per iteration, and
  plots of the posterior density of all the parameters for those models with at
  least 100 accepted parameter sets
