This is the code accompanying the paper "Differences among species in seed dispersal and conspecific neighbor effects can interact to influence coexistence".  It was all originally run on Matlab 2017a.  To run it, you will need to download the function shuffle.m, which can be found at https://www.mathworks.com/matlabcentral/fileexchange/27076-shuffle.  Otherwise, this contains everything you will need.

analyze_ML_data.m- This runs the analysis in appendix section 5.2, to calculate how seed dispersal limitation makes distance-responsive predators less effective.  It includes data from Muller-Landau et al. (2008) (full citation in main text).  This runs by itself, and displayed the outputs on screen.  
DE_shaken2.m- This function runs the difference equation version of the "shaken" approximation.  It is meant to be run on data generated by JC_invade1.  It also estimates how well the covariance approximations fit the data from JC_invade1.
dispersal_kern.m- This function is used by JC_invade1 and JC_holdAbund2.  It generates the dispersal kernels for a given plant species.
JC_holdAbund2.m- This runs the simulation.  It is similar to JC_invade1, however, instead of running an invasion analysis, it picks a focal species and hold's that species frequency steady.
JC_invade1.m- This is the main function running the simulations.  It runs an invasion analysis of the model in the main text.  It can also be used to run ordinary simulations.  There is a function for generating a movie file inside of it (see notes inside; you will need to alter one variable in the code to activate it).
make_neighbor.m- This code was modified from Stump et al. (2018, Journal of the Royal Society Interface).  It is used by JC_invade1 and JC_holdAbund2.  It makes a large matrix which lists all of the neighboring cells of each site.
manyCol.m- This generates 24 colors that I was easily distinguished (it has not been tested for people with color blindness).  It is used for drawing some of the figures.
ripleyK.m- This is used by JC_invade1 and JC_holdAbund2.  It puts a number on how aggregated each species is at the end of a simulation by calculating their pair correlation.
runAll.m- This runs all of the functions and creates all of the figures (except figures that were hand-drawn).  
run_approx_check.m- This runs a "shaken" simulation, along with the difference equation, to see how good the approximation is.
run_emergent_tradeoff.m- This creates a bunch of random communities, simulates them using JC_invade1, and tests if a dispersal-predation trade-off emerges naturally from community assembly.
run_for_distribution.m- This runs JC_holdAbund2 4 times, and saves the distribution of species.  It is used for Fig. 2b and 2c.
run_holdAbund.m- This runs JC_holdAbund2 several times.  It tracks the pair correlation function of each species, along with their growth rate, and the Delta mechanisms.
run_invasion.m- This runs JC_invade1 many times to run an invasion analysis.  It tracks the growth rate of the invader, the pair correlation of all species, and the strength of the mechanisms.  It also keeps track of the "shaken" DE approximation.

The functions "draw_figX" each draw one of the figures.

Run the following functions to generate each figure (or just run runAll.m for all of them):

Fig. 2:
run_for_distribution
run_holdAbund(11)

Fig. 3:
run_invasion(1,0)
run_invasion(2,0)
run_invasion(3,0)
run_invasion(4,0)
run_invasion(6,0)
run_invasion(7,0)

Fig. A1:
run_invasion(1,1)
run_invasion(2,1)
run_invasion(3,1)
run_invasion(4,1)
run_invasion(6,1)
run_invasion(7,1)

Fig. A2:
run_approx_check(2)
run_approx_check(3)
run_approx_check(6)

Fig. A3 and A4:
run_holdAbund(1)
run_holdAbund(2)
run_holdAbund(3)

Fig. A5:
run_emergent_tradeoff
