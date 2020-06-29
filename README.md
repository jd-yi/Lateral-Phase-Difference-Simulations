# Lateral-Phase-Difference-Simulations
Contains the simulation scripts and data necessary to reproduce all figures in the manuscript: "Lateral cortical phase differences are sufficient for the development rhythmic spatial sampling," by Justin D. Yi and  Katsushi Arisaka

To reproduce the following figures, you may run the approriate `*.m` file from the list below:

**Fig. 2**: TRAVEL_ONLY_WHEN_STIM_two_area_model_iFFT_granger_allplot_ANNOTATED.m. To adjust the phase difference, edit the variable `shift` (line 66)


**Fig. 3A/B**:attractor_dynamics_phase_bias_testing_ANNOTATED.m To adjust the phase difference, edit the variable `shift` (line 67). To edit the amplitude of the stimulation, edit the coef. in front of the variable `stimA` (lines 85 and 86 for Regions 1 and 2 respectively).


**Fig. 3C**:attractor_dynamics_phase_bias_ANNOTATED.m. Like above, `shift` is on line 74, ampl (the amplitude of the LFP) is on line 85, and stimA coefs. are on 94 and 95. 


**Fig. 4/S2**:attractor_dynamics_phase_bias_presentatphase_ANNOTATED.m. The aforementioned variables are all editable. However, the times at which to poll the state (is E3>E4?) is controlled by `tsamps` on line 57. 


**Fig. S1**:ifft_alpha_ANNOTATED.m

*Disclaimer: Only tested and ran on MATLAB version R2019b.*
