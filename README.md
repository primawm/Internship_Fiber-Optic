# Internship_Fiber-Optic
Superficial Structure of Stromboli Volcano using Fibre Optic Measurement

SCRIPTS DESCRIPTION

1. plot_imshow_channel_offset :
- Plotting seismic trace of each DAS channel vs time

2. check_availibility :
- Plot a diagram showing when the data is present
  
3. stack_all :
- Cross-correlate DAS strain data using phase cross-correlation (PCC) or phase auto-correlation (aPCC)
- Stack cross-correlogram using time-frequency phase weighted stack (tf-PWS)

4. stack_functions :
- Functions used in stack_all.py

5. similarity :
- To plot the effect of number of data towards the convergence of stack results

6. plot_station :
- To plot the station location

7. fk_transform_and_phase_velocity :
- Make a plot of correlogram gather
- Make FK-transform over several cross-correlograms/traces
- Pick highest amplitude manually if needed
- If highest amplitude is picked, will plot phase velocity

8. stockwell_transform_and_group_velocity
- Compute stockwell transform of a correlogram
- Auto-pick the group velocity
