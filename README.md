# Freq_Dept_Opac_Table

A frequency-dependent opacity table, assuming LTE and arbitrary density, temperature, and chemical composition, with free-free, bound-free, and bound-bound components.
Written in Matlab.

————————————————

How to Use:
- Run 'setup.m' once
- Produce and plot an opacity for a specific (rho,T) and composition using 'produce_plot_opac.m'
- Produce tables for several (rho, T) or (R,T) using 'produce_and_save_hires_opac_tbl….m'
- 'standard_opac_profile.m' sets which physical components are included.
- Install the 3rd party function ScaleTime for a ~25% efficiency improvement. Make sure to turn its use on in 'standard_opac_profile.m'.
ScaleTime available at https://www.mathworks.com/matlabcentral/fileexchange/25463-scaletime

————————————————

Physical Details:
