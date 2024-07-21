# NuPac - A Frequency Dependent Opacity Table Code

A frequency-dependent opacity table, assuming LTE and arbitrary density, temperature, and chemical composition, with free-free, bound-free, and bound-bound components.
Written in Matlab (Octave not supported).

————————————————

cite as: https://arxiv.org/abs/2307.05598

————————————————

How to Use:
- Download folder and run 'setup.m' once
- Produce and plot an opacity for a specific (rho,T) and composition using 'produce_plot_opac.m'
- Produce tables for multiple (rho,T) or (R,T) using 'produce_and_save_hires_opac_tbl….m', where R=T^3/rho.
- 'standard_opac_profile.m' sets which physical components are included.
- Optional: Install the 3rd party function ScaleTime for a ~25% efficiency improvement. Make sure to turn its use on in 'standard_opac_profile.m'.
ScaleTime available at https://www.mathworks.com/matlabcentral/fileexchange/25463-scaletime
