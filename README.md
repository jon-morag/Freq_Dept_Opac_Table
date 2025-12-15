# NuPac - A Frequency Dependent Opacity Table Code

A frequency-dependent opacity table, assuming LTE and arbitrary density, temperature, and chemical composition, with free-free, bound-free, and bound-bound components.
Written in Matlab (Octave not supported).

————————————————

Cite as:
DOI: 10.5281/zenodo.17935053
(release 2025 v2.1)

————————————————

How to Use:
- Download folder and cd into it (no additional setup required)
- Produce and plot an opacity for a specific (rho,T) and composition using 'produce_plot_opac.m'
- Produce tables for multiple (rho,T) or (R,T) using 'produce_save_hires_and_MG_tbls...', where R=T^3/rho.
   Can also produce low-res (MG) tables and employ the Sobolev approximation (Eastman, Pinto 1993)
- Limiting expansion opacity formalism is now available as well when producing opacities (see https://arxiv.org/pdf/2509.17003).
- Doppler shift and broaden opacity tables (approximately) using 'doppler_broaden_shift.m'
- 'standard_opac_profile.m' sets which physical components are included.

- To download and use additional line lists from Kurucz, see 'setup_old' directory (not supported)
- Optional: Install the 3rd party function ScaleTime for a ~25% efficiency improvement. Make sure to turn its use on in 'standard_opac_profile.m'.
ScaleTime available at https://www.mathworks.com/matlabcentral/fileexchange/25463-scaletime
