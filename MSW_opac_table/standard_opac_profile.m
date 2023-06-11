function [ Mode_profile ] = standard_opac_profile()

% Opac tables definitions
Mode_profile.consts = set_consts();

%% Ionization
Mode_profile.Ionization.on = 1;
Mode_profile.Ionization.Saha = 1;
Mode_profile.Ionization.use_w_hm = 1; % Hummer-Mihalas factor - only used in Saha before HnonHlike split
Mode_profile.Ionization.use_w_hm_Hlike = 1;
Mode_profile.Ionization.use_w_hm_nonHlike = 0;
Mode_profile.Ionization.levels_cutoff_frac = 1;

Mode_profile.Ionization.nmax_NIST = 1e8;
Mode_profile.Ionization.nmax = 100;
Mode_profile.Ionization.nmax_levels_Hydrogenic = 150;

%% Kappa main
Mode_profile.kappa.Hydrogen_bb_Nir = 1;
Mode_profile.kappa.Hydrogenic_bb_Nir = 1;
Mode_profile.kappa.bb_Kurucz = 1;
Mode_profile.kappa.low_res_spacing = 100;


Mode_profile.kappa.bb_species_select = 1:10;
Mode_profile.kappa.bb_Hydrogenic_species_select = 1:10;
Mode_profile.Plasma.bb_i_state_select = 0; % 0 means ignore,any other no. is selected state.
Mode_profile.kappa.bf_species_select = 1:10;

%% Bound free
Mode_profile.kappa.bf_Hydrogenic_on = 1;
Mode_profile.kappa.bf_VY_on = 1;

%% Bound Bound
Mode_profile.Hjerting_size = 'normal';
Mode_profile.kurucz_cd = 23;

Mode_profile.kappa.add_delta_func = 0;
% Mode_profile.kappa.n_FWHM_hi_res = 3;
Mode_profile.kappa.n_FWHM_hi_res = 100;
Mode_profile.kappa.n_FWHM_fcons = 2;
Mode_profile.kappa.binning_function = 'Optimal';

%% Lines - Levels - Extensions
% Mode.kappa.just_the_levels = 1; % For level counting
Mode_profile.kappa.nmax_lines_Kurucz = 1e8;
Mode_profile.kappa.extend_partition_func_n = 0;
Mode_profile.kappa.e_partition_frac_cutoff = 1e-3;
Mode_profile.kappa.extension_quantum_number = 500;
Mode_profile.kappa.levels_cutoff_frac = 1;
Mode_profile.Plasma.ion_pop_threshold = 1e-14;

%% Physics Extensions
Mode_profile.kappa.include_w_hm_Hydrogenic = 1; %Inside the kappa only. 
Mode_profile.kappa.include_w_hm_nonHydrogenic = 0; %Inside the kappa only.
% Mode_profile.kappa.include_w_hm_Kurucz = 0; %Inside the kappa only.
% Mode_profile.kappa.include_w_hm = 1; %Inside the kappa only. 
Mode_profile.kappa.hm_use_Potekhin = 0; 

Mode_profile.kappa.Hydrogenic_Dappen_HM_factor = 0; % A prescription for melding bf and bb contributions.
Mode_profile.kappa.Kurucz_Dappen_HM_factor = 0; %Currently didn't add bf contribution

Mode_profile.kappa.fix_A_decay_Kurucz = 1;
Mode_profile.kappa.add_omega4 = 1; %Seaton et. al. 1994
Mode_profile.kappa.H_e_impact_broadening = 0; % Armstrong et. al. 1994

%% Calculating / Plotting functions
Mode_profile.kappa.use_ScaleTime = 0;
Mode_profile.kappa.plot_FWHM = 0;
Mode_profile.kappa.plot_n_eff = 0;

end