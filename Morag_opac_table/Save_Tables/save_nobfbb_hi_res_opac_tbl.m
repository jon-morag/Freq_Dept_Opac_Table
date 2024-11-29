%% Best to first get Mode.temp from simulation prompt
path_settings;

clear Mode_temp;
Mode_temp.Plasma = Mode.Plasma;
Mode_temp.kappa.ff_on = 1;
Mode_temp.kappa.bf_on = 0;
Mode_temp.kappa.bb_on = 0;

[ kappa_abs_no_bf_bb,kappa_es,nu_calc ] = produce_high_res_opac_tbl(Mode_temp, N_nu );

save(OpacTableFilename , 'kappa_abs_no_bf_bb', 'kappa_es','nu_calc');