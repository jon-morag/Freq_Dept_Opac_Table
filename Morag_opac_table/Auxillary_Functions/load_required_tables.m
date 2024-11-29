function [Mode_profile] = load_required_tables( Mode_profile)
% These tables can be produced using: produce_bb_Transition_Rate_matrix_from_Kurucz

if strcmp(Mode_profile.kurucz_cd,'23')
    load('Kurucz_lines_cd23','Transition_Rate_Matrix')
elseif strcmp(Mode_profile.kurucz_cd,'23_08oct17')
    load('Kurucz_lines_cd23_08oct17','Transition_Rate_Matrix')
else
    error('Unknown line data source.')
end

load('H_Fe_ion_energies_NIST','H_Fe_ion_energies');
load('H_Fe_levels_NIST','H_Fe_levels');
load('H_Transition_Rate_Matrix_data','H_Transition_Rate_Matrix');
load('bf_VY_tbl_short','bf_VY_tbl_short');

switch Mode_profile.Hjerting_size
    case 'small'
        load('HjertingTable_data_small.mat','HjertingTable_data_small');
        Mode_profile.HjertingTable = HjertingTable_data_small;
    case 'normal'
        load('HjertingTable_data.mat','HjertingTable_data');
        Mode_profile.HjertingTable = HjertingTable_data;
    case 'large'
        load('HjertingTable_data_large.mat','HjertingTable_data_large');
        Mode_profile.HjertingTable = HjertingTable_data_large;
end

Mode_profile.Transition_Rate_Matrix = Transition_Rate_Matrix;
Mode_profile.H_Transition_Rate_Matrix = H_Transition_Rate_Matrix;
Mode_profile.kappa.bf_VY_tbl = bf_VY_tbl_short;
Mode_profile.Ionization.levels = H_Fe_levels;
Mode_profile.Ionization.energies = H_Fe_ion_energies;
Mode_profile.HoltsTable=HoltsmarkTable();

end