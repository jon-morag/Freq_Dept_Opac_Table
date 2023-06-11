function [tbl_rho , tbl_T] = set_rhoT_grid_23_66()
%% Current table settings
    tbl_rho = 10.^(-16:0.5:-5);
    tbl_T = 10.^(linspace( 0 , sqrt(2.7+0.6) , 66 ).^2)/10^0.6;
end