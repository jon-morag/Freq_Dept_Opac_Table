function [tbl_R , tbl_T] = set_RT_grid_16_66()
%% Current table RT settings
    tbl_R = logspace(-24 , -18 , 16);      
    tbl_T = 10.^(linspace( 0 , sqrt(2.7+0.6) , 66 ).^2)/10^0.6;
end