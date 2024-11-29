function [tbl_R, tbl_T] = set_RT_grid_30_120()
%% Current table RT settings
    tbl_R = logspace(-24 , -18 , 30);      
    tbl_T = 10.^(linspace( 0 , sqrt(2.7+0.6) , 120 ).^2)/10^0.6;
end