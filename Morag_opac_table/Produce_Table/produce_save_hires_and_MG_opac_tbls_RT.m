% Produce and save a high resolution opacity table binned in nu, T and R=rho/T^3
c=set_consts();

%% Commonly Changed Settings
% [tbl_R, tbl_T]=set_RT_grid_16_66;
% [tbl_R, tbl_T]=set_RT_grid_30_120;
[tbl_R,tbl_T] = deal([1e-20 1e-21],[1 4]); % test

save_dir = './Output_Tables/';
include_ff = 1;
include_bf = 1;
include_bb = 1;

sobolev = 0;% Produces a MG approximation for bb, based on EP 93 instead of hi-res
t_sobolev = 1*c.day;

make_MG = 1;

%% Nu Grid
if ~sobolev
    N_nu = 1e6; %Number of frequencies in grid
else
    N_nu = 100; %Number of frequency bins
end
Nnu_MG = 100;

%% Composition (comment out to choose or indicate aribtrary mixture, Xfrac by mass)
     
% Solar (choose metallicity below)
Xfrac_Sun = [0.7376 0.2494 0.0023865 0.0006902 0.005722 0.0012510 0.0007147 0.0006697 0.0003058 0.0012872];
A = [1.0079 4.0026 12.0107 14.0067 15.9994 20.1797 24.305 28.0855 32.065 55.845];
Z = [1 2 6 7 8 10 12 14 16 26];
Z_sun = sum(Xfrac_Sun(3:end));

mixname = 'Solar'; Z_metal = Z_sun;
% mixname = 'Solar0_1Z'; Z_metal = 0.1*Z_sun;
% mixname = 'SolarZ0'; Z_metal = 0;

% Adjust Xfrac to metalicity choice
Xfrac = Xfrac_Sun;
Z_sun = sum(Xfrac(3:end));
Xfrac(3:end) = Xfrac(3:end)*Z_metal/Z_sun;
Xfrac(1:2) = Xfrac(1:2) * (1-Z_metal)/(1-Z_sun);

     %%%%%%%%%%%%%%
     
% % C/O
% A = [12.0107 15.9994];
% nofrac = [0.5 0.5];
% Xfrac = nofrac.*A./(nofrac*A');
% Z = [6 8];
% mixname = 'CO';

     %%%%%%%%%%%%%%

% % He-C/O
% A = [4.0026 12.0107 15.9994];
% COnofrac = [0.5 0.5];
% Z_CO = 0.1;
% Xfrac = [1-Z_CO Z_CO*COnofrac.*A(2:3)./(COnofrac*A(2:3)')];
% Z = [2 6 8];
% mixname = ['HeCO' num2str(Z_CO)];

%% Saved Table name
if sobolev
    OpacTableFilename = [save_dir  mixname 'SobolevOpacTableR' num2str(length(tbl_R)) 'T' num2str(length(tbl_T)) 'Nnu' num2str(N_nu) 't_exp' num2str(t_sobolev/c.day) 'days.mat'];
else
    OpacTableFilename = [save_dir  mixname 'HiResOpacTableR' num2str(length(tbl_R)) 'T' num2str(length(tbl_T)) 'Nnu1e' num2str(log10(N_nu),2) '.mat'];
end
   
%% Produce, Save
tic
[ kappa_abs,kappa_es,nu_hi_res ] = produce_hires_opac_tbl_RT( N_nu , tbl_R , tbl_T,include_ff , include_bf ,include_bb, sobolev, t_sobolev, A , Z , Xfrac ); %Z_metal=[]->do nothing
toc

save(OpacTableFilename , 'kappa_abs', 'kappa_es','nu_hi_res','tbl_R','tbl_T','-v7.3');

%% Make accompanying MG table
if make_MG && ~ sobolev
    disp('Extracting MG table')
    [MG_tbl] = extract_MG_tbl_from_hi_res_kappa_RT(Nnu_MG,nu_hi_res,kappa_abs,kappa_es,tbl_R,tbl_T);

    MGTableFilename = [save_dir  mixname 'MGTableR' num2str(length(tbl_R)) 'T' num2str(length(tbl_T)) 'NG' num2str(Nnu_MG) '_from_Nnu1e' num2str(log10(N_nu),2) '.mat'];
    save(MGTableFilename,'MG_tbl','-v7.3')
end