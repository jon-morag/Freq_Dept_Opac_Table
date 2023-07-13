%% Grid
rho = 1e-11; %g/cc
T = 1; %eV
N_nu = 1e6; %number of frequencies
should_plot = 1;

%% Commonly changed settings
include_ff = 1;
include_bf = 1;
include_bb = 1;

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
   
%% Produce
tic
[ kappa_abs,kappa_es,nu_calc ] = produce_hires_opac_tbl_rhoT( N_nu , rho , T,include_ff , include_bf ,include_bb, A , Z , Xfrac ); %Z_metal=[]->do nothing
toc

%% Plot
if should_plot
    c=set_consts();
    disp(['kappa_es=' num2str(kappa_es)])
    
    figure
    ax=gca(); ax.ColorOrderIndex = 1; ax.FontSize = 14;
    ax.XScale = 'log'; ax.YScale = 'log';

    loglog(nu_calc/c.eV , kappa_abs , 'LineWidth',1.5);
    hold on
    xlabel('h\nu [eV]')
    ylabel('\kappa_\nu [cm^2 g^{-1}]')
end
