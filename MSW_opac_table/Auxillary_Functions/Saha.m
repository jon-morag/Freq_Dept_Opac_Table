function [ne,Y,ResPop,Zeff]=Saha(T,rho,Mode)
% [ne,Y,ResPop,Zeff]=SahaMix(T,rho,Z,A,X)
% Saha equation for a mixture descirbed by mass fractions X, atomic charge
% Z and atomic weight A, including Hummer Mihalas factor. A table describing the energy levels and the
% degeneracy for each level should be supplied. The function produces the
% total electron density ne, and the fractions of each ionized
% specie (non-ionized, single-ionized, doubly-ionized, etc.). T and
% rho are the temperature and the density and should be vectors of the
% same size. 
% This function uses a Newton-Raphson solver for the non-linear set of
% equations.
% example for usage:
% X=[0.7,0.3];
% Z=[1,2];
% A=[1,4];
% A=[1.008,4.003];
% T=1.5*c.eV;
% rho=1e-13;
% [ne,Y,ResPop,Zeff] = Saha(T,rho,Mode);
c=set_consts();

X=Mode.Plasma.Xfrac;
Z=Mode.Plasma.Z;
A=Mode.Plasma.A;

nmax = Mode.Ionization.nmax_NIST; %Number of excitation levels
ion_levels = Mode.Ionization.levels;
ion_energies =  Mode.Ionization.energies;

% initial guess for the electron density :
%near full ionization if above 1 eV. very weak ionization for low T
nef_guess = (T/c.eV).^4;
nef_guess(nef_guess>0.97) = 0.97;
ne = nef_guess.*sum(X./A.*Z).*rho/c.amu;

% the number density of each element
natom = rho*sum(X./A)/c.amu;
rion = (4*pi*natom/3).^(-1/3);

iter = 1;
non_conv_cond_vec = true(length(T),1);

%% Find ne such that Saha is satisfied
convAbsY = 1e-7; convTolX = 1e-7; convTolY = 1e-7;
% convAbsY = 1e-10; convTolX = 1e-10; convTolY = 1e-10;
while any(non_conv_cond_vec) && iter<1000
    ixtocalc = find(non_conv_cond_vec);
    T_calc = T(ixtocalc);
    ne_calc = ne(ixtocalc);
    rho_calc = rho(ixtocalc);
    sum_z_r = zeros(length(X),length(ixtocalc));
    dsum_z_rdne = zeros(length(X),length(ixtocalc));
    
    for ix = 1:length(X)
        %Partition function
        g_r = zeros(Z(ix),length(rho));
        F0 = ne_calc./rho_calc*c.amu/sum(X./A)./rion.^2;
        for ij = 1:Z(ix)
            if ij == Z(ix)
                n = 1:Mode.Ionization.nmax_levels_Hydrogenic;
                g = 2*n.^2;
                E = c.IH*Z(ix)^2*(1-1./n.^2);
                if Mode.Ionization.use_w_hm_Hlike
                    w_hm = get_w_HM(F0,n,1,1,Mode,c,ij);
                    g_r(ij,:) = sum(g.*w_hm.*exp(-E ./ T_calc));
                else
                    g_r(ij,:) = sum(g.*exp(-E ./ T_calc));
                end
            elseif Z(ix) > ij
%             if true
                g = ion_levels(Z(ix),ij).g;
                E = ion_levels(Z(ix),ij).energies;
                levels_cutoff_length = ceil(length(g)*Mode.Ionization.levels_cutoff_frac);
                g = g(1:levels_cutoff_length);
                E = E(1:levels_cutoff_length);
           
                if Mode.Ionization.use_w_hm_nonHlike
                    n_eff = ij*real(sqrt(c.IH/c.eV./(ion_energies(Z(ix),ij)-E)));
                    w_hm = get_w_HM(F0,n_eff,1,1,Mode,c,ij);
                    w_hm(n_eff==0) = 0;
                    
                    g_r(ij,:) = sum(g(1:min(nmax,length(g))).*w_hm(1:min(nmax,length(g)),:) .* exp(-E(1:min(nmax,length(g)))*c.eV./T_calc))';
                    %                 n_eff = cellfun(@(I,Z_ion,E) Z_ion*sqrt(c.IH/c.eV./(I-E)) , num2cell(ion_energies(Z(ix),1:Z(ix))) , num2cell(1:Z(ix)) , {ion_levels(Z(ix),1:Z(ix)).energies},'UniformOutput',false);
                    %                 w_hm = cellfun(@(n,Z) get_w_HM(F0,n,1,1,Mode,c,Z) , n_eff , num2cell(1:Z(ix)),'UniformOutput',false );
                    %                 g_r = cell2mat(cellfun(@(g,w,E) sum(g(1:min(nmax,length(g))).*w(1:min(nmax,length(g))) .* exp(-E(1:min(nmax,length(g)))*c.eV./Tm5n0T))',{ion_levels(Z(ix),1:Z(ix)).g},w_hm,{ion_levels(Z(ix),1:Z(ix)).energies},'UniformOutput',false))';
                else
                    g_r(ij,:) = sum(g(1:min(nmax,length(g))) .* exp(-E(1:min(nmax,length(g)))*c.eV./T_calc))';
                    %The cell format does not allow us to vary the
                    %             g_r = cell2mat(cellfun(@(g,E) sum(g(1:min(nmax,length(g))).* exp(-E(1:min(nmax,length(g)))*c.eV./Tm5n0T))',{ion_levels(Z(ix),1:Z(ix)).g},{ion_levels(Z(ix),1:Z(ix)).energies},'UniformOutput',false))';
                end
            end
        end
        
        g_r = [g_r ; ones(1,length(T_calc))];
        I_r = ion_energies(Z(ix),1:Z(ix))'*c.eV;
        %Helium example
        %             g_r = [1,2,1];
        %             I_r = [24.580*c.eV, 0.5*Z(ix).^2*c.ale^2*c.me*c.c^2];
        
        theta_r = 2*g_r(2:end,:)./g_r(1:end-1,:).*(2*pi*c.me*T_calc).^(3/2)/c.h^3.*exp(-I_r./T_calc)./ne;
        
        lam_r = cumprod(theta_r,1);
        
        dlam_rddne = -(1:Z(ix))'.*lam_r./ne_calc;
        
        z_r = lam_r./(1+sum(lam_r,1));
        
        dz_rdne = dlam_rddne./(1+sum(lam_r,1))-z_r./(1+sum(lam_r,1)).*sum(dlam_rddne,1);
        
        sum_z_r(ix,:) = (1:Z(ix))*z_r;
        dsum_z_rdne(ix,:) = (1:Z(ix))*dz_rdne;
    end

    func = ne_calc./(rho_calc/c.amu)-(X./A)*sum_z_r;
    dfunc = 1./(rho_calc/c.amu)-(X./A)*dsum_z_rdne;
    
    dnem5n0 = func./dfunc;
    ne(ixtocalc) = ne_calc-dnem5n0;
    
    non_conv_cond_vec(ixtocalc) = (max(abs(func./(ne_calc./(rho_calc/c.amu))))>convAbsY)&& ...
        (max(abs(dnem5n0./ne_calc))>convTolX);
    
    iter = iter+1;

end

%% Use derived ne to get all electron levels
n0 = (X./A)'*rho/c.amu;
for ix = 1:length(X)
    %Partition function
    g_r = zeros(Z(ix),1);
    
    F0 = ne./rho*c.amu/sum(X./A)./rion.^2;
    for ij = 1:Z(ix)
        if ij == Z(ix)
            n = 1:Mode.Ionization.nmax_levels_Hydrogenic;
            g = 2*n.^2;
            E = c.IH*Z(ix)^2*(1-1./n.^2);
            if Mode.Ionization.use_w_hm_Hlike
                w_hm = get_w_HM(F0,n,1,1,Mode,c,ij);
                g_r(ij,:) = sum(g.*w_hm.*exp(-E ./ T_calc));
            else
                g_r(ij,:) = sum(g.*exp(-E ./ T_calc));
            end
        elseif Z(ix) > ij
%         if true
            g = ion_levels(Z(ix),ij).g;
            E = ion_levels(Z(ix),ij).energies;
            levels_cutoff_length = ceil(length(g)*Mode.Ionization.levels_cutoff_frac);
            g = g(1:levels_cutoff_length);
            E = E(1:levels_cutoff_length);
            
            if Mode.Ionization.use_w_hm_nonHlike
                n_eff = ij*real(sqrt(c.IH/c.eV./(ion_energies(Z(ix),ij)-E)));
                w_hm = get_w_HM(F0,n_eff,1,1,Mode,c,ij);
                w_hm(n_eff==0) = 0;
                
                g_r(ij,:) = sum(g(1:min(nmax,length(g))).*w_hm(1:min(nmax,length(g)),:) .* exp(-E(1:min(nmax,length(g)))*c.eV./T_calc))';
                %                 n_eff = cellfun(@(I,Z_ion,E) Z_ion*sqrt(c.IH/c.eV./(I-E)) , num2cell(ion_energies(Z(ix),1:Z(ix))) , num2cell(1:Z(ix)) , {ion_levels(Z(ix),1:Z(ix)).energies},'UniformOutput',false);
                %                 w_hm = cellfun(@(n,Z) get_w_hummermihalas(F0,n,1,1,Mode,c,Z) , n_eff , num2cell(1:Z(ix)),'UniformOutput',false );
                %                 g_r = cell2mat(cellfun(@(g,w,E) sum(g(1:min(nmax,length(g))).*w(1:min(nmax,length(g))) .* exp(-E(1:min(nmax,length(g)))*c.eV./Tm5n0T))',{ion_levels(Z(ix),1:Z(ix)).g},w_hm,{ion_levels(Z(ix),1:Z(ix)).energies},'UniformOutput',false))';
            else
                g_r(ij,:) = sum(g(1:min(nmax,length(g))) .* exp(-E(1:min(nmax,length(g)))*c.eV./T_calc))';
                %The cell format does not allow us to vary the
                %             g_r = cell2mat(cellfun(@(g,E) sum(g(1:min(nmax,length(g))).* exp(-E(1:min(nmax,length(g)))*c.eV./Tm5n0T))',{ion_levels(Z(ix),1:Z(ix)).g},{ion_levels(Z(ix),1:Z(ix)).energies},'UniformOutput',false))';
            end
        end
    end
    
    g_r = [g_r ; ones(1,length(T_calc))];
    I_r = ion_energies(Z(ix),1:Z(ix))'*c.eV;
    
    theta_r = 2*g_r(2:end,:)./g_r(1:end-1,:).*(2*pi*c.me*T).^(3/2)/c.h^3.*exp(-I_r./T)./ne;
    
    lam_r = cumprod(theta_r,1);
    
    ResPop(ix).n0 = n0(ix,:);
    ResPop(ix).z = [ones(1,size(lam_r,2)) ; lam_r]./(1+sum(lam_r,1));
end
Y = ne./sum((X.*Z./A)'*rho/c.amu,1);

ni = sum((X./A)'*rho/c.amu);
Zeff = ne./ni; % this is what appears in TOPS tables
% ci = (X./A)'/c.amu.*rho./ni; number fractions