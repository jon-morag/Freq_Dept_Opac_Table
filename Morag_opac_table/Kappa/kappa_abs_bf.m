function k_bf = kappa_abs_bf(n, T, nu, Mode,e_pop,Z,num_frac,F0,c)
k_bf = 0;
for ix =  Mode.kappa.bf_species_select % select specific species
    if Mode.kappa.plot_n_eff && Z(ix)>=2; figure(ix) ; end
    for ij = 1:Z(ix)
        if Z(ix) == ij && Mode.kappa.bf_Hydrogenic_on %H-like
            k_bf = k_bf + kappa_abs_bf_ij_Hlike(n, T, nu, Mode,e_pop(ix),Z(ix),num_frac(ix),F0,c);
        elseif Z(ix) > ij && Mode.kappa.bf_VY_on
            k_bf = k_bf + kappa_abs_bf_ij_VY(T,nu, Mode,e_pop(ix),Z(ix),num_frac(ix),F0,c,ij);
        end
    end
    if Mode.kappa.plot_n_eff && Z(ix)>=2
        ax = gca(); ax.FontSize = 14;
        xlabel('State Number'); ylabel('State Energy');
        title(['NIST listing vs. Effective "H-like" energy - Z = ' num2str(Z(ix))]);
        legend(p_En,cellfun(@(ij) ['Ionization ' num2str(ij)] , num2cell(2:Z(ix)),'UniformOutput',0))
    end
end
end

function k_bf_ix_Hlike = kappa_abs_bf_ij_Hlike(n, T, nu, Mode,e_pop_i , Z_i , num_frac_i,F0,c)
IH=c.IH;
gn_lo = 2*n.^2;

En = IH*Z_i^2./n.^2;
if Mode.kappa.include_w_hm_Hydrogenic
    w_hm = get_w_HM(F0,n,1,1,Mode,c,Z_i);
    Zr = gn_lo.*w_hm.*exp(-(IH*Z_i^2-En)./T');
else
    Zr = gn_lo.*exp(-(IH*Z_i^2-En)./T');
end
Gr = sum(Zr,2);
n_n = Zr./Gr;
nu_m_En = nu-En'; %+Icm5n0';
En_o_nu = En'./nu;
ind_allowed = (nu_m_En>0);

if Mode.kappa.Hydrogenic_Dappen_HM_factor
    % The equivalent n state
    n_DAM = real(( 1./n.^2 - nu/IH/Z_i.^2 ).^(-1/2));
    w_hm_DAM_bf = get_w_HM(F0,n_DAM,1,1,Mode,c,Z_i);

    nu_I = (nu_m_En).*(ind_allowed)+1.*(~ind_allowed);
    z = sqrt(En'./nu_I).*(ind_allowed)+1.*(~ind_allowed);
    g_nf = 8*pi*sqrt(3)*En_o_nu.*exp(-4*z.*acot(z))./(1-exp(-2*pi*z));
    sig_bf = n'./Z_i^2*8/sqrt(3)*c.ale^-3*c.sigT.*En_o_nu.^3;
    sig_bf = sig_bf .* (g_nf.*(nu_m_En>0)+8*pi*sqrt(3)*En_o_nu*exp(-4).*(nu_m_En<=0)); %This is a cheat, I use the gaunt factor value at zero
    sig_bf = sig_bf.*((1-w_hm_DAM_bf).*(nu_m_En<0) + 1.*(nu_m_En>=0));
else
    nu_I = (nu_m_En).*(ind_allowed)+1.*(~ind_allowed);
    z = sqrt(En'./nu_I).*(ind_allowed)+1.*(~ind_allowed);
    g_nf = 8*pi*sqrt(3)*En_o_nu.*exp(-4*z.*acot(z))./(1-exp(-2*pi*z));
    sig_bf = n'./Z_i^2*8/sqrt(3)*c.ale^-3*c.sigT.*En_o_nu.^3.*(g_nf.*(nu_m_En>0)+8*pi*sqrt(3)*En_o_nu*exp(-4).*(nu_m_En==0)+0.*(nu_m_En<0));
end

x = nu./T;
indlow = x<1e-8;

k_bf_ix_Hlike = e_pop_i.z(Z_i)*sum(n_n'.*sig_bf.*((1-exp(-x)).*(~indlow)+x.*(indlow)))*num_frac_i/c.amu;
end

function k_bf_ix_VY = kappa_abs_bf_ij_VY(T,nu, Mode,e_pop_i,Z_i,num_frac_i,F0,c,ij)
IH=c.IH;

[sigma_VY_bf_arr, n_VY] = sigma_VY_bf(nu,Z_i,ij,Mode.kappa.bf_VY_tbl,c);
if ~isempty(sigma_VY_bf_arr)
    In = Mode.Ionization.energies(Z_i,ij)*c.eV;
    gn = Mode.Ionization.levels(Z_i,ij).g; En=Mode.Ionization.levels(Z_i,ij).energies*c.eV;

    % Cut the levels by half to check if we are far from convergence
    levels_cutoff_length = ceil(length(gn)*Mode.kappa.levels_cutoff_frac);
    gn = gn(1:levels_cutoff_length);
    En = En(1:levels_cutoff_length);

    if Mode.kappa.include_w_hm_nonHydrogenic
        n_eff= real(ij*sqrt(IH./(In-En)));
        w_hm = get_w_HM(F0,n_eff,1,1,Mode,c,ij);
        Zr = gn.*w_hm.*exp(-En./T');
        Zr(find(n_eff==0)) = 0; %this means it was an imaginary result.
    else
        Zr = gn.*exp(-En./T');
    end

    Gr = sum(Zr,1);

    if Mode.kappa.extend_partition_func_n % && Zr(find(Zr,1,'last'))/Gr > Mode.kappa.e_partition_frac_cutoff %Turns out that this criterion is not always right. So always extend the partition function.
        if ~Mode.kappa.include_w_hm_nonHydrogenic; n_eff = ij * real(sqrt(IH./(In-En))); end % Don't calculate twice
        n_extend = round(max(n_eff)) + 1 : Mode.kappa.extension_quantum_number;
        %                         n_check = 1:round(max(n_eff));

        % For consistent addition assume the NIST ground state is zero.
        En_extend = In - ij^2*c.IH./n_extend.^2;

        if Mode.kappa.include_w_hm_nonHydrogenic
            w_hm_extend = get_w_HM(F0,n_extend,1,1,Mode,c,ij);
            Gr_extend = sum(2*n_extend.^2.*w_hm_extend.*exp(-En_extend./T'));
        else
            Gr_extend = sum(2*n_extend.^2.*exp(-En_extend./T'));
        end

        Gr = Gr + Gr_extend;
    end

    n_ground = sum(Zr(En<=0.3*c.eV))/Gr;
    %                     n_n = Zr./Gr;
    %All cross sections are given in terms of the ground
    %states.
    k_bf_ix_VY = e_pop_i.z(ij)*n_ground*sum(sigma_VY_bf_arr,1)*num_frac_i/c.amu;  

    if Mode.kappa.plot_n_eff
        % Examine the n_eff prescription
        p_En(ij) = plot(En/c.eV);
        hold on
        ax = gca(); ax.ColorOrderIndex = mod(ax.ColorOrderIndex-2,7)+1;
        plot(In*(1-1./round(n_eff).^2)/c.eV,'--')
        ax = gca(); ax.ColorOrderIndex = mod(ax.ColorOrderIndex-2,7)+1;
        len_n_eff = length(n_eff(n_eff>0));
        plot([len_n_eff len_n_eff+10] , [In/c.eV In/c.eV] ,'LineWidth',2.5)
    end
end

end