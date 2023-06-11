function w_hm=get_w_HM(F0,n,nionm5n0,Tm5n0,Mode,c,Za)
    if ~isfield(Mode,'kappa') ||~isfield(Mode.kappa,'hm_use_Potekhin') || ~Mode.kappa.hm_use_Potekhin
        w_hm=get_w_hm_original(F0,n,Mode.HoltsTable,c.a0,Za);
    else
        %Hydrogen only! Must adjust nion otherwise.
        w_hm=get_w_hm_Potekhin(F0,n,nionm5n0 , Tm5n0 , nionm5n0 ,Mode,c);
    end
end

function [w_hm] = get_w_hm_original(F0,n,HoltsTable,a0,Za)

Fcrit = Za^3/3 ./ (n+1).^2 ./ n.^3/a0^2;

Fcrit(imag(Fcrit)~=0)=0;
beta_crit = Fcrit./F0;

w_hm = HoltsTable.F(beta_crit);

w_hm(beta_crit>max(HoltsTable.beta)) = 1;
w_hm(beta_crit<min(HoltsTable.beta) | isinf(beta_crit)) = 0;
% w_hm(beta_crit<min(HoltsTable.beta)) = 0;
end

function [w_hm] = get_w_hm_Potekhin(F0,n,nionm5n0 , Tm5n0 , nem5n0 ,Mode,c)
%% Currently only specified for single species plasmas!
rion = rion/c.a0; % in atomic units

%% From Potekhin et. al. 
% Note that Znim5n0 = nem5n0;
Fcrit = (1)^3/3.*(n+1/2)./(n.^2.*(n+1).^2.*(n.^2+n+1/2));
Fcrit(imag(Fcrit)~=0)=0;
beta_crit = Fcrit./F0;

%% Code currently only working for Hydrogen!!!
a = (4*pi*nionm5n0/3).^(-1/3);
Gamma = (Mode.Plasma.Z(1)*c.q)^2./a./Tm5n0;
% clear chi k_s
% for i = 1:length(Tm5n0)
%     chi(i) = fzero(@(x) fermi_integral(1/2,x)-pi*pi*(c.hbar)^3*(c.me*Tm5n0(i)).^(-3/2).*nem5n0/sqrt(2),1);
%     k_s(i) = sqrt(c.q^2/pi/c.hbar^3*(2*c.me)^(3/2)*sqrt(Tm5n0(i))*fermi_integral(-1/2,chi));
% end
% 
% s = a*k_s;

% disp(['nim5n0 ' num2str(nionm5n0) ' T(eV) ' num2str(Tm5n0/c.eV) ' Gam ' num2str(Gamma) ' s ' num2str(s)]);

%% From Potekhin choose the unscreened s=0 fit for now. Eq. 17
alpha0=14.600; alpha1=103.20; alpha2=11.127; alpha3=16.178;
beta0=0.41; beta1=1.54; beta2=0.58; beta3=0.60;
gamma0=0.707; gamma1=1.64; gamma2=0.572; gamma3=0.915;

q0 = alpha0*(1+beta0*sqrt(Gamma)).^(-gamma0);
q1 = alpha1*(1+beta1*sqrt(Gamma)).^(-gamma1);
q2 = alpha2*(1+beta2*sqrt(Gamma)).^(-gamma2);
q3 = alpha3*(1+beta3*sqrt(Gamma)).^(-gamma3);

%Q_beta is the cumulative probability of an atom experiencing up to beta or
%less.
w_hm = ( q0.*beta_crit.^3-1.33.*beta_crit.^(9/2)+beta_crit.^6 )./( q1+q2.*beta_crit.^2+q3.*beta_crit.^3-beta_crit.^(9/2)/3+beta_crit.^6 );
end