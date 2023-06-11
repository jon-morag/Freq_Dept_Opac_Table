function [km5g0n0]=kappa_abs_Brems(Tm5n0,rhom5n0,nu,Mode,Ym5n0,ResPop,~)
c=Mode.consts;
NG=length(nu);
% bremsstrahlung. 

% plasma properties
A=Mode.Plasma.A; % Atomic mass
Z=Mode.Plasma.Z; % Atomic number
% Y=Mode.Plasma.Y; % Degree of ionization (ne/np)

if length(A)>1
    Xfrac = Mode.Plasma.Xfrac;
else
    Xfrac = 1;
end

mec2=c.me*c.c^2;
hc=c.h*c.c;

% fac=4/3*pi^(-3/2)*(c.ale*hc)^3/mec2*(2/3/mec2)^(1/2)*(Z)^3/(A*c.amu)^2/c.h;

fac=4/3*pi^(-3/2)*(c.ale*hc)^3/mec2*(2/3/mec2)^(1/2)/c.amu^2/c.h;

% normalization by Z and A, through Y

ni_Zsq = 0;
for ix = 1:length(Xfrac)
    if Mode.Ionization.on
        jz = 1 : Mode.Plasma.Z(ix);
        ni_Zsq = ni_Zsq + Xfrac(ix)/A(ix).* ( jz.*jz * ResPop(ix).z(2:end,:) );
    else %Assume fully ionized
        ni_Zsq = ni_Zsq + Xfrac(ix)/A(ix).* Mode.Plasma.Z(ix)^2;
    end
end
ZAYm5n0 = Ym5n0*(Xfrac*(Z./A)') .* ni_Zsq;

% the emissivity is in units (cm^-3 s^-1) (as it is radiated power per unit
% photon energy) note that the emissivity is the radiated power in 4 pi
% steradians.

% nuM=repmat(nu,max(length(Tm5n0),length(rhom5n0)),1);
% if length(Tm5n0)>1
%     TMm5n0=repmat(Tm5n0',1,NG);
% else
%     TMm5n0=Tm5n0;
% end
% if length(rhom5n0)>1
%     rhoMm5n0=repmat(rhom5n0',1,NG);
% else
%     rhoMm5n0=rhom5n0;
% end   
% if length(Ym5n0)>1
%     YMm5n0=repmat(Ym5n0',1,NG);
% else
%     YMm5n0=Ym5n0;
% end 
xM=nu./Tm5n0';

% old Gaunt factor
% geff=3^(1/2)/pi*(log(2.25./xM).*(xM<1)+log(2.25).*(xM).^(-1/2).*(xM>=1));

% accurate gaunt factor in the Born approximation (Maxon 1972)
% geff=3^(1/2)/pi*cat(1,(log(4/1.78./xM(xM<=0.4)).*(1+3.52*(xM(xM<=0.4)/2/3.6).^2)+0.62*(xM(xM<=0.4)/4).^2).*exp(xM(xM<=0.4)/2), besselk(0,xM(xM>0.4&xM<6.96)/2).*exp(xM(xM>0.4&xM<6.96)/2), (2*1.57./xM(xM>=6.96)).^(1/2).*(1-2*0.125./xM(xM>=6.96)+4*0.07./xM(xM>=6.96).^2));
% geff=reshape(geff,size(xM));

% new and accurate Gaunt factor
geff=3^(1/2)/pi*besselk(0,xM/2,1); % matlab's function for besselk(0,x)*exp(x)

% geff=1;

fac = fac.*(c.h*c.c)^3/(8*pi)/c.c;
indlow = xM<1e-8;
kk = (fac.*ZAYm5n0.*rhom5n0./Tm5n0.^(0.5))'./nu.^2.*geff;
km5g0n0 = kk.*((1-exp(-xM))./nu.*(~indlow)+1./Tm5n0'.*(indlow));

km5g0n0(isnan(km5g0n0)==1)=0;

% km5g0n0=0*km5g0n0;
end