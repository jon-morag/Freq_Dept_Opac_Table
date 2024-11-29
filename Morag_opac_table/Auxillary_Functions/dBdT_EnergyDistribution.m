%% Derivative of the planck spectrum - (dimensionsless): energy density per unit frequency (in energy)
function [dBdTm5g0n0]=dBdT_EnergyDistribution(Tm5n0,nu)

c=set_consts;
hc=c.h*c.c;

% nuM=repmat(nu,length(Tm5n0),1);
% TM=repmat(Tm5n0',1,length(nu));

x=(nu./Tm5n0);
indhigh=x>300;
indlow=x<1e-8;

hc3=power_int(hc,3);
dBdT1m5g0n0=8*pi/hc3*(nu).^4./Tm5n0.^2;
dBdTm5g0n0=dBdT1m5g0n0.*((indhigh).*exp(-x)+(~indlow&~indhigh).*exp(x)./(exp(x)-1).^2+(indlow).*exp(x).*(Tm5n0./nu).^2);

dBdTm5g0n0=max(dBdTm5g0n0,1e-300);