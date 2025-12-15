%% Planck spectrum - (dimensionsless): energy density per unit frequency (in energy)
function [Bm5g0n0]=PlanckEnergyDistribution(Tm5n0,nu)

c=set_consts;
hc=c.h*c.c;

x=(nu./Tm5n0);
indhigh=x>300;
indlow=x<1e-8;

hc3=hc^3;
B1m5g0n0=8*pi/hc3*(nu).^2;
Bm5g0n0=B1m5g0n0.*((indhigh).*nu.*exp(-x)+(~indlow&~indhigh).*nu./(exp(x)-1)+(indlow).*Tm5n0);

Bm5g0n0=max(Bm5g0n0,1e-300);