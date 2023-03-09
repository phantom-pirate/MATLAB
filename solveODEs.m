function ODEs
k1a=10*exp(4000*(1/300-1/T));
K2a=0.9*exp(9000*(1/300-1/T));
ct0=0.1;
T0=423;
Ft=Fa+Fb+Fc;
ca=ct0*(Fa/Ft)*(T0/T);
cb=ct0*(Fb/Ft)*(T0/T);
cc=ct0*(Fc/Ft)*(T0/T);
r1=-k1a*ca;
r2=-k2a*ca^2;
dFadV=r1+r2;
dFbdV=-r1;
dFcdV=-r2/2;
dTdV=(4000*(373-T)+(-r1)*20000+(-r2)*60000/(90*Fa+90*Fb+180*Fc));
end
