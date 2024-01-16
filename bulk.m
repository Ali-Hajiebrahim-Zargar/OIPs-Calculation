function [energy_level,H]=bulk(kx,ky,kz,a,Esa,Esc,Epa,Epc,Essa,Essc,Esasc,Esaxc,Exasc,Essaxc,Exassc,Exaxc,Exayc,deltaa,deltac)


k=[kx ky kz];


T1=(a/4)*[1 1 1];
T2=(a/4)*[-1 -1 1];
T3=(a/4)*[1 -1 -1];
T4=(a/4)*[-1 1 -1];

x1=dot(k,T1);
x2=dot(k,T2);
x3=dot(k,T3);
x4=dot(k,T4);

g0=(exp(1i*x1)+exp(1i*x2)+exp(1i*x3)+exp(1i*x4))/4;
g1=(exp(1i*x1)-exp(1i*x2)+exp(1i*x3)-exp(1i*x4))/4;
g2=(exp(1i*x1)-exp(1i*x2)-exp(1i*x3)+exp(1i*x4))/4;
g3=(exp(1i*x1)+exp(1i*x2)-exp(1i*x3)-exp(1i*x4))/4;

H=[Esa           0               0               0              0              g0*Esasc      g1*Esaxc        g2*Esaxc        g3*Esaxc          0
   0             Epa             -1i*(deltaa/3)  (deltaa/3)     0              -g1*Exasc     g0*Exaxc        g3*Exayc        g2*Exayc          -g1*Exassc
   0             1i*(deltaa/3)   Epa             -1i*(deltaa/3) 0              -g2*Exasc     g3*Exayc        g0*Exaxc        g1*Exayc          -g2*Exassc
   0             (deltaa/3)      1i*(deltaa/3)   Epa            0              -g3*Exasc     g2*Exayc        g1*Exayc        g0*Exaxc          -g3*Exassc
   0             0               0               0              Essa           0             g1*Essaxc       g2*Essaxc       g3*Essaxc         0
   g0'*Esasc     -g1'*Exasc      -g2'*Exasc      -g3'*Exasc     0              Esc           0               0               0                 0
   g1'*Esaxc     g0'*Exaxc       g3'*Exayc       g2'*Exayc      g1'*Essaxc      0             Epc            -1i*(deltac/3)  (deltac/3)        0
   g2'*Esaxc     g3'*Exayc       g0'*Exaxc       g1'*Exayc      g2'*Essaxc      0             1i*(deltac/3)  Epc             -1i*(deltac/3)    0 
   g3'*Esaxc     g2'*Exayc       g1'*Exayc       g0'*Exaxc      g3'*Essaxc      0             (deltac/3)     1i*(deltac/3)   Epc               0
   0             -g1'*Exassc     -g2'*Exassc     -g3'*Exassc    0               0             0              0               0                 Essc];
 
energy_level=eig(H);



% Eg =abs(energy_level(5,1)-energy_level(4,1));
end