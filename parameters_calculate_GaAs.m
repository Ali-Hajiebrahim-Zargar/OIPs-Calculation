
function [FitnessFunction,Epa,Epc,Exaxc,Esasc]=parameters_calculate_GaAs(Esa,Esc,Essa,Essc,Exayc,Esaxc,Exasc,Essaxc,Exassc)



% m_star_Gamma_tensor_massE100   m_star_Gamma_tensor_massHH100    m_star_Gamma_tensor_massLH100   m_star_Gamma_tensor_massSO100     m_star_Gamma_tensor_massE110    m_star_Gamma_tensor_massHH110   m_star_Gamma_tensor_massLH110      m_star_Gamma_tensor_massSO110        m_star_Gamma_tensor_massE111   m_star_Gamma_tensor_massHH111     m_star_Gamma_tensor_massLH111_out     m_star_Gamma_tensor_massSO111        Eg_gamma_out      Eg_gamma_6_valence    Eg_gamma_7_conduction   Eg_gamma_8_conduction         Ehh_gamma    Elh_gamma    ESO_gamma    m_star_X_transvers   m_star_X_longitude   Eg_X_7_conduction    Eg_X_6_conduction    Eg_X_7_valence   Eg_X_6_valence     Eg_X_5_valence   X_level    minimomOFenergy_X    Eg_L       m_star_L_transvers        m_star_L_longitude   Eg_L_7_conduction  Eg_L_7_valence   Eg_L_6_valence   Eg_L_5_valence_out
 Fin=[100000                               10000                                  100000                       10000                              100000                    10000                             100000                                 10000                            1                           10000                              1                                    1                                1000                  10                    10                     10                         10000      10000          1000                1                     1                     10              100                  10                10                100        100              100            100                1                       1                   1                  1                 1                    100          ];

% Fin=ones(1,35);



% err=zeros(35,1);
ali=0;
s=1000000;
h=1.05458e-34;
minimomOFenergy_X=20;
X_level=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%fitting parameters %%%%%%%%%%%%%%%%%%%%%%
%at GAMA point
m_star_Gamma_tensor_massE100_define=0.067;
m_star_Gamma_tensor_massHH100_define=-0.403;
m_star_Gamma_tensor_massLH100_define=-0.087;
m_star_Gamma_tensor_massSO100_define=-0.150;

m_star_Gamma_tensor_massE110_define=0.067;
m_star_Gamma_tensor_massHH110_define=-0.660;
m_star_Gamma_tensor_massLH110_define=-0.080;
m_star_Gamma_tensor_massSO110_define=-0.150;

m_star_Gamma_tensor_massE111_define=0.067;
m_star_Gamma_tensor_massHH111_define=-0.813;
m_star_Gamma_tensor_massLH111_define=-0.079;
m_star_Gamma_tensor_massSO111_define=-0.150;

Eg_gamma_define=1.424;
Ehh_gamma_define=0.0000;
Elh_gamma_define=0.0000;
ESO_gamma_define=-0.340;
Eg_gamma_6_valence_define=-13.100;
Eg_gamma_7_conduction_define=4.530;
Eg_gamma_8_conduction_define=4.716;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%at X point
X_level_define=0.90;
minimomOFenergy_X_define=1.900;
m_star_X_transvers_define=0.230;
m_star_X_longitude_define=1.300;
Eg_X_7_conduction_define=2.320;
Eg_X_6_conduction_define=1.980;
Eg_X_7_valence_define=-2.800;
Eg_X_6_valence_define=-2.880;
Eg_X_5_valence_define=-6.800;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%at L point
m_star_L_tensor_transvers_define=0.075;
m_star_L_longitude_define=1.900;
Eg_L_define=1.708;
Eg_L_7_valence_define=-1.200;
Eg_L_6_valence_define=-1.42;
Eg_L_5_valence_define=-8.000;
Eg_L_7_conduction_define=5.470;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a=5.6660;
deltaa=0.421;
deltac=0.174;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Yajon parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta=-ESO_gamma_define;
Epa=(Eg_gamma_7_conduction_define*(deltaa-delta)-deltaa*delta+(2/3)*deltaa^2+deltaa*delta/3)/(deltaa-deltac);
Epc=(delta^2+delta*(Epa-(2*deltaa)/3-(2*deltac)/3)+deltaa*deltac/3-Epa*deltac)/(deltaa-delta);
Exaxc=sqrt((Epa+(deltaa/3))*(Epc+(deltac/3)));
Esasc=-sqrt(Eg_gamma_define^2-Eg_gamma_define*(Esa+Esc)+Esa*Esc);


%%%%%%%%%%%%%%%%%%%%%%%%% Gamma point %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%effective masses
%%%%%100
kx=0;
ky=0;
kz=0;
[energy_level_Gamma1,~]=bulk(kx,ky,kz,a,Esa,Esc,Epa,Epc,Essa,Essc,Esasc,Esaxc,Exasc,Essaxc,Exassc,Exaxc,Exayc,deltaa,deltac);
%
kx=0-(1/s);
ky=0;
kz=0;
[energy_level_Gamma2,~]=bulk(kx,ky,kz,a,Esa,Esc,Epa,Epc,Essa,Essc,Esasc,Esaxc,Exasc,Essaxc,Exassc,Exaxc,Exayc,deltaa,deltac);
%
kx=0+(1/s);
ky=0;
kz=0;
[energy_level_Gamma3,~]=bulk(kx,ky,kz,a,Esa,Esc,Epa,Epc,Essa,Essc,Esasc,Esaxc,Exasc,Essaxc,Exassc,Exaxc,Exayc,deltaa,deltac);

%%%%%%%%%%%%%%%%%%%%% effective mass of electron at Gamma point
m_star_Gamma_tensor1=((energy_level_Gamma3(5,1)-2*energy_level_Gamma1(5,1)+energy_level_Gamma2(5,1))/(((1/s))^2))*1.60218e-39;
m_star_Gamma_tensor_massE100=(h^2)/(m_star_Gamma_tensor1*(0.91095e-30));
%%%%%%%%%%%%%%%%%%%%% effective mass of heavy hole at Gamma point
m_star_Gamma_tensor2=((energy_level_Gamma3(4,1)-2*energy_level_Gamma1(4,1)+energy_level_Gamma2(4,1))/(((1/s))^2))*1.60218e-39;
m_star_Gamma_tensor_massHH100=(h^2)/(m_star_Gamma_tensor2*(0.91095e-30));
%%%%%%%%%%%%%%%%%%%%% effective mass of light hole at Gamma point
m_star_Gamma_tensor3=((energy_level_Gamma3(3,1)-2*energy_level_Gamma1(3,1)+energy_level_Gamma2(3,1))/(((1/s))^2))*1.60218e-39;
m_star_Gamma_tensor_massLH100=(h^2)/(m_star_Gamma_tensor3*(0.91095e-30));
%%%%%%%%%%%%%%%%%%%%% effective mass of spin orbit at Gamma point
m_star_Gamma_tensor4=((energy_level_Gamma3(2,1)-2*energy_level_Gamma1(2,1)+energy_level_Gamma2(2,1))/(((1/s))^2))*1.60218e-39;
m_star_Gamma_tensor_massSO100=(h^2)/(m_star_Gamma_tensor4*(0.91095e-30));
%
%
%
%
%
%%%%%110
kx=0;
ky=0-(1/s);
kz=0-(1/s);
[energy_level_Gamma4,~]=bulk(kx,ky,kz,a,Esa,Esc,Epa,Epc,Essa,Essc,Esasc,Esaxc,Exasc,Essaxc,Exassc,Exaxc,Exayc,deltaa,deltac);
%
kx=0;
ky=0+(1/s);
kz=0+(1/s);
[energy_level_Gamma5,~]=bulk(kx,ky,kz,a,Esa,Esc,Epa,Epc,Essa,Essc,Esasc,Esaxc,Exasc,Essaxc,Exassc,Exaxc,Exayc,deltaa,deltac);
%
%%%%%%%%%%%%%%%%%%%%% effective mass of electron at Gamma point
m_star_Gamma_tensor5=((energy_level_Gamma4(5,1)-2*energy_level_Gamma1(5,1)+energy_level_Gamma5(5,1))/(2*((1/s))^2))*1.60218e-39;
m_star_Gamma_tensor_massE110=(h^2)/(m_star_Gamma_tensor5*(0.91095e-30));
%%%%%%%%%%%%%%%%%%%%% effective mass of heavy hole at Gamma point
m_star_Gamma_tensor5=((energy_level_Gamma4(4,1)-2*energy_level_Gamma1(4,1)+energy_level_Gamma5(4,1))/(2*((1/s))^2))*1.60218e-39;
m_star_Gamma_tensor_massHH110=(h^2)/(m_star_Gamma_tensor5*(0.91095e-30));
%%%%%%%%%%%%%%%%%%%%% effective mass of light hole at Gamma point
m_star_Gamma_tensor5=((energy_level_Gamma4(3,1)-2*energy_level_Gamma1(3,1)+energy_level_Gamma5(3,1))/(2*((1/s))^2))*1.60218e-39;
m_star_Gamma_tensor_massLH110=(h^2)/(m_star_Gamma_tensor5*(0.91095e-30));
%%%%%%%%%%%%%%%%%%%%% effective mass of spin orbit at Gamma point
m_star_Gamma_tensor6=((energy_level_Gamma4(2,1)-2*energy_level_Gamma1(2,1)+energy_level_Gamma5(2,1))/(2*((1/s))^2))*1.60218e-39;
m_star_Gamma_tensor_massSO110=(h^2)/(m_star_Gamma_tensor6*(0.91095e-30));
%
%
%
%
%
%%%%%111
kx=0-(1/s);
ky=0-((1/s));
kz=0-((1/s));
[energy_level_Gamma6,~]=bulk(kx,ky,kz,a,Esa,Esc,Epa,Epc,Essa,Essc,Esasc,Esaxc,Exasc,Essaxc,Exassc,Exaxc,Exayc,deltaa,deltac);
%
kx=0+(1/s);
ky=0+((1/s));
kz=0+((1/s));
[energy_level_Gamma7,~]=bulk(kx,ky,kz,a,Esa,Esc,Epa,Epc,Essa,Essc,Esasc,Esaxc,Exasc,Essaxc,Exassc,Exaxc,Exayc,deltaa,deltac);
%
%%%%%%%%%%%%%%%%%%%%% effective mass of electron at Gamma point
m_star_Gamma_tensor7=((energy_level_Gamma6(5,1)-2*energy_level_Gamma1(5,1)+energy_level_Gamma7(5,1))/(3*((1/s))^2))*1.60218e-39;
m_star_Gamma_tensor_massE111=(h^2)/(m_star_Gamma_tensor7*(0.91095e-30));
%%%%%%%%%%%%%%%%%%%%% effective mass of heavy hole at Gamma point
m_star_Gamma_tensor8=((energy_level_Gamma6(4,1)-2*energy_level_Gamma1(4,1)+energy_level_Gamma7(4,1))/(3*((1/s))^2))*1.60218e-39;
m_star_Gamma_tensor_massHH111=(h^2)/(m_star_Gamma_tensor8*(0.91095e-30));
%%%%%%%%%%%%%%%%%%%%% effective mass of light hole at Gamma point
m_star_Gamma_tensor9=((energy_level_Gamma6(3,1)-2*energy_level_Gamma1(3,1)+energy_level_Gamma7(3,1))/(3*((1/s))^2))*1.60218e-39;
m_star_Gamma_tensor_massLH111=(h^2)/(m_star_Gamma_tensor9*(0.91095e-30));
%%%%%%%%%%%%%%%%%%%%% effective mass of spin orbit at Gamma point
m_star_Gamma_tensor10=((energy_level_Gamma6(2,1)-2*energy_level_Gamma1(2,1)+energy_level_Gamma7(2,1))/(3*((1/s))^2))*1.60218e-39;
m_star_Gamma_tensor_massSO111=(h^2)/(m_star_Gamma_tensor10*(0.91095e-30));
%
%%%%%%%%Energy levels at gamma point
Eg_gamma=energy_level_Gamma1(5,1);
Ehh_gamma=energy_level_Gamma1(4,1);
Elh_gamma=energy_level_Gamma1(3,1);
ESO_gamma=energy_level_Gamma1(2,1);
Eg_gamma_6_valence=energy_level_Gamma1(1,1);
Eg_gamma_7_conduction=energy_level_Gamma1(6,1);
Eg_gamma_8_conduction=energy_level_Gamma1(7,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% X point %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find minima in X path
for i=0.5:0.001:1
kx=(i*2*pi/a);
ky=0;
kz=0;
[energy_level_GammaToX,~]=bulk(kx,ky,kz,a,Esa,Esc,Epa,Epc,Essa,Essc,Esasc,Esaxc,Exasc,Essaxc,Exassc,Exaxc,Exayc,deltaa,deltac);
if minimomOFenergy_X > energy_level_GammaToX(5,1)
    minimomOFenergy_X=energy_level_GammaToX(5,1);
    X_level=i;
end
end
%
%%%%%%%%%%%%%%%%%%%%%% for X point in Transvers effective mass
kx=(X_level*2*pi/a);
ky=0;
kz=0;
[energy_level_X1,~]=bulk(kx,ky,kz,a,Esa,Esc,Epa,Epc,Essa,Essc,Esasc,Esaxc,Exasc,Essaxc,Exassc,Exaxc,Exayc,deltaa,deltac);
%
kx=(X_level*2*pi/a);
ky=0-(1/s);
kz=0;
[energy_level_X2,~]=bulk(kx,ky,kz,a,Esa,Esc,Epa,Epc,Essa,Essc,Esasc,Esaxc,Exasc,Essaxc,Exassc,Exaxc,Exayc,deltaa,deltac);
%
kx=(X_level*2*pi/a);
ky=0+(1/s);
kz=0;
[energy_level_X3,~]=bulk(kx,ky,kz,a,Esa,Esc,Epa,Epc,Essa,Essc,Esasc,Esaxc,Exasc,Essaxc,Exassc,Exaxc,Exayc,deltaa,deltac);
%

m_star_X_tensor_transvers=((energy_level_X3(5,1)-2*energy_level_X1(5,1)+energy_level_X2(5,1))/((1/s)^2))*1.60218e-39;


%%%%%%%%%%%%%%%%%%%%%% for X point in longitude effective mass

kx=(X_level*2*pi/a)-(1/s);
ky=0;
kz=0;
[energy_level_X4,~]=bulk(kx,ky,kz,a,Esa,Esc,Epa,Epc,Essa,Essc,Esasc,Esaxc,Exasc,Essaxc,Exassc,Exaxc,Exayc,deltaa,deltac);
%
kx=(X_level*2*pi/a)+(1/s);
ky=0;
kz=0;
[energy_level_X5,~]=bulk(kx,ky,kz,a,Esa,Esc,Epa,Epc,Essa,Essc,Esasc,Esaxc,Exasc,Essaxc,Exassc,Exaxc,Exayc,deltaa,deltac);
%

m_star_X_tensor_longitude=((energy_level_X4(5,1)-2*energy_level_X1(5,1)+energy_level_X5(5,1))/((1/s)^2))*1.60218e-39;

m_star_X_transvers=(h^2)/(m_star_X_tensor_transvers*(0.91095e-30));

m_star_X_longitude=(h^2)/(m_star_X_tensor_longitude*(0.91095e-30));

kx=2*pi/a;
ky=0;
kz=0;
[energy_level_X6,~]=bulk(kx,ky,kz,a,Esa,Esc,Epa,Epc,Essa,Essc,Esasc,Esaxc,Exasc,Essaxc,Exassc,Exaxc,Exayc,deltaa,deltac);

Eg_X_7_conduction=energy_level_X6(6,1);
Eg_X_6_conduction=energy_level_X6(5,1);
Eg_X_7_valence=energy_level_X6(4,1);
Eg_X_6_valence=energy_level_X6(3,1);
Eg_X_5_valence=energy_level_X6(2,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% L point %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%transvers
kx=pi/a;
ky=pi/a;
kz=pi/a;
[energy_level_L1,~]=bulk(kx,ky,kz,a,Esa,Esc,Epa,Epc,Essa,Essc,Esasc,Esaxc,Exasc,Essaxc,Exassc,Exaxc,Exayc,deltaa,deltac);
%
kx=(pi/a)-(1/s);
ky=(pi/a);
kz=(pi/a)+(1/s);
[energy_level_L2,~]=bulk(kx,ky,kz,a,Esa,Esc,Epa,Epc,Essa,Essc,Esasc,Esaxc,Exasc,Essaxc,Exassc,Exaxc,Exayc,deltaa,deltac);
%
kx=(pi/a)+(1/s);
ky=(pi/a);
kz=(pi/a)-(1/s);
[energy_level_L3,~]=bulk(kx,ky,kz,a,Esa,Esc,Epa,Epc,Essa,Essc,Esasc,Esaxc,Exasc,Essaxc,Exassc,Exaxc,Exayc,deltaa,deltac);
%
m_star_L_tensor1=((energy_level_L3(5,1)-2*energy_level_L1(5,1)+energy_level_L2(5,1))/(2*((1/s)^2)))*1.60218e-39;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%longitude
kx=(pi/a)-(1/s);
ky=(pi/a)-(1/s);
kz=(pi/a)-(1/s);
[energy_level_L4,~]=bulk(kx,ky,kz,a,Esa,Esc,Epa,Epc,Essa,Essc,Esasc,Esaxc,Exasc,Essaxc,Exassc,Exaxc,Exayc,deltaa,deltac);
%
kx=(pi/a)-(2/s);
ky=(pi/a)-(2/s);
kz=(pi/a)-(2/s);
[energy_level_L5,~]=bulk(kx,ky,kz,a,Esa,Esc,Epa,Epc,Essa,Essc,Esasc,Esaxc,Exasc,Essaxc,Exassc,Exaxc,Exayc,deltaa,deltac);
%

m_star_L_tensor2=((energy_level_L1(5,1)-2*energy_level_L4(5,1)+energy_level_L5(5,1))/(3*((1/s)^2)))*1.60218e-39;

m_star_L_tensor_transvers=(h^2)/(m_star_L_tensor1*(0.91095e-30));
m_star_L_longitude=(h^2)/(m_star_L_tensor2*(0.91095e-30));

Eg_L_7_conduction=energy_level_L1(6,1);
Eg_L=energy_level_L1(5,1);
Eg_L_7_valence=energy_level_L1(4,1);
Eg_L_6_valence=energy_level_L1(3,1);
Eg_L_5_valence=energy_level_L1(2,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fitting to zero %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m_star_Gamma_tensor_massE100_out=abs(m_star_Gamma_tensor_massE100-m_star_Gamma_tensor_massE100_define);
m_star_Gamma_tensor_massHH100_out=abs(m_star_Gamma_tensor_massHH100-m_star_Gamma_tensor_massHH100_define);
m_star_Gamma_tensor_massLH100_out=abs(m_star_Gamma_tensor_massLH100-m_star_Gamma_tensor_massLH100_define);
m_star_Gamma_tensor_massSO100_out=abs(m_star_Gamma_tensor_massSO100-m_star_Gamma_tensor_massSO100_define);

m_star_Gamma_tensor_massE110_out=abs(m_star_Gamma_tensor_massE110-m_star_Gamma_tensor_massE110_define);
m_star_Gamma_tensor_massHH110_out=abs(m_star_Gamma_tensor_massHH110-m_star_Gamma_tensor_massHH110_define);
m_star_Gamma_tensor_massLH110_out=abs(m_star_Gamma_tensor_massLH110-m_star_Gamma_tensor_massLH110_define);
m_star_Gamma_tensor_massSO110_out=abs(m_star_Gamma_tensor_massSO110-m_star_Gamma_tensor_massSO110_define);

m_star_Gamma_tensor_massE111_out=abs(m_star_Gamma_tensor_massE111-m_star_Gamma_tensor_massE111_define);
m_star_Gamma_tensor_massHH111_out=abs(m_star_Gamma_tensor_massHH111-m_star_Gamma_tensor_massHH111_define);
m_star_Gamma_tensor_massLH111_out=abs(m_star_Gamma_tensor_massLH111-m_star_Gamma_tensor_massLH111_define);
m_star_Gamma_tensor_massSO111_out=abs(m_star_Gamma_tensor_massSO111-m_star_Gamma_tensor_massSO111_define);

Eg_gamma_out=abs(Eg_gamma-Eg_gamma_define);
Ehh_gamma_out=abs(Ehh_gamma-Ehh_gamma_define);
Elh_gamma_out=abs(Elh_gamma-Elh_gamma_define);
ESO_gamma_out=abs(ESO_gamma-ESO_gamma_define);
Eg_gamma_6_valence_out=abs(Eg_gamma_6_valence-Eg_gamma_6_valence_define);
Eg_gamma_7_conduction_out=abs(Eg_gamma_7_conduction-Eg_gamma_7_conduction_define);
Eg_gamma_8_conduction_out=abs(Eg_gamma_8_conduction-Eg_gamma_8_conduction_define);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m_star_X_transvers_out=abs(m_star_X_transvers-m_star_X_transvers_define);
m_star_X_longitude_out=abs(m_star_X_longitude-m_star_X_longitude_define);

Eg_X_7_conduction_out=abs(Eg_X_7_conduction-Eg_X_7_conduction_define);
Eg_X_6_conduction_out=abs(Eg_X_6_conduction-Eg_X_6_conduction_define);
Eg_X_7_valence_out=abs(Eg_X_7_valence-Eg_X_7_valence_define);
Eg_X_6_valence_out=abs(Eg_X_6_valence-Eg_X_6_valence_define);
Eg_X_5_valence_out=abs(Eg_X_5_valence-Eg_X_5_valence_define);


X_level_out=abs(X_level-X_level_define);
minimomOFenergy_X_out=abs(minimomOFenergy_X-minimomOFenergy_X_define);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m_star_L_tensor_transvers_out=abs(m_star_L_tensor_transvers-m_star_L_tensor_transvers_define);
m_star_L_longitude_out=abs(m_star_L_longitude-m_star_L_longitude_define);

Eg_L_7_conduction_out=abs(Eg_L_7_conduction-Eg_L_7_conduction_define);
Eg_L_out=abs(Eg_L-Eg_L_define);
Eg_L_7_valence_out=abs(Eg_L_7_valence-Eg_L_7_valence_define);
Eg_L_6_valence_out=abs(Eg_L_6_valence-Eg_L_6_valence_define);
Eg_L_5_valence_out=abs(Eg_L_5_valence-Eg_L_5_valence_define);


% %sign consideration
% 1
if ~((m_star_Gamma_tensor_massE100>0.060)&&(0.070>m_star_Gamma_tensor_massE100))
    
    ali=100000+ali;
%     err(1,1)=1;
end
%2
if ~((m_star_Gamma_tensor_massHH100>-0.5)&&(-0.3>m_star_Gamma_tensor_massHH100))
    
    ali=100000+ali;
%     err(2,1)=1;
end
%3
if ~((m_star_Gamma_tensor_massLH100>-0.095)&&(-0.075>m_star_Gamma_tensor_massLH100))
    
    ali=100000+ali;
%     err(3,1)=1;
end
%4
if ~((m_star_Gamma_tensor_massSO100>-0.165)&&(-0.145>m_star_Gamma_tensor_massSO100))
    
    ali=100000+ali;
%     err(4,1)=1;
end
%5
if ~((m_star_Gamma_tensor_massE110>0.060)&&(0.070>m_star_Gamma_tensor_massE110))
%     err(5,1)=1;
    ali=100000+ali;
end
%6
if ~((m_star_Gamma_tensor_massHH110>-0.670)&&(-0.650>m_star_Gamma_tensor_massHH110))
%     err(6,1)=1;
    ali=100000+ali;
end
%7
if ~((m_star_Gamma_tensor_massLH110>-0.090)&&(-0.070>m_star_Gamma_tensor_massLH110))
%     err(7,1)=1;
    ali=100000+ali;
end
%8
if ~((m_star_Gamma_tensor_massSO110>-0.169)&&(-0.140>m_star_Gamma_tensor_massSO110))
%     err(8,1)=1;
    ali=100000+ali;
end
%9
% if ~((m_star_Gamma_tensor_massE111>0.060)&&(0.073>m_star_Gamma_tensor_massE111))
%     err(9,1)=1;
%     ali=100000+ali;
% end
%10
if ~((m_star_Gamma_tensor_massHH111>-0.9)&&(-0.75>m_star_Gamma_tensor_massHH111))
%     err(10,1)=1;
    ali=100000+ali;
end
%11
% if ~((m_star_Gamma_tensor_massLH111>-0.09)&&(-0.070>m_star_Gamma_tensor_massLH111))
%     err(11,1)=1;
%     ali=100000+ali;
% end
%12
% if ~((m_star_Gamma_tensor_massSO111>-0.16)&&(-0.14>m_star_Gamma_tensor_massSO111))
%     err(12,1)=1;
%     ali=100000+ali;
% end
% %13
if ~((Eg_gamma>1.4240)&&(1.4249>Eg_gamma))
%     err(13,1)=1;
    ali=100000+ali;
end
%14
if ~((Eg_gamma_6_valence>-14)&&(-12>Eg_gamma_6_valence))
%     err(14,1)=1;
    ali=100000+ali;
end
%15
if ~((Eg_gamma_7_conduction>3.5)&&(5.5>Eg_gamma_7_conduction))
%     err(15,1)=1;
    ali=100000+ali;
end
%16
if ~((Eg_gamma_8_conduction>3.7)&&(5.7>Eg_gamma_8_conduction))
%     err(16,1)=1;
    ali=100000+ali;
end
%17
if ~((0.0001>Ehh_gamma))
%     err(17,1)=1;
    ali=100000+ali;
end
%18
if ~((0.0001>Elh_gamma))
%     err(18,1)=1;
    ali=100000+ali;
end
% %19
if ~((ESO_gamma<-0.30)&&(-0.36<ESO_gamma))
%     err(19,1)=1;
    ali=100000+ali;
end
%20
if ~((m_star_X_transvers>0)&&(5>m_star_X_transvers))
%     err(20,1)=1;
    ali=100000+ali;
end
%21
if ~((m_star_X_longitude>1)&&(2>m_star_X_longitude))
%     err(21,1)=1;
    ali=100000+ali;
end
%22
if ~((Eg_X_7_conduction>1.5)&&(3.5>Eg_X_7_conduction))
%     err(22,1)=1;
    ali=100000+ali;
end
%23
if ~((Eg_X_6_conduction>1)&&(3>Eg_X_6_conduction))
%     err(23,1)=1;
    ali=100000+ali;
end
%24
if ~((Eg_X_7_valence>-4)&&(-1>Eg_X_7_valence))
%     err(24,1)=1;
    ali=100000+ali;
end
%25
if ~((Eg_X_6_valence>-4)&&(-0.5>Eg_X_6_valence))
%     err(25,1)=1;
    ali=100000+ali;
end
%26
if ~((Eg_X_5_valence>-8)&&(-4>Eg_X_5_valence))
%     err(28,1)=1;
    ali=100000+ali;
end
%27
if ~((X_level>0.8)&&(1>X_level))
%     err(27,1)=1;
    ali=100000+ali;
end
%28
if ~((minimomOFenergy_X>1.7)&&(2>minimomOFenergy_X))
%     err(28,1)=1;
    ali=100000+ali;
end
%29
if ~((Eg_L>1)&&(3>Eg_L))
%     err(29,1)=1;
    ali=100000+ali;
end
%30
if ~((m_star_L_tensor_transvers>0)&&(1>m_star_L_tensor_transvers))
%     err(30,1)=1;
    ali=100000+ali;
end
%31
if ~((m_star_L_longitude>1)&&(2.2>m_star_L_longitude))
%     err(31,1)=1;
    ali=100000+ali;
end
%32
if ~((Eg_L_7_conduction>3)&&(8>Eg_L_7_conduction))
%     err(32,1)=1;
    ali=100000+ali;
end
%33
if ~((Eg_L_7_valence>-3)&&(-0.5>Eg_L_7_valence))
%     err(33,1)=1;
    ali=100000+ali;
end
%34
if ~((Eg_L_6_valence>-3)&&(-0.5>Eg_L_6_valence))
%     err(34,1)=1;
    ali=100000+ali;
end
%35
if ~((Eg_L_5_valence>-11)&&(-4>Eg_L_5_valence))
%     err(35,1)=1;
    ali=100000+ali;
end


  FitnessFunction=(Fin(1,1)*m_star_Gamma_tensor_massE100_out+Fin(1,2)*m_star_Gamma_tensor_massHH100_out+Fin(1,3)*m_star_Gamma_tensor_massLH100_out+Fin(1,4)*m_star_Gamma_tensor_massSO100_out+Fin(1,5)*m_star_Gamma_tensor_massE110_out+Fin(1,6)*m_star_Gamma_tensor_massHH110_out+Fin(1,7)*m_star_Gamma_tensor_massLH110_out+Fin(1,8)*m_star_Gamma_tensor_massSO110_out+Fin(1,9)*m_star_Gamma_tensor_massE111_out+Fin(1,10)*m_star_Gamma_tensor_massHH111_out+Fin(1,11)*m_star_Gamma_tensor_massLH111_out+Fin(1,12)*m_star_Gamma_tensor_massSO111_out+Fin(1,13)*Eg_gamma_out+Fin(1,14)*Eg_gamma_6_valence_out+Fin(1,15)*Eg_gamma_7_conduction_out+Fin(1,16)*Eg_gamma_8_conduction_out+Fin(1,17)*Ehh_gamma_out+Fin(1,18)*Elh_gamma_out+Fin(1,19)*ESO_gamma_out+Fin(1,20)*m_star_X_transvers_out+Fin(1,21)*m_star_X_longitude_out+Fin(1,22)*Eg_X_7_conduction_out+Fin(1,23)*Eg_X_6_conduction_out+Fin(1,24)*Eg_X_7_valence_out+Fin(1,25)*Eg_X_6_valence_out+Fin(1,26)*Eg_X_5_valence_out+Fin(1,27)*X_level_out+Fin(1,28)*minimomOFenergy_X_out+Fin(1,29)*Eg_L_out+Fin(1,30)*m_star_L_tensor_transvers_out+Fin(1,31)*m_star_L_longitude_out+Fin(1,32)*Eg_L_7_conduction_out+Fin(1,33)*Eg_L_7_valence_out+Fin(1,34)*Eg_L_6_valence_out+Fin(1,35)*Eg_L_5_valence_out+ali);

%FitnessFunction=(1/35)*(m_star_Gamma_tensor_massE100_out+m_star_Gamma_tensor_massHH100_out+m_star_Gamma_tensor_massLH100_out+m_star_Gamma_tensor_massSO100_out+m_star_Gamma_tensor_massE110_out+m_star_Gamma_tensor_massHH110_out+m_star_Gamma_tensor_massLH110_out+m_star_Gamma_tensor_massSO110_out+m_star_Gamma_tensor_massE111_out+m_star_Gamma_tensor_massHH111_out+m_star_Gamma_tensor_massLH111_out+m_star_Gamma_tensor_massSO111_out+Eg_gamma_out+Eg_gamma_6_valence_out+Eg_gamma_7_conduction_out+Eg_gamma_8_conduction_out+Ehh_gamma_out+Elh_gamma_out+ESO_gamma_out+m_star_X_transvers_out+m_star_X_longitude_out+Eg_X_7_conduction_out+Eg_X_6_conduction_out+Eg_X_7_valence_out+Eg_X_6_valence_out+Eg_X_5_valence_out+X_level_out+minimomOFenergy_X_out+Eg_L_out+m_star_L_tensor_transvers_out+m_star_L_longitude_out+Eg_L_7_conduction_out+Eg_L_7_valence_out+Eg_L_6_valence_out+Eg_L_5_valence_out+ali);



end





