function [Hac]= Hac_gen(f1,f2,beta,cosx,cosy,cosz,Esasc,Esaxc,Exasc,Exaxc,Exayc,Exassc,Essaxc)


 hac=[ beta*f1*Esasc                      beta*sqrt(3)*f2*Esaxc*cosx             beta*sqrt(3)*f2*Esaxc*cosy                     beta*sqrt(3)*f1*Esaxc*cosz                 0
       -beta*sqrt(3)*f2*Exasc*cosx        beta*f1*(Exaxc+(3*(cosx^2)-1)*Exayc)   beta*3*f1*Exayc*cosy*cosx                      beta*3*f2*Exayc*cosz*cosx                  -sqrt(3)*beta*f2*Exassc*cosx
       -beta*sqrt(3)*f2*Exasc*cosy        beta*3*f1*Exayc*cosy*cosx              beta*f1*(Exaxc+((3*(cosy^2))-1)*Exayc)         beta*3*f2*Exayc*cosz*cosy                  -sqrt(3)*beta*f2*Exassc*cosy
       -beta*sqrt(3)*f1*Exasc*cosz        beta*3*f2*Exayc*cosz*cosx              beta*3*f2*Exayc*cosz*cosy                      beta*f1*(Exaxc+((3*(cosz^2))-1)*Exayc)     -sqrt(3)*beta*f1*Exassc*cosz
       0                                  sqrt(3)*beta*f2*Essaxc*cosx            sqrt(3)*beta*f2*Essaxc*cosy                    sqrt(3)*beta*f1*Essaxc*cosz                0   ];
   
   
 zz=zeros(5);
  Hac=[hac    zz
       zz      hac];
end
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   