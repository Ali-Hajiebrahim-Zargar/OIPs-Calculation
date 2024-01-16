function [Hca]= Hca_gen(g1,g2,beta,cosx,cosy,cosz,Esasc,Esaxc,Exasc,Exaxc,Exayc,Exassc,Essaxc)


 hca=[ beta*g1*Esasc                      -beta*sqrt(3)*g2*Exasc*cosx             beta*sqrt(3)*g2*Exasc*cosy                     beta*sqrt(3)*g1*Exasc*cosz                 0
       beta*sqrt(3)*g2*Esaxc*cosx         beta*g1*(Exaxc+(3*(cosx^2)-1)*Exayc)    -beta*3*g1*Exayc*cosy*cosx                     -beta*3*g2*Exayc*cosz*cosx                 sqrt(3)*beta*g2*Essaxc*cosx
       -beta*sqrt(3)*g2*Esaxc*cosy        -beta*3*g1*Exayc*cosy*cosx              beta*g1*(Exaxc+((3*(cosy^2))-1)*Exayc)         beta*3*g2*Exayc*cosz*cosy                  -sqrt(3)*beta*g2*Essaxc*cosy
       -beta*sqrt(3)*g1*Esaxc*cosz        -beta*3*g2*Exayc*cosz*cosx              beta*3*g2*Exayc*cosz*cosy                      beta*g1*(Exaxc+((3*(cosz^2))-1)*Exayc)     -sqrt(3)*beta*g1*Essaxc*cosz
       0                                  -sqrt(3)*beta*g2*Exassc*cosx            sqrt(3)*beta*g2*Exassc*cosy                    sqrt(3)*beta*g1*Exassc*cosz                0   ];
   
   
 zz=zeros(5);
  Hca=[hca    zz
       zz      hca];
end
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   