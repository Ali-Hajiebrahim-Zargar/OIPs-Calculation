function [Hsl,E]=SL_Hamiltonian(ct,inputss,atoms,made1chiye,m,kx,ky,kz,bandlineup)
k=[kx ky kz];
% Hsl=zeros(atoms*10);
%%%%%%%%%%%%%%%%%% mater1 hamiltonian construction %%%%%%%%%%%%%%%%%%%%%%%%
[satr1,satr2,satr3,satr4]= satr_gen(inputss(2,1));
x=inputss(2,2);
y=inputss(2,3);
[a,D,Esa,Esc,Epa,Epc,Essa,Essc,Esasc,Esaxc,Exasc,Essaxc,Exassc,Exaxc,Exayc,Deltaa,Deltac,bandlineup_new]=energy_generator(x,y,inputss,bandlineup,satr1,satr2,satr3,satr4,ct);
[f1,f2,g1,g2,cosx,cosy,cosz,beta]=tao_gen(a,D,inputss,ct,k);
[Ha1]=Ha_gen(Esa,Epa,Essa,Deltaa,bandlineup_new);
[Hc1]=Hc_gen(Esc,Epc,Essc,Deltac,bandlineup_new);
[Hac1]= Hac_gen(f1,f2,beta,cosx,cosy,cosz,Esasc,Esaxc,Exasc,Exaxc,Exayc,Exassc,Essaxc);
[Hca1]= Hca_gen(g1,g2,beta,cosx,cosy,cosz,Esasc,Esaxc,Exasc,Exaxc,Exayc,Exassc,Essaxc);
%%%%%%%%%%%%%%%%%% interface1 hamiltonian construction %%%%%%%%%%%%%%%%%%%%
[satr1,satr2,satr3,satr4]= satr_gen(inputss(3,1));
x=inputss(3,2);
y=inputss(3,3);
[a,D,Esa,Esc,Epa,Epc,Essa,Essc,Esasc,Esaxc,Exasc,Essaxc,Exassc,Exaxc,Exayc,Deltaa,Deltac,bandlineup_new]=energy_generator(x,y,inputss,bandlineup,satr1,satr2,satr3,satr4,ct);
[f1,f2,g1,g2,cosx,cosy,cosz,beta]=tao_gen(a,D,inputss,ct,k);
[Ha2]=Ha_gen(Esa,Epa,Essa,Deltaa,bandlineup_new);
[Hc2]=Hc_gen(Esc,Epc,Essc,Deltac,bandlineup_new);
[Hac2]= Hac_gen(f1,f2,beta,cosx,cosy,cosz,Esasc,Esaxc,Exasc,Exaxc,Exayc,Exassc,Essaxc);
[Hca2]= Hca_gen(g1,g2,beta,cosx,cosy,cosz,Esasc,Esaxc,Exasc,Exaxc,Exayc,Exassc,Essaxc);
%%%%%%%%%%%%%%%%%% mater2 hamiltonian construction %%%%%%%%%%%%%%%%%%%%%%%%
[satr1,satr2,satr3,satr4]= satr_gen(inputss(4,1));
x=inputss(4,2);
y=inputss(4,3);
[a,D,Esa,Esc,Epa,Epc,Essa,Essc,Esasc,Esaxc,Exasc,Essaxc,Exassc,Exaxc,Exayc,Deltaa,Deltac,bandlineup_new]=energy_generator(x,y,inputss,bandlineup,satr1,satr2,satr3,satr4,ct);
[f1,f2,g1,g2,cosx,cosy,cosz,beta]=tao_gen(a,D,inputss,ct,k);
[Ha3]=Ha_gen(Esa,Epa,Essa,Deltaa,bandlineup_new);
[Hc3]=Hc_gen(Esc,Epc,Essc,Deltac,bandlineup_new);
[Hac3]= Hac_gen(f1,f2,beta,cosx,cosy,cosz,Esasc,Esaxc,Exasc,Exaxc,Exayc,Exassc,Essaxc);
[Hca3]= Hca_gen(g1,g2,beta,cosx,cosy,cosz,Esasc,Esaxc,Exasc,Exaxc,Exayc,Exassc,Essaxc);
%%%%%%%%%%%%%%%%%% interface2 hamiltonian construction %%%%%%%%%%%%%%%%%%%%
[satr1,satr2,satr3,satr4]= satr_gen(inputss(5,1));
x=inputss(5,2);
y=inputss(5,3);
[a,D,Esa,Esc,Epa,Epc,Essa,Essc,Esasc,Esaxc,Exasc,Essaxc,Exassc,Exaxc,Exayc,Deltaa,Deltac,bandlineup_new]=energy_generator(x,y,inputss,bandlineup,satr1,satr2,satr3,satr4,ct);
[f1,f2,g1,g2,cosx,cosy,cosz,beta]=tao_gen(a,D,inputss,ct,k);
[Ha4]=Ha_gen(Esa,Epa,Essa,Deltaa,bandlineup_new);
[Hc4]=Hc_gen(Esc,Epc,Essc,Deltac,bandlineup_new);
[Hac4]= Hac_gen(f1,f2,beta,cosx,cosy,cosz,Esasc,Esaxc,Exasc,Exaxc,Exayc,Exassc,Essaxc);
[Hca4]= Hca_gen(g1,g2,beta,cosx,cosy,cosz,Esasc,Esaxc,Exasc,Exaxc,Exayc,Exassc,Essaxc);
%%%%%%%%%%%%%%%%%%%%%%%%%%% Ha & Hc for interface %%%%%%%%%%%%%%%%%%%%%%%%%
  Ha2first=(0.5)*(Ha2+Ha1);
  Ha2second=(0.5)*(Ha2+Ha3);
  Hc2first=(0.5)*(Hc2+Hc1);
  Hc2second=(0.5)*(Hc2+Hc3);
  
  Ha4first=(0.5)*(Ha4+Ha3);
  Ha4second=(0.5)*(Ha4+Ha1);
  Hc4first=(0.5)*(Hc4+Hc3);
  Hc4second=(0.5)*(Hc4+Hc1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Hamiltonian matrix construction %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% for IF1=IF2 and made1=anion %%%%%%%%%%%%%%%%%%%%%%%%%%
 if (made1chiye==1)&&(inputss(3,1)==inputss(5,1))
     
     %%%%%%%%%%%%Ha,Hc
     Hsl(1:10,1:10)=Ha4second;
   
   
     for i=1:1:2*m-1
         if rem(i,2)==0
             Hsl(10*i+1:10*i+10,10*i+1:10*i+10)=Ha1;
         else
             Hsl(10*i+1:10*i+10,10*i+1:10*i+10)=Hc1;
         end
     end
     for i=2*m
         Hsl(10*i+1:10*i+10,10*i+1:10*i+10)=Ha2first;
     end
     for i=2*m+1
         Hsl(10*i+1:10*i+10,10*i+1:10*i+10)=Hc2second;
     end
      for i=2*m+2:1:atoms-2
         if rem(i,2)==0
             Hsl(10*i+1:10*i+10,10*i+1:10*i+10)=Ha3;
         else
             Hsl(10*i+1:10*i+10,10*i+1:10*i+10)=Hc3;
         end
      end
     for i=atoms-1
         Hsl(10*i+1:10*i+10,10*i+1:10*i+10)=Hc4first;
     end
     %%%%%%%%%%%% Hac1,Hca1
     for i=0:1:2*m-1
         if rem(i,2)==0
             Hsl(10*i+1:10*i+10,10*i+11:10*i+20)=Hac1;
             Hsl(10*i+11:10*i+20,10*i+1:10*i+10)=Hac1';
         else
             Hsl(10*i+1:10*i+10,10*i+11:10*i+20)=Hca1;
             Hsl(10*i+11:10*i+20,10*i+1:10*i+10)=Hca1'; 
         end
     end
     %%%%%%%%%%%% Hac2
     for i=2*m
         Hsl(10*i+1:10*i+10,10*i+11:10*i+20)=Hac2;
         Hsl(10*i+11:10*i+20,10*i+1:10*i+10)=Hac2';
     end

     %%%%%%%%%% Hac3,Hca3
     for i=2*m+1:1:atoms-2
         if rem(i,2)==0
             Hsl(10*i+1:10*i+10,10*i+11:10*i+20)=Hac3;
             Hsl(10*i+11:10*i+20,10*i+1:10*i+10)=Hac3';
         else
             Hsl(10*i+1:10*i+10,10*i+11:10*i+20)=Hca3;
             Hsl(10*i+11:10*i+20,10*i+1:10*i+10)=Hca3'; 
         end
     end
     %%%%%%%%% Hca4
     Hsl(1:10,10*(atoms-1)+1:10*atoms)=Hca4';
     Hsl(10*(atoms-1)+1:10*atoms,1:10)=Hca4;
 end

%%%%%%%%%%%%%%%%%%%% for IF1=IF2 and made1=cation %%%%%%%%%%%%%%%%%%%%%%%%%
  if (made1chiye==0)&&(inputss(3,1)==inputss(5,1))
     
     %%%%%%%%%%%%Ha,Hc
     Hsl(1:10,1:10)=Hc4second;
   
     for i=1:1:2*m-1
         if rem(i,2)==0
             Hsl(10*i+1:10*i+10,10*i+1:10*i+10)=Hc1;
         else
             Hsl(10*i+1:10*i+10,10*i+1:10*i+10)=Ha1;
         end
     end
     for i=2*m
         Hsl(10*i+1:10*i+10,10*i+1:10*i+10)=Hc2first;
     end
     for i=2*m+1
         Hsl(10*i+1:10*i+10,10*i+1:10*i+10)=Ha2second;
     end
      for i=2*m+2:1:atoms-2
         if rem(i,2)==0
             Hsl(10*i+1:10*i+10,10*i+1:10*i+10)=Hc3;
         else
             Hsl(10*i+1:10*i+10,10*i+1:10*i+10)=Ha3;
         end
      end
     for i=atoms-1
         Hsl(10*i+1:10*i+10,10*i+1:10*i+10)=Ha4first;
     end
     %%%%%%%%%%%% Hac1,Hca1
     for i=0:1:2*m-1
         if rem(i,2)==0
             Hsl(10*i+1:10*i+10,10*i+11:10*i+20)=Hca1;
             Hsl(10*i+11:10*i+20,10*i+1:10*i+10)=Hca1';
         else
             Hsl(10*i+1:10*i+10,10*i+11:10*i+20)=Hac1;
             Hsl(10*i+11:10*i+20,10*i+1:10*i+10)=Hac1'; 
         end
     end
     %%%%%%%%%%%% Hac2
     for i=2*m
         Hsl(10*i+1:10*i+10,10*i+11:10*i+20)=Hca2;
         Hsl(10*i+11:10*i+20,10*i+1:10*i+10)=Hca2';
     end

     %%%%%%%%%% Hac3,Hca3
     for i=2*m+1:1:atoms-2
         if rem(i,2)==0
             Hsl(10*i+1:10*i+10,10*i+11:10*i+20)=Hca3;
             Hsl(10*i+11:10*i+20,10*i+1:10*i+10)=Hca3';
         else
             Hsl(10*i+1:10*i+10,10*i+11:10*i+20)=Hac3;
             Hsl(10*i+11:10*i+20,10*i+1:10*i+10)=Hac3'; 
         end
     end
     %%%%%%%%% Hca4
     Hsl(1:10,10*(atoms-1)+1:10*atoms)=Hac4';
     Hsl(10*(atoms-1)+1:10*atoms,1:10)=Hac4;
 end
%%%%%%%%%%%%%%%%%%%% for IF1~=IF2 and made1=anion %%%%%%%%%%%%%%%%%%%%%%%%%
 if (made1chiye==1)&&(inputss(3,1)~=inputss(5,1))
     
     %%%%%%%%%%%%Ha,Hc
     Hsl(1:10,1:10)=Ha4second;
   
     for i=1:1:2*m
         if rem(i,2)==0
             Hsl(10*i+1:10*i+10,10*i+1:10*i+10)=Ha1;
         else
             Hsl(10*i+1:10*i+10,10*i+1:10*i+10)=Hc1;
         end
     end
     for i=2*m+1
         Hsl(10*i+1:10*i+10,10*i+1:10*i+10)=Hc2first;
     end
     for i=2*m+2
         Hsl(10*i+1:10*i+10,10*i+1:10*i+10)=Ha2second;
     end
      for i=2*m+3:1:atoms-2
         if rem(i,2)==0
             Hsl(10*i+1:10*i+10,10*i+1:10*i+10)=Ha3;
         else
             Hsl(10*i+1:10*i+10,10*i+1:10*i+10)=Hc3;
         end
      end
     for i=atoms-1
         Hsl(10*i+1:10*i+10,10*i+1:10*i+10)=Hc4first;
     end
     %%%%%%%%%%%% Hac1,Hca1
     for i=0:1:2*m
         if rem(i,2)==0
             Hsl(10*i+1:10*i+10,10*i+11:10*i+20)=Hac1;
             Hsl(10*i+11:10*i+20,10*i+1:10*i+10)=Hac1';
         else
             Hsl(10*i+1:10*i+10,10*i+11:10*i+20)=Hca1;
             Hsl(10*i+11:10*i+20,10*i+1:10*i+10)=Hca1'; 
         end
     end
     %%%%%%%%%%%% Hac2
     for i=2*m+1
         Hsl(10*i+1:10*i+10,10*i+11:10*i+20)=Hca2;
         Hsl(10*i+11:10*i+20,10*i+1:10*i+10)=Hca2';
     end

     %%%%%%%%%% Hac3,Hca3
     for i=2*m+2:1:atoms-2
         if rem(i,2)==0
             Hsl(10*i+1:10*i+10,10*i+11:10*i+20)=Hac3;
             Hsl(10*i+11:10*i+20,10*i+1:10*i+10)=Hac3';
         else
             Hsl(10*i+1:10*i+10,10*i+11:10*i+20)=Hca3;
             Hsl(10*i+11:10*i+20,10*i+1:10*i+10)=Hca3'; 
         end
     end
     %%%%%%%%% Hca4
     Hsl(1:10,10*(atoms-1)+1:10*atoms)=Hca4';
     Hsl(10*(atoms-1)+1:10*atoms,1:10)=Hca4;
 end
 %%%%%%%%%%%%%%%%%%%% for IF1~=IF2 and made1=cation %%%%%%%%%%%%%%%%%%%%%%%%%
 if (made1chiye==0)&&(inputss(3,1)~=inputss(5,1))
     
     %%%%%%%%%%%%Ha1,HC1
     Hsl(1:10,1:10)=Hc4second;
   
     for i=1:1:2*m
         if rem(i,2)==0
             Hsl(10*i+1:10*i+10,10*i+1:10*i+10)=Hc1;
         else
             Hsl(10*i+1:10*i+10,10*i+1:10*i+10)=Ha1;
         end
     end
     for i=2*m+1
         Hsl(10*i+1:10*i+10,10*i+1:10*i+10)=Ha2first;
     end
     for i=2*m+2
         Hsl(10*i+1:10*i+10,10*i+1:10*i+10)=Hc2second;
     end
      for i=2*m+3:1:atoms-2
         if rem(i,2)==0
             Hsl(10*i+1:10*i+10,10*i+1:10*i+10)=Hc3;
         else
             Hsl(10*i+1:10*i+10,10*i+1:10*i+10)=Ha3;
         end
      end
     for i=atoms-1
         Hsl(10*i+1:10*i+10,10*i+1:10*i+10)=Ha4first;
     end
     %%%%%%%%%%%% Hac1,Hca1
     for i=0:1:2*m
         if rem(i,2)==0
             Hsl(10*i+1:10*i+10,10*i+11:10*i+20)=Hca1;
             Hsl(10*i+11:10*i+20,10*i+1:10*i+10)=Hca1';
         else
             Hsl(10*i+1:10*i+10,10*i+11:10*i+20)=Hac1;
             Hsl(10*i+11:10*i+20,10*i+1:10*i+10)=Hac1'; 
         end
     end
     %%%%%%%%%%%% Hac2
     for i=2*m+1
         Hsl(10*i+1:10*i+10,10*i+11:10*i+20)=Hac2;
         Hsl(10*i+11:10*i+20,10*i+1:10*i+10)=Hac2';
     end

     %%%%%%%%%% Hac3,Hca3
     for i=2*m+2:1:atoms-2
         if rem(i,2)==0
             Hsl(10*i+1:10*i+10,10*i+11:10*i+20)=Hca3;
             Hsl(10*i+11:10*i+20,10*i+1:10*i+10)=Hca3';
         else
             Hsl(10*i+1:10*i+10,10*i+11:10*i+20)=Hac3;
             Hsl(10*i+11:10*i+20,10*i+1:10*i+10)=Hac3'; 
         end
     end
     %%%%%%%%% Hca4
     Hsl(1:10,10*(atoms-1)+1:10*atoms)=Hac4';
     Hsl(10*(atoms-1)+1:10*atoms,1:10)=Hac4;
 end
% zzz=zeros(10);
%  Hsl1=[ Hc4second     Hca1     zzz     zzz     zzz     zzz          zzz           zzz      zzz      zzz      zzz      Hac4'
%         Hca1'         Ha1      Hac1    zzz     zzz     zzz          zzz           zzz      zzz      zzz      zzz      zzz
%         zzz           Hac1'    Hc1     Hca1    zzz     zzz          zzz           zzz      zzz      zzz      zzz      zzz
%         zzz           zzz      Hca1'   Ha1     Hac1    zzz          zzz           zzz      zzz      zzz      zzz      zzz
%         zzz           zzz      zzz     Hac1'   Hc1     Hca1         zzz           zzz      zzz      zzz      zzz      zzz
%         zzz           zzz      zzz     zzz     Hca1'   Ha2first     Hac2          zzz      zzz      zzz      zzz      zzz
%         zzz           zzz      zzz     zzz     zzz     Hac2'        Hc2second     Hca3     zzz      zzz      zzz      zzz
%         zzz           zzz      zzz     zzz     zzz     zzz          Hca3'         Ha3      Hac3     zzz      zzz      zzz
%         zzz           zzz      zzz     zzz     zzz     zzz          zzz           Hac3'    Hc3      Hca3     zzz      zzz
%         zzz           zzz      zzz     zzz     zzz     zzz          zzz           zzz      Hca3'    Ha3      Hac3     zzz
%         zzz           zzz      zzz     zzz     zzz     zzz          zzz           zzz      zzz      Hac3'    Hc3      Hca3
%         Hac4          zzz      zzz     zzz     zzz     zzz          zzz           zzz      zzz      zzz      Hca3'    Ha4first];
%         
%     if Hsl1==Hsl
%         ppt=1;
%     else
%         ppt=0;
%     end
E=eig(Hsl);
end





























