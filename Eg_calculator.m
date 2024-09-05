
function [Eg,wave_lenght,mismatch_Sl]=Eg_calculator(ct,m,n,inputss,kx,ky,kz)
% tic
[sub,dir,matter1,cation1,anion1,interface1,cation_interface1,anion_interface1,matter2,cation2,anion2,interface2,cation_interface2,anion_interface2]= type_generator(inputss);
%%%%%%%%%%%%%%%%%%%%%%%%%%% number of atoms in SLS %%%%%%%%%%%%%%%%%%%%%%%
if (inputss(3,1)==inputss(5,1))
    atoms=(m+n+1)*2;
else
    atoms=(m+n+2)*2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%% determine starter matter %%%%%%%%%%%%%%%%%%%%%
if anion1==anion_interface2
    made1chiye=1;
elseif cation1==cation_interface2
    made1chiye=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% constants matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% calculation of bandline up %%%%%%%%%%%%%%%%%%%%%
    bandlineup=zeros(17,1);
for i=1:1:17
    bandlineup(i,1)=ct(i,6)+(ct(i,5)/3);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% mismatch and a prependecular calculation%%%%%%%%%%%%%%%%%
    a_prependecular=zeros(1,5);
for i=2:1:5
    a_prependecular(1,i)= (ct(inputss(i,1),1))*(1-ct(inputss(i,1),inputss(1,2))*((ct(inputss(1,1),1)/ct(inputss(i,1),1)-1)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
    if inputss(3,1)== inputss(5,1)
        type=1;
    else
        type=0;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (type)
                        
        SL_const= m*a_prependecular(1,2)+0.5*a_prependecular(1,3)+n*a_prependecular(1,4)+0.5* a_prependecular(1,5);
        a_SL=SL_const/(m+n+1);
        mismatch_Sl= (ct(inputss(1,1),1)- a_SL)/ct(inputss(1,1),1);

    else

        SL_const= (m+0.5)*a_prependecular(1,2)+0.5*a_prependecular(1,3)+(n+0.5)*a_prependecular(1,4)+0.5* a_prependecular(1,5);
        a_SL=SL_const/(m+n+2);
        mismatch_Sl= (ct(inputss(1,1),1)- a_SL)/ct(inputss(1,1),1);
    end
    az=SL_const/2;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    
% for j=0:1:s
%     kx = 0;
%     ky = 0 ;
%     kz = 0 ;
    [~,E]=SL_Hamiltonian(ct,inputss,atoms,made1chiye,m,kx,ky,kz,bandlineup); 
    Ec=E((4*atoms)+1,1);
    Ev=E(4*atoms,1);
%     x1 = sqrt(((pi/ct(inputss(1,1),1))-kx)^2+((pi/ct(inputss(1,1),1))-ky)^2+((pi/az)-kz)^2) ;
%     Y(j+1,1) = x1 ;
%     Y(j+1,2:(10*atoms+1)) = E' ;
% end
% 
% 
% for j=s+1:1:2*s
%     kx = ((pi)/ct(inputss(1,1),1))*((j-s)/s) ;
%     ky = 0 ;
%     kz=0 ;
%     [~,E]=SL_Hamiltonian(ct,inputss,atoms,made1chiye,m,kx,ky,kz,bandlineup);
%     Ec(j+1,1)=E((4*atoms)+1,1);
%     Ev(j+1,1)=E(4*atoms,1);
%     x2 = sqrt((kx)^2)+x1 ;
%     Y(j+1,1) = x2 ;
%     Y(j+1,2:(10*atoms+1)) = E' ;
% end


% for j=2*s+1:1:3*s
%     kx = (pi)/ct(inputss(1,1),1) ;
%     ky = ((pi)/ct(inputss(1,1),1))*((j-2*s)/s) ;
%     kz = 0;
%     [~,E]=SL_Hamiltonian(ct,inputss,atoms,made1chiye,m,kx,ky,kz,bandlineup);
%     Ec(j+1,1)=E((4*atoms)+1,1);
%     Ev(j+1,1)=E(4*atoms,1);
%     x3 = sqrt(((ky)^2))+ x2 ;
%     Y(j+1,1) = x3 ;
%     Y(j+1,2:(10*atoms+1)) = E' ;
% end

% 
% for j=3*s+1:1:4*s
%     kx = ((pi)/(ct(inputss(1,1),1)))- ((pi)/(ct(inputss(1,1),1)))*((j-3*s)/s) ;
%     ky = ((pi)/(ct(inputss(1,1),1)))- ((pi)/(ct(inputss(1,1),1)))*((j-3*s)/s) ;
%     kz = 0 ;
%     [~,E]=SL_Hamiltonian(ct,inputss,atoms,made1chiye,m,kx,ky,kz,bandlineup);
%     Ec(j+1,1)=E((4*atoms)+1,1);
%     Ev(j+1,1)=E(4*atoms,1);
%     x4 = sqrt(((((pi)/(ct(inputss(1,1),1)))-kx)^2)+((((pi)/(ct(inputss(1,1),1)))-ky)^2))+x3 ;
%     Y(j+1,1) = x4 ;
%     Y(j+1,2:(10*atoms+1)) = E' ;
% end

% k=[0 0 0];
% [Hsl,E]=SL_Hamiltonian(ct,inputss,atoms,made1chiye,m,kx,ky,kz,bandlineup);
% csvwrite('ali.csv',Hsl);
% for j = 4*s+1:1:5*s
%     kx = 0 ;
%     ky = 0 ;
%     kz = ((pi)/az)*((j-4*s)/s) ;
%     [~,E]=SL_Hamiltonian(ct,inputss,atoms,made1chiye,m,kx,ky,kz,bandlineup);
%     Ec(j+1,1)=E((4*atoms)+1,1);
%     Ev(j+1,1)=E(4*atoms,1);
%     x5 = sqrt((kz)^2)+x4 ;
%     Y(j+1,1) = x5 ;
%     Y(j+1,2:(10*atoms+1)) = E' ;
% end



% Ecc=min(Ec);
% Evv=max(Ev);
Eg=abs(Ec-Ev); 

wave_lenght=(1.239/Eg);
% toc

end
