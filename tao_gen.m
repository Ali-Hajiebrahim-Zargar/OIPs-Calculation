
function [f1,f2,g1,g2,cosx,cosy,cosz,beta]=tao_gen(a,D,inputss,ct,k)
exx=(ct(inputss(1,1),1)/a)-1;
eyy=exx;
ezz=-D*exx;
beta= 3/((1+exx)^2+(1+eyy)^2+(1+ezz)^2);
cosx=sqrt(beta/3)*(1+exx);
cosy=sqrt(beta/3)*(1+eyy);
cosz=sqrt(beta/3)*(1+ezz);
tao1=a/4 *[1+exx,1+eyy,1+ezz];
tao2=a/4 *[-1-exx,-1-eyy,1+ezz];
tao3=a/4 *[1+exx,-1-eyy,-1-ezz];
tao4=a/4 *[-1-exx,1+eyy,-1-ezz];
x1=sum(k.*tao1);
x2=sum(k.*tao2);
x3=sum(k.*tao3);
x4=sum(k.*tao4);
g1=(exp(-1i*x3)+exp(-1i*x4))/4;
g2=(exp(-1i*x3)-exp(-1i*x4))/4;
f1=(exp(1i*x1)+exp(1i*x2))/4;
f2=(exp(1i*x1)-exp(1i*x2))/4;
end

