function [a,D,Esa,Esc,Epa,Epc,Essa,Essc,Esasc,Esaxc,Exasc,Essaxc,Exassc,Exaxc,Exayc,Deltaa,Deltac,bandlineup_new]=energy_generator(x,y,inputss,bandlineup,satr1,satr2,satr3,satr4,ct)
%a Gax In(1-x) Sby As(1-y)

a= x*y*ct(satr1,1)+ x*(1-y)*ct(satr2,1)+ (1-x)*y*ct(satr3,1)+ (1-x)*(1-y)*ct(satr4,1);
D= x*y*ct(satr1,inputss(1,2))+ x*(1-y)*ct(satr2,inputss(1,2))+ (1-x)*y*ct(satr3,inputss(1,2))+ (1-x)*(1-y)*ct(satr4,inputss(1,2));
Esa=x*y*ct(satr1,7)+ x*(1-y)*ct(satr2,7)+ (1-x)*y*ct(satr3,7)+ (1-x)*(1-y)*ct(satr4,7);
Esc=x*y*ct(satr1,8)+ x*(1-y)*ct(satr2,8)+ (1-x)*y*ct(satr3,8)+ (1-x)*(1-y)*ct(satr4,8);
Epa=x*y*ct(satr1,9)+ x*(1-y)*ct(satr2,9)+ (1-x)*y*ct(satr3,9)+ (1-x)*(1-y)*ct(satr4,9);
Epc=x*y*ct(satr1,10)+ x*(1-y)*ct(satr2,10)+ (1-x)*y*ct(satr3,10)+ (1-x)*(1-y)*ct(satr4,10);
Essa=x*y*ct(satr1,11)+ x*(1-y)*ct(satr2,11)+ (1-x)*y*ct(satr3,11)+ (1-x)*(1-y)*ct(satr4,11);
Essc=x*y*ct(satr1,12)+ x*(1-y)*ct(satr2,12)+ (1-x)*y*ct(satr3,12)+ (1-x)*(1-y)*ct(satr4,12);
Esasc=x*y*ct(satr1,13)+ x*(1-y)*ct(satr2,13)+ (1-x)*y*ct(satr3,13)+ (1-x)*(1-y)*ct(satr4,13);
Esaxc=x*y*ct(satr1,14)+ x*(1-y)*ct(satr2,14)+ (1-x)*y*ct(satr3,14)+ (1-x)*(1-y)*ct(satr4,14);
Exasc=x*y*ct(satr1,15)+ x*(1-y)*ct(satr2,15)+ (1-x)*y*ct(satr3,15)+ (1-x)*(1-y)*ct(satr4,15);
Essaxc=x*y*ct(satr1,16)+ x*(1-y)*ct(satr2,16)+ (1-x)*y*ct(satr3,16)+ (1-x)*(1-y)*ct(satr4,16);
Exassc=x*y*ct(satr1,17)+ x*(1-y)*ct(satr2,17)+ (1-x)*y*ct(satr3,17)+ (1-x)*(1-y)*ct(satr4,17);
Exaxc=x*y*ct(satr1,18)+ x*(1-y)*ct(satr2,18)+ (1-x)*y*ct(satr3,18)+ (1-x)*(1-y)*ct(satr4,18);
Exayc=x*y*ct(satr1,19)+ x*(1-y)*ct(satr2,19)+ (1-x)*y*ct(satr3,19)+ (1-x)*(1-y)*ct(satr4,19);
Deltaa=x*y*ct(satr1,20)+ x*(1-y)*ct(satr2,20)+ (1-x)*y*ct(satr3,20)+ (1-x)*(1-y)*ct(satr4,20);
Deltac=x*y*ct(satr1,21)+ x*(1-y)*ct(satr2,21)+ (1-x)*y*ct(satr3,21)+ (1-x)*(1-y)*ct(satr4,21);
bandlineup_new=x*y*bandlineup(satr1,1)+ x*(1-y)*bandlineup(satr2,1)+ (1-x)*y*bandlineup(satr3,1)+ (1-x)*(1-y)*bandlineup(satr4,1);
end


