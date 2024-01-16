
function [sub,dir,matter1,cation1,anion1,interface1,cation_interface1,anion_interface1,matter2,cation2,anion2,interface2,cation_interface2,anion_interface2]= type_generator(inputss)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% type of substrate %%%%%%%%%%%%%%%%%%%%%%
if (inputss(1,1))==1
    sub='Si';
    elseif (inputss(1,1))==2
    sub='Ge';
    elseif (inputss(1,1))==3
    sub='GaAs';
    elseif (inputss(1,1))==4
    sub='AlAs';
    elseif (inputss(1,1))==5
    sub='InAs';
    elseif (inputss(1,1))==6
    sub='GaP';
    elseif (inputss(1,1))==7
    sub='AlP';
    elseif (inputss(1,1))==8
    sub='InP';
    elseif (inputss(1,1))==9
    sub='GaSb';
    elseif (inputss(1,1))==10
    sub='AlSb';
    elseif (inputss(1,1))==11
    sub='InSb';
    elseif (inputss(1,1))==12
    sub='ZnSe';
    elseif (inputss(1,1))==13
    sub='ZnS';
    elseif (inputss(1,1))==14
    sub='ZnTe';
    elseif (inputss(1,1))==15
    sub='CdTe';
    elseif (inputss(1,1))==16
    sub='HgTe';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% type of direction%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if inputss(1,2)==2
    dir='D001';
elseif inputss(1,2)==3
    dir='D110';
elseif inputss(1,2)==4
    dir='D111';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%type of matter 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%                    
if (inputss(2,1))==1
    matter1='Si';
    elseif (inputss(2,1))==2
    matter1='Ge';
    elseif (inputss(2,1))==3
    matter1='GaAs';
    cation1='Ga';
    anion1='As';
%     satr1=17;
%     satr2=17;
%     satr3=17;
%     satr4=3;
    elseif (inputss(2,1))==4
    matter1='AlAs';
    cation1='Al';
    anion1='As';
%     satr1=17;
%     satr2=17;
%     satr3=17;
%     satr4=4;
    elseif (inputss(2,1))==5
    matter1='InAs';
    cation1='In';
    anion1='As';
%     satr1=17;
%     satr2=17;
%     satr3=17;
%     satr4=5;
    elseif (inputss(2,1))==6
    matter1='GaP';
    cation1='Ga';
    anion1='P';
%     satr1=17;
%     satr2=17;
%     satr3=17;
%     satr4=6;
    elseif (inputss(2,1))==7
    matter1='AlP';
    cation1='Al';
    anion1='P';
%     satr1=17;
%     satr2=17;
%     satr3=17;
%     satr4=7;
    elseif (inputss(2,1))==8
    matter1='InP';
    cation1='In';
    anion1='P';
%     satr1=17;
%     satr2=17;
%     satr3=17;
%     satr4=8;
    elseif (inputss(2,1))==9
    matter1='GaSb';
    cation1='Ga';
    anion1='Sb';
%     satr1=17;
%     satr2=17;
%     satr3=17;
%     satr4=9;
    elseif (inputss(2,1))==10
    matter1='AlSb';
    cation1='Al';
    anion1='Sb';
%     satr1=17;
%     satr2=17;
%     satr3=17;
%     satr4=10;
    elseif (inputss(2,1))==11
    matter1='InSb';
    cation1='In';
    anion1='Sb';
%     satr1=17;
%     satr2=17;
%     satr3=17;
%     satr4=11;
    elseif (inputss(2,1))==12
    matter1='ZnSe';
    elseif (inputss(2,1))==13
    matter1='ZnS';
    elseif (inputss(2,1))==14
    matter1='ZnTe';
    elseif (inputss(2,1))==15
    matter1='CdTe';
    elseif (inputss(2,1))==16
    matter1='HgTe';
    elseif (inputss(2,1))==17
    matter1='AlGaP';
    cation1='AlGa';
    anion1='P';
%     satr1=17;
%     satr2=7;
%     satr3=17;
%     satr4=6;
    elseif (inputss(2,1))==18
    matter1='AlGaAs';
    cation1='AlGa';
    anion1='As';
%     satr1=17;
%     satr2=4;
%     satr3=17;
%     satr4=3;
    elseif (inputss(2,1))==19
    matter1='AlGaSb';
    cation1='AlGa';
    anion1='Sb';
%     satr1=17;
%     satr2=10;
%     satr3=17;
%     satr4=9;
    elseif (inputss(2,1))==20
    matter1='AlInP';
    cation1='AlIn';
    anion1='P';
%     satr1=17;
%     satr2=7;
%     satr3=17;
%     satr4=8;
    elseif (inputss(2,1))==21
    matter1='AlInAs';
    cation1='AlIn';
    anion1='As';
%     satr1=17;
%     satr2=4;
%     satr3=17;
%     satr4=5;
    elseif (inputss(2,1))==22
    matter1='AlInSb';
    cation1='AlIn';
    anion1='Sb';
%     satr1=17;
%     satr2=10;
%     satr3=17;
%     satr4=11;
    elseif (inputss(2,1))==23
    matter1='GaInP';
    cation1='GaIn';
    anion1='P';
%     satr1=17;
%     satr2=6;
%     satr3=17;
%     satr4=8;
    elseif (inputss(2,1))==24
    matter1='GaInAs';
    cation1='GaIn';
    anion1='As';
%     satr1=17;
%     satr2=3;
%     satr3=17;
%     satr4=5;
    elseif (inputss(2,1))==25
    matter1='GaInSb';
    cation1='GaIn';
    anion1='Sb';
%     satr1=17;
%     satr2=9;
%     satr3=17;
%     satr4=11;
    elseif (inputss(2,1))==26
    matter1='PAsAl';
    cation1='Al';
    anion1='PAs';
%     satr1=17;
%     satr2=7;
%     satr3=17;
%     satr4=4;
    elseif (inputss(2,1))==27
    matter1='PAsGa';
    cation1='Ga';
    anion1='PAs';
%     satr1=17;
%     satr2=6;
%     satr3=17;
%     satr4=3;
    elseif (inputss(2,1))==28
    matter1='PAsIn';
    cation1='In';
    anion1='PAs';
%     satr1=17;
%     satr2=8;
%     satr3=17;
%     satr4=5;
    elseif (inputss(2,1))==29
    matter1='PSbAl';
    cation1='Al';
    anion1='PSb';
%     satr1=17;
%     satr2=7;
%     satr3=17;
%     satr4=10;
    elseif (inputss(2,1))==30
    matter1='PSbGa';
    cation1='Ga';
    anion1='PSb';
%     satr1=17;
%     satr2=6;
%     satr3=17;
%     satr4=9;
    elseif (inputss(2,1))==31
    matter1='PSbIn';
    cation1='In';
    anion1='PSb';
%     satr1=17;
%     satr2=8;
%     satr3=17;
%     satr4=11;
    elseif (inputss(2,1))==32
    matter1='AsSbAl';
    cation1='Al';
    anion1='AsSb';
%     satr1=17;
%     satr2=4;
%     satr3=17;
%     satr4=10;
    elseif (inputss(2,1))==33
    matter1='AsSbGa';
    cation1='Ga';
    anion1='AsSb';
%     satr1=17;
%     satr2=3;
%     satr3=17;
%     satr4=9;
    elseif (inputss(2,1))==34
    matter1='AsSbIn';
    cation1='In';
    anion1='AsSb';
%     satr1=17;
%     satr2=5;
%     satr3=17;
%     satr4=11;
    elseif (inputss(2,1))==35
    matter1='AlGaPAs';
    cation1='AlGa';
    anion1='PAs';
%     satr1=7;
%     satr2=4;
%     satr3=6;
%     satr4=3;
    elseif (inputss(2,1))==36
    matter1='AlGaPSb';
    cation1='AlGa';
    anion1='PSb';
%     satr1=7;
%     satr2=10;
%     satr3=6;
%     satr4=9;
    elseif (inputss(2,1))==37
    matter1='AlGaAsSb';
    cation1='AlGa';
    anion1='AsSb';
%     satr1=4;
%     satr2=10;
%     satr3=3;
%     satr4=9;
    elseif (inputss(2,1))==38
    matter1='AlInPAs';
    cation1='AlIn';
    anion1='PAs';
%     satr1=7;
%     satr2=4;
%     satr3=8;
%     satr4=5;
    elseif (inputss(2,1))==39
    matter1='AlInPSb';
    cation1='AlIn';
    anion1='PSb';
%     satr1=7;
%     satr2=10;
%     satr3=8;
%     satr4=11;
    elseif (inputss(2,1))==40
    matter1='AlInAsSb';
    cation1='AlIn';
    anion1='AsSb';
%     satr1=4;
%     satr2=10;
%     satr3=5;
%     satr4=11;
    elseif (inputss(2,1))==41
    matter1='GaInPAs';
    cation1='GaIn';
    anion1='PAs';
%     satr1=6;
%     satr2=3;
%     satr3=8;
%     satr4=5;
    elseif (inputss(2,1))==42
    matter1='GaInPSb';
    cation1='GaIn';
    anion1='PSb';
%     satr1=6;
%     satr2=9;
%     satr3=8;
%     satr4=11;
    elseif (inputss(2,1))==43
    matter1='GaInAsSb';
    cation1='GaIn';
    anion1='AsSb';
%     satr1=3;
%     satr2=9;
%     satr3=5;
%     satr4=11;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%type of interface 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%                    
if (inputss(3,1))==1
    interface1='Si';
    elseif (inputss(3,1))==2
    interface1='Ge';
    elseif (inputss(3,1))==3
    interface1='GaAs';
    cation_interface1='Ga';
    anion_interface1='As';
    elseif (inputss(3,1))==4
    interface1='AlAs';
    cation_interface1='Al';
    anion_interface1='As';
    elseif (inputss(3,1))==5
    interface1='InAs';
    cation_interface1='In';
    anion_interface1='As';
    elseif (inputss(3,1))==6
    interface1='GaP';
    cation_interface1='Ga';
    anion_interface1='P';
    elseif (inputss(3,1))==7
    interface1='AlP';
    cation_interface1='Al';
    anion_interface1='P';
    elseif (inputss(3,1))==8
    interface1='InP';
    cation_interface1='In';
    anion_interface1='P';
    elseif (inputss(3,1))==9
    interface1='GaSb';
    cation_interface1='Ga';
    anion_interface1='Sb';
    elseif (inputss(3,1))==10
    interface1='AlSb';
    cation_interface1='Al';
    anion_interface1='Sb';
    elseif (inputss(3,1))==11
    interface1='InSb';
    cation_interface1='In';
    anion_interface1='Sb';
    elseif (inputss(3,1))==12
    interface1='ZnSe';
    elseif (inputss(3,1))==13
    interface1='ZnS';
    elseif (inputss(3,1))==14
    interface1='ZnTe';
    elseif (inputss(3,1))==15
    interface1='CdTe';
    elseif (inputss(3,1))==16
    interface1='HgTe';
    elseif (inputss(3,1))==17
    interface1='AlGaP';
    cation_interface1='AlGa';
    anion_interface1='P';
    elseif (inputss(3,1))==18
    interface1='AlGaAs';
    cation_interface1='AlGa';
    anion_interface1='As';
    elseif (inputss(3,1))==19
    interface1='AlGaSb';
    cation_interface1='AlGa';
    anion_interface1='Sb';
    elseif (inputss(3,1))==20
    interface1='AlInP';
    cation_interface1='AlIn';
    anion_interface1='P';
    elseif (inputss(3,1))==21
    interface1='AlInAs';
    cation_interface1='AlIn';
    anion_interface1='As';
    elseif (inputss(3,1))==22
    interface1='AlInSb';
    cation_interface1='AlIn';
    anion_interface1='Sb';
    elseif (inputss(3,1))==23
    interface1='GaInP';
    cation_interface1='GaIn';
    anion_interface1='P';
    elseif (inputss(3,1))==24
    interface1='GaInAs';
    cation_interface1='GaIn';
    anion_interface1='As';
    elseif (inputss(3,1))==25
    interface1='GaInSb';
    cation_interface1='GaIn';
    anion_interface1='Sb';
    elseif (inputss(3,1))==26
    interface1='PAsAl';
    cation_interface1='Al';
    anion_interface1='PAs';
    elseif (inputss(3,1))==27
    interface1='PAsGa';
    cation_interface1='Ga';
    anion_interface1='PAs';
    elseif (inputss(3,1))==28
    interface1='PAsIn';
    cation_interface1='In';
    anion_interface1='PAs';
    elseif (inputss(3,1))==29
    interface1='PSbAl';
    cation_interface1='Al';
    anion_interface1='PSb';
    elseif (inputss(3,1))==30
    interface1='PSbGa';
    cation_interface1='Ga';
    anion_interface1='PSb';
    elseif (inputss(3,1))==31
    interface1='PSbIn';
    cation_interface1='In';
    anion_interface1='PSb';
    elseif (inputss(3,1))==32
    interface1='AsSbAl';
    cation_interface1='Al';
    anion_interface1='AsSb';
    elseif (inputss(3,1))==33
    interface1='AsSbGa';
    cation_interface1='Ga';
    anion_interface1='AsSb';
    elseif (inputss(3,1))==34
    interface1='AsSbIn';
    cation_interface1='In';
    anion_interface1='AsSb';
    elseif (inputss(3,1))==35
    interface1='AlGaPAs';
    cation_interface1='AlGa';
    anion_interface1='PAs';
    elseif (inputss(3,1))==36
    interface1='AlGaPSb';
    cation_interface1='AlGa';
    anion_interface1='PSb';
    elseif (inputss(3,1))==37
    interface1='AlGaAsSb';
    cation_interface1='AlGa';
    anion_interface1='AsSb';
    elseif (inputss(3,1))==38
    interface1='AlInPAs';
    cation_interface1='AlIn';
    anion_interface1='PAs';
    elseif (inputss(3,1))==39
    interface1='AlInPSb';
    cation_interface1='AlIn';
    anion_interface1='PSb';
    elseif (inputss(3,1))==40
    interface1='AlInAsSb';
    cation_interface1='AlIn';
    anion_interface1='AsSb';
    elseif (inputss(3,1))==41
    interface1='GaInPAs';
    cation_interface1='GaIn';
    anion_interface1='PAs';
    elseif (inputss(3,1))==42
    interface1='GaInPSb';
    cation_interface1='GaIn';
    anion_interface1='PSb';
    elseif (inputss(3,1))==43
    interface1='GaInAsSb';
    cation_interface1='GaIn';
    anion_interface1='AsSb';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%type of matter 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%                    
if (inputss(4,1))==1
    matter2='Si';
    elseif (inputss(4,1))==2
    matter2='Ge';
    elseif (inputss(4,1))==3
    matter2='GaAs';
    cation2='Ga';
    anion2='As';
    elseif (inputss(4,1))==4
    matter2='AlAs';
    cation2='Al';
    anion2='As';
    elseif (inputss(4,1))==5
    matter2='InAs';
    cation2='In';
    anion2='As';
    elseif (inputss(4,1))==6
    matter2='GaP';
    cation2='Ga';
    anion2='P';
    elseif (inputss(4,1))==7
    matter2='AlP';
    cation2='Al';
    anion2='P';
    elseif (inputss(4,1))==8
    matter2='InP';
    cation2='In';
    anion2='P';
    elseif (inputss(4,1))==9
    matter2='GaSb';
    cation2='Ga';
    anion2='Sb';
    elseif (inputss(4,1))==10
    matter2='AlSb';
    cation2='Al';
    anion2='Sb';
    elseif (inputss(4,1))==11
    matter2='InSb';
    cation2='In';
    anion2='Sb';
    elseif (inputss(4,1))==12
    matter2='ZnSe';
    elseif (inputss(4,1))==13
    matter2='ZnS';
    elseif (inputss(4,1))==14
    matter2='ZnTe';
    elseif (inputss(4,1))==15
    matter2='CdTe';
    elseif (inputss(4,1))==16
    matter2='HgTe';
    elseif (inputss(4,1))==17
    matter2='AlGaP';
    cation2='AlGa';
    anion2='P';
    elseif (inputss(4,1))==18
    matter2='AlGaAs';
    cation2='AlGa';
    anion2='As';
    elseif (inputss(4,1))==19
    matter2='AlGaSb';
    cation2='AlGa';
    anion2='Sb';
    elseif (inputss(4,1))==20
    matter2='AlInP';
    cation2='AlIn';
    anion2='P';
    elseif (inputss(4,1))==21
    matter2='AlInAs';
    cation2='AlIn';
    anion2='As';
    elseif (inputss(4,1))==22
    matter2='AlInSb';
    cation2='AlIn';
    anion2='Sb';
    elseif (inputss(4,1))==23
    matter2='GaInP';
    cation2='GaIn';
    anion2='P';
    elseif (inputss(4,1))==24
    matter2='GaInAs';
    cation2='GaIn';
    anion2='As';
    elseif (inputss(4,1))==25
    matter2='GaInSb';
    cation2='GaIn';
    anion2='Sb';
    elseif (inputss(4,1))==26
    matter2='PAsAl';
    cation2='Al';
    anion2='PAs';
    elseif (inputss(4,1))==27
    matter2='PAsGa';
    cation2='Ga';
    anion2='PAs';
    elseif (inputss(4,1))==28
    matter2='PAsIn';
    cation2='In';
    anion2='PAs';
    elseif (inputss(4,1))==29
    matter2='PSbAl';
    cation2='Al';
    anion2='PSb';
    elseif (inputss(4,1))==30
    matter2='PSbGa';
    cation2='Ga';
    anion2='PSb';
    elseif (inputss(4,1))==31
    matter2='PSbIn';
    cation2='In';
    anion2='PSb';
    elseif (inputss(4,1))==32
    matter2='AsSbAl';
    cation2='Al';
    anion2='AsSb';
    elseif (inputss(4,1))==33
    matter2='AsSbGa';
    cation2='Ga';
    anion2='AsSb';
    elseif (inputss(4,1))==34
    matter2='AsSbIn';
    cation2='In';
    anion2='AsSb';
    elseif (inputss(4,1))==35
    matter2='AlGaPAs';
    cation2='AlGa';
    anion2='PAs';
    elseif (inputss(4,1))==36
    matter2='AlGaPSb';
    cation2='AlGa';
    anion2='PSb';
    elseif (inputss(4,1))==37
    matter2='AlGaAsSb';
    cation2='AlGa';
    anion2='AsSb';
    elseif (inputss(4,1))==38
    matter2='AlInPAs';
    cation2='AlIn';
    anion2='PAs';
    elseif (inputss(4,1))==39
    matter2='AlInPSb';
    cation2='AlIn';
    anion2='PSb';
    elseif (inputss(4,1))==40
    matter2='AlInAsSb';
    cation2='AlIn';
    anion2='AsSb';
    elseif (inputss(4,1))==41
    matter2='GaInPAs';
    cation2='GaIn';
    anion2='PAs';
    elseif (inputss(4,1))==42
    matter2='GaInPSb';
    cation2='GaIn';
    anion2='PSb';
    elseif (inputss(4,1))==43
    matter2='GaInAsSb';
    cation2='GaIn';
    anion2='AsSb';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%type of interface 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%                    
if (inputss(5,1))==1
    interface2='Si';
    elseif (inputss(5,1))==2
    interface2='Ge';
    elseif (inputss(5,1))==3
    interface2='GaAs';
    cation_interface2='Ga';
    anion_interface2='As';
    elseif (inputss(5,1))==4
    interface2='AlAs';
    cation_interface2='Al';
    anion_interface2='As';
    elseif (inputss(5,1))==5
    interface2='InAs';
    cation_interface2='In';
    anion_interface2='As';
    elseif (inputss(5,1))==6
    interface2='GaP';
    cation_interface2='Ga';
    anion_interface2='P';
    elseif (inputss(5,1))==7
    interface2='AlP';
    cation_interface2='Al';
    anion_interface2='P';
    elseif (inputss(5,1))==8
    interface2='InP';
    cation_interface2='In';
    anion_interface2='P';
    elseif (inputss(5,1))==9
    interface2='GaSb';
    cation_interface2='Ga';
    anion_interface2='Sb';
    elseif (inputss(5,1))==10
    interface2='AlSb';
    cation_interface2='Al';
    anion_interface2='Sb';
    elseif (inputss(5,1))==11
    interface2='InSb';
    cation_interface2='In';
    anion_interface2='Sb';
    elseif (inputss(5,1))==12
    interface2='ZnSe';
    elseif (inputss(5,1))==13
    interface2='ZnS';
    elseif (inputss(5,1))==14
    interface2='ZnTe';
    elseif (inputss(5,1))==15
    interface2='CdTe';
    elseif (inputss(5,1))==16
    interface2='HgTe';
    elseif (inputss(5,1))==17
    interface2='AlGaP';
    cation_interface2='AlGa';
    anion_interface2='P';
    elseif (inputss(5,1))==18
    interface2='AlGaAs';
    cation_interface2='AlGa';
    anion_interface2='As';
    elseif (inputss(5,1))==19
    interface2='AlGaSb';
    cation_interface2='AlGa';
    anion_interface2='Sb';
    elseif (inputss(5,1))==20
    interface2='AlInP';
    cation_interface2='AlIn';
    anion_interface2='P';
    elseif (inputss(5,1))==21
    interface2='AlInAs';
    cation_interface2='AlIn';
    anion_interface2='As';
    elseif (inputss(5,1))==22
    interface2='AlInSb';
    cation_interface2='AlIn';
    anion_interface2='Sb';
    elseif (inputss(5,1))==23
    interface2='GaInP';
    cation_interface2='GaIn';
    anion_interface2='P';
    elseif (inputss(5,1))==24
    interface2='GaInAs';
    cation_interface2='GaIn';
    anion_interface2='As';
    elseif (inputss(5,1))==25
    interface2='GaInSb';
    cation_interface2='GaIn';
    anion_interface2='Sb';
    elseif (inputss(5,1))==26
    interface2='PAsAl';
    cation_interface2='Al';
    anion_interface2='PAs';
    elseif (inputss(5,1))==27
    interface2='PAsGa';
    cation_interface2='Ga';
    anion_interface2='PAs';
    elseif (inputss(5,1))==28
    interface2='PAsIn';
    cation_interface2='In';
    anion_interface2='PAs';
    elseif (inputss(5,1))==29
    interface2='PSbAl';
    cation_interface2='Al';
    anion_interface2='PSb';
    elseif (inputss(5,1))==30
    interface2='PSbGa';
    cation_interface2='Ga';
    anion_interface2='PSb';
    elseif (inputss(5,1))==31
    interface2='PSbIn';
    cation_interface2='In';
    anion_interface2='PSb';
    elseif (inputss(5,1))==32
    interface2='AsSbAl';
    cation_interface2='Al';
    anion_interface2='AsSb';
    elseif (inputss(5,1))==33
    interface2='AsSbGa';
    cation_interface2='Ga';
    anion_interface2='AsSb';
    elseif (inputss(5,1))==34
    interface2='AsSbIn';
    cation_interface2='In';
    anion_interface2='AsSb';
    elseif (inputss(5,1))==35
    interface2='AlGaPAs';
    cation_interface2='AlGa';
    anion_interface2='PAs';
    elseif (inputss(5,1))==36
    interface2='AlGaPSb';
    cation_interface2='AlGa';
    anion_interface2='PSb';
    elseif (inputss(5,1))==37
    interface2='AlGaAsSb';
    cation_interface2='AlGa';
    anion_interface2='AsSb';
    elseif (inputss(5,1))==38
    interface2='AlInPAs';
    cation_interface2='AlIn';
    anion_interface2='PAs';
    elseif (inputss(5,1))==39
    interface2='AlInPSb';
    cation_interface2='AlIn';
    anion_interface2='PSb';
    elseif (inputss(5,1))==40
    interface2='AlInAsSb';
    cation_interface2='AlIn';
    anion_interface2='AsSb';
    elseif (inputss(5,1))==41
    interface2='GaInPAs';
    cation_interface2='GaIn';
    anion_interface2='PAs';
    elseif (inputss(5,1))==42
    interface2='GaInPSb';
    cation_interface2='GaIn';
    anion_interface2='PSb';
    elseif (inputss(5,1))==43
    interface2='GaInAsSb';
    cation_interface2='GaIn';
    anion_interface2='AsSb';
end



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%