function [Hc]=Hc_gen(Esc,Epc,Essc,Deltac,bandlineup_new)
hcc=[Esc+bandlineup_new         0                              0                         0                     0
       0                        Epc+bandlineup_new             -1i*(Deltac/3)            0                     0
       0                        1i*(Deltac/3)                  Epc+bandlineup_new        0                     0
       0                        0                              0                         Epc+bandlineup_new    0
       0                        0                              0                         0                     Essc+bandlineup_new];
       
      
      
      
  hcso=[0     0                0                    0                 0
        0     0                0                    (Deltac/3)        0
        0     0                0                    -1i*(Deltac/3)    0
        0     -(Deltac/3)      1i*(Deltac/3)        0                 0
        0     0                0                    0                 0];
        
      Hc=[hcc     hcso
          hcso'   conj(hcc)];
end







