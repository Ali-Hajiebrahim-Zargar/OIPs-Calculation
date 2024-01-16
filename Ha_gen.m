function [Ha]=Ha_gen(Esa,Epa,Essa,Deltaa,bandlineup_new)
haa=[Esa+bandlineup_new         0                              0                         0                     0
       0                        Epa+bandlineup_new             -1i*(Deltaa/3)            0                     0
       0                        1i*(Deltaa/3)                  Epa+bandlineup_new        0                     0
       0                        0                              0                         Epa+bandlineup_new    0
       0                        0                              0                         0                     Essa+bandlineup_new];
       
      
      
      
  haso=[0     0                0                    0                 0
        0     0                0                    (Deltaa/3)        0
        0     0                0                    -1i*(Deltaa/3)    0
        0     -(Deltaa/3)      1i*(Deltaa/3)        0                 0
        0     0                0                    0                 0];
        
      Ha=[haa     haso
          haso'   conj(haa)];
end







