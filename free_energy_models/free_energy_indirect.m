% Free energy function used for maltose chemotaxis model developed in
%     Neumann, et al. The EMBO Journal 29, (2010).
%
% Noele Norris
%
%
function [f_total] = free_energy_indirect(parameters, L)

    KI = parameters(1);
    KA = parameters(2);
    K_BP = parameters(3);
    p0 = parameters(4);
    N2 = parameters(5);
   
    BP_bound = (p0+L./(L+K_BP));
    
    ratio1 = BP_bound/KI;
    ratio2 = BP_bound/KA;
        
    f_total =  N2*log((1+ratio1)./(1+ratio2));

end