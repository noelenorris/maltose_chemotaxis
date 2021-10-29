function [f_total] = free_energy_direct(parameters, L)

    KI = parameters(1); 
    KA = parameters(2);
    N1 = parameters(3);
    
    f_MeAsp = log((1+L/KI)./(1+L/KA));
   
    f_total = N1*f_MeAsp;

end