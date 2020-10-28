% Free energy function used in baseline transport-and-sensing chemotaxis model
%
% Noele Norris
%
function [f_total] = free_energy_TS_linear(parameters, L)


    KI = parameters(1); 
    KA = parameters(2);
    R = parameters(3);
    M = parameters(4);
    
    LBP = L; %L./(alpha+L);
    
    % concentration of bound receptors are solutions to quadratic equation
        RL_I = ((R+LBP+KI)-sqrt((R+LBP+KI).^2 - 4*R*LBP))/2;
        C_I = RL_I./(R-RL_I);

        RL_A = ((R+LBP+KA)-sqrt((R+LBP+KA).^2 - 4*R*LBP))/2;
        C_A = RL_A./(R-RL_A);
    
    
    f_total = M*log((1+C_I)./(1+C_A));


end