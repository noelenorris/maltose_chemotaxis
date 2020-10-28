% Free energy function used in baseline transport-and-sensing chemotaxis model
%
% Noele Norris
%
function [f_total] = free_energy_TS_not_porin_limited(parameters, L)

    KI = parameters(1); 
    KA = parameters(2);
    R = parameters(3);
    BP = parameters(4);

    Vc = parameters(5);
    
    M = parameters(6);
    
    Kc = 100;
    Kbp = 2;
    Kp = 10^4;
    Vp = 1;
    %Vc = 10e-4;
    
    Lp = L;
    %LBP = -(BP*(BP*Kp*Vc - (BP^2*Kp^2*Vc^2 + 2*BP^2*Kp*Lext*Vc^2 - 2*BP^2*Kp*Lext*Vc*Vp + BP^2*Lext.^2*Vc^2 - 2*BP^2*Lext.^2*Vc*Vp + BP^2*Lext.^2*Vp^2 + 2*BP*Kbp*Kc*Kp*Vc*Vp + 6*BP*Kbp*Kc*Lext*Vc*Vp + 2*BP*Kbp*Kc*Lext*Vp^2 + 2*BP*Kc*Kp^2*Vc^2 + 4*BP*Kc*Kp*Lext*Vc^2 - 4*BP*Kc*Kp*Lext*Vc*Vp + 2*BP*Kc*Lext.^2*Vc^2 - 4*BP*Kc*Lext.^2*Vc*Vp + 2*BP*Kc*Lext.^2*Vp^2 + Kbp^2*Kc^2*Vp^2 + 2*Kbp*Kc^2*Kp*Vc*Vp + 6*Kbp*Kc^2*Lext*Vc*Vp + 2*Kbp*Kc^2*Lext*Vp^2 + Kc^2*Kp^2*Vc^2 + 2*Kc^2*Kp*Lext*Vc^2 - 2*Kc^2*Kp*Lext*Vc*Vp + Kc^2*Lext.^2*Vc^2 - 2*Kc^2*Lext.^2*Vc*Vp + Kc^2*Lext.^2*Vp^2).^(1/2) + BP*Lext*Vc - BP*Lext*Vp + Kbp*Kc*Vp + Kc*Kp*Vc + Kc*Lext*Vc - Kc*Lext*Vp))./((BP^2*Kp^2*Vc^2 + 2*BP^2*Kp*Lext*Vc^2 - 2*BP^2*Kp*Lext*Vc*Vp + BP^2*Lext.^2*Vc^2 - 2*BP^2*Lext.^2*Vc*Vp + BP^2*Lext.^2*Vp^2 + 2*BP*Kbp*Kc*Kp*Vc*Vp + 6*BP*Kbp*Kc*Lext*Vc*Vp + 2*BP*Kbp*Kc*Lext*Vp^2 + 2*BP*Kc*Kp^2*Vc^2 + 4*BP*Kc*Kp*Lext*Vc^2 - 4*BP*Kc*Kp*Lext*Vc*Vp + 2*BP*Kc*Lext.^2*Vc^2 - 4*BP*Kc*Lext.^2*Vc*Vp + 2*BP*Kc*Lext.^2*Vp^2 + Kbp^2*Kc^2*Vp^2 + 2*Kbp*Kc^2*Kp*Vc*Vp + 6*Kbp*Kc^2*Lext*Vc*Vp + 2*Kbp*Kc^2*Lext*Vp^2 + Kc^2*Kp^2*Vc^2 + 2*Kc^2*Kp*Lext*Vc^2 - 2*Kc^2*Kp*Lext*Vc*Vp + Kc^2*Lext.^2*Vc^2 - 2*Kc^2*Lext.^2*Vc*Vp + Kc^2*Lext.^2*Vp^2).^(1/2) + 2*BP*Kbp*Vc + 2*BP*Kbp*Vp - BP*Kp*Vc - BP*Lext*Vc + BP*Lext*Vp + 2*Kbp*Kc*Vc + Kbp*Kc*Vp - Kc*Kp*Vc - Kc*Lext*Vc + Kc*Lext*Vp);
    LBP = BP*Lp./(Kbp+Lp);
    
    % concentration of bound receptors are solutions to quadratic equation
        RL_I = ((R+LBP+KI)-sqrt((R+LBP+KI).^2 - 4*R*LBP))/2;
        C_I = RL_I./(R-RL_I);

        RL_A = ((R+LBP+KA)-sqrt((R+LBP+KA).^2 - 4*R*LBP))/2;
        C_A = RL_A./(R-RL_A);
    
    
    f_total = M*log((1+C_I)./(1+C_A));


end