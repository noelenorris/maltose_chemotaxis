
%% 
%
%  Fitting transport-and-sensing model to FRET data: 
%
%
%       !IMPORTANT! PLEASE READ: This code will not run  without 
%       MATLAB's Global Optimization Toolbox. To install the toolbox, 
%       click: Home --> Add-Ons --> Get Add-Ons. 
%       Search for "Global Optimization Toolbox" and click on toolbox 
%       then "Install". Follow the prompts to finish installation. 
%       (In older versions of MATLAB, you can obtain the toolbox 
%        by re-running the installer.)
%
%
%       This code uses MATLAB's fmincon to determine the parameter values
%       that best fit the experimental FRET reporter assays from:
%       Neumann, et al. "Differences in signalling by directly and 
%       indirectly binding ligands in bacterial chemotaxis." The EMBO 
%       Journal 29, 3484–3495 (2010).
%
%
%       To reduce run time, decrease number of iterations.
%
%  Noele Norris
% 
%%

clear all;
close all;
clc;

seed = rng('shuffle');
addpath ../free_energy_models ../data
load('Neumann_FRET_data.mat')


% Inputs to fmincon optimization program
    % Pass values of FRET assay dose response and dynamic range data
    % (Figures 1A and 2A in Neumann et al.)
    fun = @(x)objfun(x, maltose_dose_response, maltose_range, measp_dose_response, measp_range);
    b = [0; 0];
    Aeq = [];
    beq = [];
    nonlcon = [];

   % lower and upper bounds
   lb = [20;   200;   0;   0;  10;   500;  1e-4;    1;   6];
   ub = [40;   600; 100; 100;  50;  1000;  1e-3;    1;   6];

    
    %K_I must be smaller than K_A
    A = [0 0 1 -1 0 0 0 0 0;
         1 -1 0 0 0 0 0 0 0];

    options = optimoptions('fmincon','Algorithm','interior-point','Display','off');
        
f_save = 1000
x_save = [];
           
           
for it=1:1:200
    it
    
    x0 = rand(length(ub),1).*[ub-lb]+lb;

   [x,f,eflag,outpt] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon, options);
   outpt;

   %save all found solutions
   if(eflag > 0)
       f_save = f;
       x_save = [x_save; [f x']];
       eflag_save = eflag;
       outpt_save = outpt;
   end
   %toc
end
           
           

function y = objfun(x, maltose_dose_response, maltose_range, measp_dose_response, measp_range)


    param = x;
    measp_param = param([1:2 end-1]);
    mal_param = param(3:end-1);
    n = param(end);

    c = 2;
    A0 = 1/(1+c);

    L = 10^3*[0 1e-6 1e-5 3e-5 1e-4 3e-4 1e-3 3e-3 1e-2 3e-2 1e-1 3e-1 1];

    f_mal = free_energy_TS(mal_param, L);
    f_mal_0 = free_energy_TS(mal_param, 0);
    f_mal_3 = free_energy_TS(mal_param, L/3);
    f_measp = free_energy_direct(measp_param, L);
    f_measp_3 = free_energy_direct(measp_param, L/3);
    f_measp_0 = free_energy_direct(measp_param, 0);
    f_measp_100 = free_energy_direct(measp_param, 100);
        
    A_sat = 1/(1+c*exp(n*(f_measp_100 - f_measp_0)));
    A_mal = 1./(1+c*exp(n*(f_mal - f_mal_0)));
    A_measp = 1./(1+c*exp(n*(f_measp - f_measp_0)));
    A_mal_3 = 1./(1+c*exp(n*(f_mal - f_mal_3)));
    A_measp_3 = 1./(1+c*exp(n*(f_measp - f_measp_3)));

    
    % activity and range values obtained from model
        a_mal = (A_mal-A_sat)/(A0-A_sat);
        a_measp = (A_measp-A_sat)/(A0-A_sat);

        d_mal = (A0-A_mal_3)/(A0-A_sat);
        d_measp = (A0-A_measp_3)/(A0-A_sat);
        
    % activity and range values from FRET data
        indx_a_mal = [1 2 3 4 5 6 7 9 11];
        a_mal_exp = maltose_dose_response(:,2)';
        
        indx_a_measp = [1 5 6 7 8 9 11];
        a_measp_exp = measp_dose_response(:,2)';

        indx_d_mal = [1 4 5 6 7 8 9 10];
        d_mal_exp = maltose_range(:,2)';
        
        indx_d_measp = [1 5 6 7 8 9 10 11 12 13];
        d_measp_exp = measp_range(:,2)';
        
        
    % objective function: minimize sum of squares of dynamic range
        y = sum((d_mal(indx_d_mal)-d_mal_exp).^2) + sum((d_measp(indx_d_measp)-d_measp_exp).^2);
    

end
