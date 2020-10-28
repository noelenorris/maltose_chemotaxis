clear all;
close all;
clc;

seed = rng('shuffle');

%% 
%
%  SYSTEM IDENTIFICATION PROTOCOL: 
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
%       that best fit the experimental data. 
%
%       The code can be modified to find best fits for different chemotaxis
%       models.
%
%       To reduce run time, decrease number of iterations on line 80. Each
%       iteration takes approximately 10 seconds. 
%
%  Noele Norris
% 
%%

    % save all solutions
    name = strcat(datestr(now, 'yyyymmddHHMMSS'),'_fminconresults_TS')
    filename = strcat(name, '.txt');
    
    % access free energy functions
    addpath ../free_energy_models ../data
    
    
    % fmincon settings
        fun = @(x)objfun(x);
        b = [0;0];
        Aeq = [];
        beq = [];
        Aeq = [];
        beq = [];
        nonlcon = [];

        
       % TS model
           lb = [0; 0; 0; 0; 0; 0; 0; 0];
           ub = [500; 5000; 1000; 1000; 1000; 10000; 1; 50];


            %K_I must be smaller than K_A
           A = [0 0 1 -1 0 0 0 0; ...
                1 -1 0 0 0 0 0 0];
            
%        % Direct model
%               lb = [0; 0; 0; 0; 0];
%               ub = [500; 5000; 500; 5000; 50];
%               A = [0 0 1 -1 0; ...
%                    1 -1 0 0 0];

%         % Indirect model
%                lb = [0; 0; 0; 0; 0; 0; 0];
%                ub = [500; 5000; 100; 1000; 20; 1; 50];
%                A = [0 0 1 -1 0 0 0; ...
%                     1 -1 0 0 0 0 0];

%            % Linear approximation
%                 lb = [0; 0; 0; 0; 0; 0];
%                 ub = [500; 5000; 100; 100; 100; 50];
%                 A = [0 0 1 -1 0 0; ...
%                      1 -1 0 0 0 0];

       

        options = optimoptions('fmincon','Algorithm','interior-point','Display','off');


  
    f_save = 1000; %initialize goodness of fit
    x_save = [];
    eflag_save = [];
    outpt_save = [];

    %for manuscript, optimizations were run 1000 times
    for it=1:1:25
        
        it
        
        x0 = rand(length(ub),1).*[ub-lb]+lb;
       [x,f,eflag,outpt] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon, options);
       outpt;
       
       %save each optimization solution to file
       dlmwrite(filename, [f [x']], '-append', ...
                                             'delimiter', ' ');
       %track optimal solution
       if(f <= f_save)
           f_save = f
           x_save = x
           eflag_save = eflag;
           outpt_save = outpt;
       end
       
    end
    
    
    %%
    
    save(name);
       
                                 
% objective function to be minimized by fmincon                             
function y = objfun(x)
    load('data/data_gradients.mat')
    distributions= who('dist_data*');
    
    %Specify data bins
    bin_number = 100;
    n = 600/bin_number;
    X = (n/2):n:600;
    X = X(4:97);

    
    % optimization measure of fit J
    measure = 0;
    
    for i=1:1:length(distributions)

        % determine chemical environment used from data file
            dist_string = char(distributions{i});
            in = strfind(dist_string, '_');
            in2 = strfind(dist_string, 'maltose_')-1;
            in3 = strfind(dist_string, 'measp')-1;

            mal_conc_str = strrep(dist_string(in(2)+1:in2), 'p', '.');
            asp_conc_str = strrep(dist_string(in(3)+1:in3), 'p', '.');

            mal_conc = str2num(mal_conc_str);
            asp_conc = str2num(asp_conc_str);

        % construct chemical gradients using MeAsp and maltose 
        % concentrations specified from data file
            mal_conc = str2num(mal_conc_str);
            mal_left = (1000/1400)*(mal_conc);
            mal_right = (400/1400)*(mal_conc);
            mal_x = mal_left + (mal_right - mal_left)*(X/600);

            asp_conc = str2num(asp_conc_str);
            asp_left = (400/1400)*(asp_conc);
            asp_right = (1000/1400)*(asp_conc);
            asp_x = asp_left + (asp_right - asp_left)*(X/600);

        
        % Calculate analytical approximation of cell distributions 
        %using inputted parameter values
        
            parameters_maltose = [x(3:end)];
            parameters_measp = [x(1:2); x(end)];

            %free energy terms: 
            %MAKE MODIFICATIONS HERE TO TEST VARIOUS MODELS
            f_maltose = free_energy_TS(parameters_maltose, mal_x);
            %f_maltose = free_energy_direct(parameters_maltose, mal_x);
            
            f_measp = free_energy_direct(parameters_measp, asp_x);
            
            %predicted cell distribution
            f_total = f_maltose + f_measp;
            p = exp(f_total);
            p = p/sum(p);

            % goodness of fit measure,
            % J = summation of |p(x)-f(x)|^2,
            % where f(x) is the smoothed empirical distribution
            measure = measure + square_measure(eval(dist_string), p');
                                            
    end
    
    
    y = measure;
    
end
