clear all;
close all;
clc;

%% 
%
%  AGENT-BASED SIMULATIONS: 
%
%       This code calls sim_HMWC to run simulations corresponding to the
%       various experimental conditions.
%
%       Each suite of simulations executes in approximately 5 minutes.
%
%  Noele Norris
% 
%%

    % access free energy functions
    addpath ../free_energy_models ../data
    
    bin_number = 100;
    n = 600/bin_number;
    X = (n/2):n:600;
    X2 = X(4:97);
    
    % specify model and parameteters to simulate in text file
    to_sim_fileID = fopen('simulations_to_do.txt');
    sim_parameters = fgetl(to_sim_fileID);
    
     
   % specify which simulations to run using corresponding
   % experimental data files
   load('data/data_gradients.mat')
   distributions= who('dist_data*');

tic
while(sim_parameters ~= -1)
    toc
    
        % use parameter values from .txt file
        parameters = str2num(sim_parameters);
        parameter_string = sim_parameters;
        
        FileName = strcat('sim_distributions', parameter_string,'.mat');


        for i=1:1:length(distributions)

            dist_string = char(distributions{i})

            
            % determine chemical environment used from data file
                dist_string = char(distributions{i});
                in = strfind(dist_string, '_');
                in2 = strfind(dist_string, 'maltose_')-1;
                in3 = strfind(dist_string, 'measp')-1;

                mal_conc_str = strrep(dist_string(in(2)+1:in2), 'p', '.');
                asp_conc_str = strrep(dist_string(in(3)+1:in3), 'p', '.');

                mal_conc = str2num(mal_conc_str);
                asp_conc = str2num(asp_conc_str);

            % run simulation
            [ss_distribution, dis] = sim_HMWC(parameters, 0, asp_conc, mal_conc, 0, 1000, 24000);

            % plot simulation steady-state distribution vs. emperical
            % distribution
            close all;
            figure;
            hold on
            h{i} = plot(X, ss_distribution);
            plot(X2, eval(dist_string), '--');
            title(strrep(strcat(dist_string, ' ', parameter_string), '_', ' '));
            hold off
            drawnow

            % save simulated steady-state distribution to file
            v = genvarname(strrep(dist_string, 'data', 'sim'));
            eval([v ' = ss_distribution;']);
            
            if exist(FileName, 'file')
                save(FileName, v,'-append');
            else
                save(FileName,v);
            end
        end
       

    sim_parameters = fgetl(to_sim_fileID);
end
        
                                    
        

                      