clear all
close all
clc

%% 
%
%  AGENT-BASED SIMULATIONS: 
%
%       This code calls sim_HMWC to run simulations corresponding to the
%       gradient conditions and model specified in text file.
%
%
%  Noele Norris
% 
%%
        

    % access free energy functions
    addpath ../free_energy_models ../data

    % specify model and parameteters to simulate in text file
    file_sims = 'simulations_to_do_BP.txt';

    fid=fopen(file_sims, 'r');

    bin_number = 100;
    n = 600/bin_number;
    X = (n/2):n:600;
    X2 = X(4:97);
   
   
    fid=fopen(file_sims, 'r');
    FileName = string(fgetl(fid));
    sim_parameters = fgetl(fid)
    
    while(sim_parameters ~= -1)

        parameters = str2num(sim_parameters);

        num_sims = str2num(fgetl(fid));

        for i=1:1:num_sims

            tic

            concentrations = str2num(string(fgetl(fid)))

            mal_conc = concentrations(1);
            asp_conc = concentrations(2);
            
            mal_string = num2str(mal_conc)
            asp_string = num2str(asp_conc)
            var_name = strcat('sim_', mal_string,'_',asp_string)

            % run simulation
            [ss_distribution, dis] = sim_HMWC(parameters, 0, asp_conc, mal_conc, 0, 2*1000, 2*2*24000);
            
            % save simulated steady-state distribution to file
            v = genvarname(var_name);
            eval([v ' = ss_distribution;']);

            toc

        end
        
        save(strcat(FileName, '_',sim_parameters,'.mat'))
        
        sim_parameters = fgetl(fid)


    end
        
                                    
        

                      