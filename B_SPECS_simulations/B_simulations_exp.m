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
function [] = B_simulations_exp(start_sim, end_sim)


        % access free energy functions
        addpath ../free_energy_models ../data
        
        file_sims = 'simulations_to_do_Fits.txt';
        
        fid=fopen(file_sims, 'r');

        bin_number = 100;
        n = 600/bin_number;
        X = (n/2):n:600;
        X2 = X(4:97);


       % specify which simulations to run using corresponding
       % experimental data files
       
       load('data/data_maltose.mat');
       %load('data/data_growth_sim.mat');
       
       distributions= who('dist_data*');

       sim_results = cell(length(distributions));

    tic

    min_measure = 1000;

    counter = start_sim
    
    list_results = [];

    while(counter <= end_sim)
        toc
        
            fid=fopen(file_sims, 'r');
            sim_parameters = textscan(fid,'%s',1,'delimiter','\n', 'headerlines',counter-1);
            sim_parameters = cell2mat(sim_parameters{1});
            
            % use parameter values from .txt file
            parameters = str2num(sim_parameters);
            parameter_string = sim_parameters;

            FileName = strcat('sim_distributions', parameter_string,'.mat');


            measure = 0;

            for i=1:1:length(distributions)

                dist_string = char(distributions{i});


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
                [ss_distribution, dis] = sim_HMWC(parameters, 0, asp_conc, mal_conc, 0, 2*1000, 2*2*24000);

                %%FIX!!!!
                if(asp_conc == 0)
                    CMC = -1;
                else
                    CMC = 1;
                end
                
                if(CMC > 0)
                    f = fit(X(2:end-1)', ss_distribution(2:end-1), 'power2');
                    fitted_distribution = f.a*X.^f.b+f.c;
                else
                    f = fit(600-X(2:end-1)', ss_distribution(2:end-1), 'power2');
                    fitted_distribution = f.a*(600-X).^f.b+f.c;
                end
                


                measure = measure + square_measure(eval(dist_string), fitted_distribution(4:97)');


%                 % plot simulation steady-state distribution vs. emperical
%                 % distribution
%                 %close all;
%                 figure()
%                 hold on
%                 h{i} = plot(X(2:end-1), ss_distribution(2:end-1));
%                 plot(X2, eval(dist_string), '--');
%                 plot(X2, fitted_distribution(4:97))
%                 title(strrep(strcat(dist_string, ' ', parameter_string), '_', ' '));
%                 hold off
%                 drawnow

                % save simulated steady-state distribution to file
                v = genvarname(strrep(dist_string, 'data', 'sim'));
                eval([v ' = ss_distribution;']);

            end
            
            list_results = [list_results; strcat(num2str(measure), {' '}, parameter_string)];
             save(FileName)
             measure
            
            if(measure < min_measure)
                min_measure = measure
                min_sim = parameter_string
                
            end

            counter = counter + 1
   
    end
    
    save(strcat('BestFits_', num2str(start_sim), '_', num2str(end_sim)), 'list_results');
end
        
                                    
        

                      