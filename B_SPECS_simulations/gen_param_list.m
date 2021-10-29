% generate parameter list

name = 'param_list';
filename = strcat(name, '.txt');
fid = fopen(filename, 'wt');

%% Fit TS model, Fit 1
        for K_I = [0.5 1 2:2:20 50 100] 
            for K_A = [0.5 1 2:2:20 50 100]
                if(K_I < K_A)
                    for R = [1 10:10:100 200 500]
                        for BP = [10 100 200 500 1000 5000 10000]
                            for Vc = [1e-5 5e-4 1e-4 5e-3 1e-3 5e-2 1e-2]
                                param_string = strcat('3 1 30 500 1', {' '}, num2str(K_I), {' '}, num2str(K_A), {' '}, num2str(R), {' '}, num2str(BP), ...
                                    {' '},num2str(Vc), {' '},"1");
                                fprintf(fid, strcat(param_string,'\n'));  
                            end
                        end
                    end
                end
            end
        end



%% Fit Neumann model, Fit 1
% for K_BP = [0.5:0.5:3]
%     for p0 = [0:0.1:0.3]
%         for K_I = 0.1*2.^[0:1:10]
%             for K_A = 0.1*2.^[0:1:10]
% 
%                 if(K_I < K_A)
%                     param_string = strcat('3 0 30 500 1', {' '}, num2str(K_I), {' '}, num2str(K_A), {' '}, num2str(K_BP), {' '}, num2str(p0), {' '},"1");
%                     fprintf(fid, strcat(param_string,'\n'));
%                 end
%             end
%         end
%     end
% end

fclose(fid);
    
    