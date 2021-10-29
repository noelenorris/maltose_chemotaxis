
%% 
%
%  SIMULATION FUNCTION: 
%
%       This function runs agent-based simulations using a modified version
%   of SPECS that incorporates HMWC model and transport-and-sensing model for
%   maltose chemotaxis. It returns the steady-state cell distribution.
%
%       Original SPECS simulation model from: 
%           Jiang, et al. "Quantitative Modeling of Escherichia coli 
%           Chemotactic Motion in Environments Varying in Space and Time".
%           PLoS Computational Biology, 2010.
%
%
%
%  Noele Norris, Filippo Menolascina
%
%%


function [ss_distribution, cell_dis] = sim_HMWC(parameters, asp_left, asp_right, ...
    maltose_left, maltose_right, num_cells, num_iterations)
   
    
    %% SPECS parameters
        num_receptors = parameters(1);
        model = parameters(2);
        measp_parameters = parameters(3:5);
        maltose_parameters = parameters(6:end);
        
        % observed cell characteristics
        speed = 20;                 % cell speed (um/s)
        tau_0 = 0.2;                % tumble time (s)
        delta_rot_diff = 30*pi/180;        % average directional change in 1s
    
        % chemotaxis parameters
        n_Tar = num_receptors;      % number of Tar receptors used in free energy calculation
        m0 = 1.45;                   % initial methylation level 
        m_saturation = 4;           % max methylation level 
        alpha = 1.875;              % methylation factor
        k_r = 0.005;                % methylation rate (1/s)
        k_b = 2*k_r;                % demehylation rate (1/s)
        H = 10.3;                   % Hill coefficient
        a_0 =k_r/(k_b+k_r);			% steady-state activity
        c = 2;
    
    %% Simulation run parameters
            

        dt = 0.1;                   %time step (s)
        dt_save = 1;                %time step for saving (s)

        total_x = 600;              %channel length (um)
        dx = 6;                

        total_y = 1000;             %channel width (um)
        dy = 5;
        
        total_z = 100;              % channel depth (um)
        dz = 5;

        tumble_time = tau_0/dt;     %time spent tumbling 
        
        % save cell distribution at specified time interval dt_save
        num_save = num_iterations*(dt/dt_save);
        cell_dis = zeros(total_x/dx, total_y/dy, num_save); 
        

    
        %% Environment parameters 
        
        %substrate concentrations at boundaries of test channel:
            %aspartate
            asp_ch_left=asp_left+(400/1400)*(asp_right-asp_left);
            asp_ch_right=asp_left+(1000/1400)*(asp_right-asp_left);

            %maltose
            mal_ch_left = maltose_left+(400/1400)*(maltose_right-maltose_left);
            mal_ch_right = maltose_left+(1000/1400)*(maltose_right-maltose_left);

        
        

    %% Cell initialization

        %initialize each cell to random position in channel
        x = total_x*rand(num_cells, 1);                         
        y = total_y*rand(num_cells, 1);
        z = total_z*rand(num_cells, 1);
        
        %initialize methylation level
        m = m0*ones(num_cells, 1);
        
        %initialize all cells to running (0: tumbling, 1: running)
        runners = ones(num_cells, 1);                   
        
        %initialize tumble counter of cells to zero
        tumble_counter = zeros(num_cells, 1);   
        
        %initialize cells to random Euler orientation
        phi = pi*rand(num_cells, 1);
        theta = 2*pi*rand(num_cells, 1);
        



    %% Main Time Loop
    for it = 1:1:num_iterations
        
        % concentration of substrate experienced by each cell
         L_asp = (asp_ch_right-asp_ch_left)*x/total_x + asp_ch_left;
         L_mal = (mal_ch_right-mal_ch_left)*x/total_x + mal_ch_left;
         
        
        %Free energy terms
        f_methylation_tar = n_Tar.*alpha.*(m0 - m);

        if(model == 0) %indirect-binding model
            f_MeAsp = free_energy_direct(measp_parameters, L_asp);
            f_MeAsp_0 = free_energy_direct(measp_parameters, 0);
            f_maltose = free_energy_indirect(maltose_parameters, L_mal);
            f_maltose_0 = free_energy_indirect(maltose_parameters, 0);
        elseif(model == 1) %transport-and-sensing chemotaxis model
           f_MeAsp = free_energy_direct(measp_parameters, L_asp);
           f_MeAsp_0 = free_energy_direct(measp_parameters, 0);
           f_maltose= free_energy_TS(maltose_parameters, L_mal);
           f_maltose_0 = free_energy_TS(maltose_parameters, 0);
        elseif(model == 2) %linear approx of transport-and-sensing model
           f_MeAsp = free_energy_direct(measp_parameters, L_asp);
           f_MeAsp_0 = free_energy_direct(measp_parameters, 0);
           f_maltose= free_energy_TS_linear(maltose_parameters, L_mal);
           f_maltose_0 = free_energy_TS_linear(maltose_parameters, 0);
        elseif(model == 3) %transport-and-sensing, not porin-limited
           f_MeAsp = free_energy_direct(measp_parameters, L_asp);
           f_MeAsp_0 = free_energy_direct(measp_parameters, 0);
           f_maltose= free_energy_TS_not_porin_limited(maltose_parameters, L_mal);
           f_maltose_0 = free_energy_TS_not_porin_limited(maltose_parameters, 0);
        else
           error('Unknown model.');
        end
        

        %current activity levels
        a = 1 ./ (1 + c*exp(f_methylation_tar + n_Tar*(f_MeAsp + f_maltose-f_MeAsp_0-f_maltose_0)));
    
        %methylation level must be bestween 0 and m_saturation
        m = min(max(m +(k_r.*(1-a)-k_b.*a).*dt,0),m_saturation);
    
        %probability of tumbling
        probability_of_tumble = 0.25*(dt/tau_0).*(a./a_0).^H;
    
        %determine cells that switch from run to tumble
        idx_r = find(runners);
        idx_t = find(rand(size(idx_r))<probability_of_tumble(idx_r));
        runners(idx_r(idx_t)) = 0;
        tumble_counter(idx_r(idx_t)) = 0;
        
        %update location of running cells
        idx_r = find(runners);

        
        x(idx_r) = x(idx_r) + speed*dt*cos(phi(idx_r)).*cos(theta(idx_r));
        y(idx_r) = y(idx_r) + speed*dt*cos(phi(idx_r)).*sin(theta(idx_r));
        z(idx_r) = z(idx_r) + speed*dt*sin(phi(idx_r));
        
         theta(idx_r) = mod(theta(idx_r)+randn(size(idx_r))*delta_rot_diff*dt.^0.5,2*pi);
         phi(idx_r) = mod(phi(idx_r)+randn(size(idx_r))*delta_rot_diff*dt.^0.5,pi);
        
        % check boundaries, changing orientation of cells at boundary 
            idx_b = find(x(idx_r)<0);
            x(idx_r(idx_b)) = 0;
            %orientation(idx_r(idx_b)) = -90+randn(size(idx_r(idx_b)))*180;
            
            idx_b = find(x(idx_r)>total_x); 
            x(idx_r(idx_b)) = total_x;
            %orientation(idx_r(idx_b)) = 90+randn(size(idx_r(idx_b)))*180;
            
            idx_b = find(y(idx_r)<0); 
            y(idx_r(idx_b)) = 0;
            %orientation(idx_r(idx_b)) = 0+randn(size(idx_r(idx_b)))*180;
            
            idx_b = find(y(idx_r)>total_y);
            y(idx_r(idx_b)) = total_y;
            %orientation(idx_r(idx_b)) = 180+randn(size(idx_r(idx_b)))*180;
            
            idx_b = find(z(idx_r)<0); 
            z(idx_r(idx_b)) = 0;
            
            idx_b = find(y(idx_r)>total_z);
            z(idx_r(idx_b)) = total_z;
            
        %update counter of tumbling cells
            % cells still tumbling
            idx_t=find(and(~runners,tumble_counter < tumble_time));
            tumble_counter(idx_t) = tumble_counter(idx_t) + 1; 
            
            % cells finished tumbling
            idx_r=find(and(~runners,tumble_counter >= tumble_time));
            runners(idx_r)=1;
            
            phi(idx_r) = pi*rand(size(idx_r));
            theta(idx_r) = 2*pi*rand(size(idx_r));
            %orientation(idx_r)=rand(size(idx_r))*360;

        %Save current cell distribution
        if(mod(dt*it, dt_save) == 0)
                  idx = find((x <= 600).*(x>=0));
                  x_save = x(idx);
                  y_save = y(idx);
                 [dis,~]=hist3([x_save,y_save], {[dx/2:dx:(total_x-dx/2)]',[dy/2:dy:(total_y-dy/2)]'});
                 cell_dis(:, :, it*dt/dt_save) = dis;
        end
             
    end

    % return averaged steady-state distribution
    dist = cell_dis(:,:,(num_save/2):num_save);
    ss_distribution = sum(sum(dist,3),2)/sum(sum(sum(dist)));
    
end