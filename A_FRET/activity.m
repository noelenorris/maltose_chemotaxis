%% 
%
%  Plotting FRET data and fits: 
%
%
%
%       This code plots the FRET data from  Neumann, et al.  along with
%       model predictions of both Neumann's fitted indirect-binding model
%       and our fitted transport-and-sensing model.
%
%  Noele Norris
% 
%%


%close all
%clear all
%clc

addpath ../free_energy_models ../data

load('Neumann_FRET_data.mat')

cmap = ...
      [0.27 0.47 0.67; ...
       0.40 0.80 0.93; ...
       0.13 0.53 0.20; ...
       0.80 0.73 0.27; ...
       0.93 0.40 0.47; ...
       0.93 0.69 0.13; ...
       0.67 0.20 0.47;
       0.73 0.73 0.73; ...
       ];

   
%% BEST FITS

%%
% optimal parameter set for TS model, obtained from activity_optimize.m
% Because overdetermined, constrained n_Tar = 6.
param1 = [27.5 365	14.4	49.7	12.6	101	1.09E-05	1	6];
param2 = [27.4777237600431,363.181970211153,394.343183411292,2043.28496935524,21.0178084717777,1286.05277153775,1.00007592510641e-05,1,6]
param3 = [27.5 365	14.4	49.7	12.6	101*4	1.09E-05	1	6];


measp_param1 = param1([1:2 end-1]);
mal_param1 = param1(3:end-1);
n1 = param1(end);

measp_param2 = param2([1:2 end-1]);
mal_param2 = param2(3:end-1);
n2 = param2(end);

measp_param3 = param3([1:2 end-1]);
mal_param3 = param3(3:end-1);
n3 = param3(end);


% Neumann's indirect-binding model parameter values
mal_param_N = [0.4 6 2 0.1 1];
measp_param_N = [30 500 1];
n_N = 6;



c = 2;
A0 = 1/(1+c);


L = [10.^[-3:0.1:4]];

% Free energy and activity levels for indirect-binding model
     f_mal_N = free_energy_indirect(mal_param_N, L);
     f_mal_3_N = free_energy_indirect(mal_param_N, L/3);
     f_mal_0_N = free_energy_indirect(mal_param_N, 0);
    
    f_measp_N = free_energy_direct(measp_param_N, L);
    f_measp_3_N = free_energy_direct(measp_param_N, L/3);
    f_measp_0_N = free_energy_direct(measp_param_N, 0);
    f_measp_100_N = free_energy_direct(measp_param_N, 100);
    A_sat_N = 1/(1+c*exp(n_N*(f_measp_100_N - f_measp_0_N)));
    A_mal_N = 1./(1+c*exp(n_N*(f_mal_N - f_mal_0_N)));
    A_mal_3_N = 1./(1+c*exp(n_N*(f_mal_N-f_mal_3_N)));
    A_measp_3_N = 1./(1+c*exp(n_N*(f_measp_N - f_measp_3_N)));
    A_measp_N = 1./(1+c*exp(n_N*(f_measp_N - f_measp_0_N)));

% Free energy and activity levels for TS model, Fit 1
    f_mal = free_energy_TS(mal_param1, L);
    f_mal_0 = free_energy_TS(mal_param1, 0);
    f_mal_3 = free_energy_TS(mal_param1, L/3);
    f_measp = free_energy_direct(measp_param1, L);
    f_measp_3 = free_energy_direct(measp_param1, L/3);
    f_measp_0 = free_energy_direct(measp_param1, 0);
    f_measp_100 = free_energy_direct(measp_param1, 100);

    A_sat = 1/(1+c*exp(n1*(f_measp_100 - f_measp_0)));
    A_mal = 1./(1+c*exp(n1*(f_mal - f_mal_0)));
    A_measp = 1./(1+c*exp(n1*(f_measp - f_measp_0)));
    A_mal_3 = 1./(1+c*exp(n1*(f_mal - f_mal_3)));
    A_measp_3 = 1./(1+c*exp(n1*(f_measp - f_measp_3)));
    
% Free energy and activity levels for TS model, Fit 2
    f_mal2 = free_energy_TS(mal_param2, L);
    f_mal_0_2 = free_energy_TS(mal_param2, 0);
    f_mal_3_2 = free_energy_TS(mal_param2, L/3);
    f_measp2 = free_energy_direct(measp_param2, L);
    f_measp_3_2 = free_energy_direct(measp_param2, L/3);
    f_measp_0_2 = free_energy_direct(measp_param2, 0);
    f_measp_100_2 = free_energy_direct(measp_param2, 100);

    A_sat2 = 1/(1+c*exp(n2*(f_measp_100_2 - f_measp_0_2)));
    A_mal2 = 1./(1+c*exp(n2*(f_mal2 - f_mal_0_2)));
    A_measp2 = 1./(1+c*exp(n2*(f_measp2 - f_measp_0_2)));
    A_mal_3_2 = 1./(1+c*exp(n2*(f_mal2 - f_mal_3_2)));
    A_measp_3_2 = 1./(1+c*exp(n2*(f_measp2 - f_measp_3_2)));
    
% Free energy and activity levels for TS model, Fit 3
    f_mal3 = free_energy_TS(mal_param3, L);
    f_mal_0_3 = free_energy_TS(mal_param3, 0);
    f_mal_3_3 = free_energy_TS(mal_param3, L/3);
    f_measp3 = free_energy_direct(measp_param3, L);
    f_measp_3_3 = free_energy_direct(measp_param3, L/3);
    f_measp_0_3 = free_energy_direct(measp_param3, 0);
    f_measp_100_3 = free_energy_direct(measp_param3, 100);

    A_sat3 = 1/(1+c*exp(n3*(f_measp_100_3 - f_measp_0_3)));
    A_mal3 = 1./(1+c*exp(n3*(f_mal3 - f_mal_0_3)));
    A_measp3 = 1./(1+c*exp(n3*(f_measp3 - f_measp_0_3)));
    A_mal_3_3 = 1./(1+c*exp(n3*(f_mal3 - f_mal_3_3)));
    A_measp_3_3 = 1./(1+c*exp(n3*(f_measp3 - f_measp_3_3)));



figure()
semilogx(10^(-3)*L, (A_mal_N-A_sat_N)/(A0-A_sat_N), ':', 'LineWidth', 8, 'Color', 'black')
hold on
semilogx(10^(-3)*L, (A_measp_N-A_sat_N)/(A0-A_sat_N), ':', 'LineWidth', 8, 'Color', cmap(8,:))
semilogx(10^(-3)*L, (A_mal-A_sat)/(A0-A_sat), 'LineWidth', 6, 'Color', cmap(2,:))
semilogx(10^(-3)*L, (A_measp-A_sat)/(A0-A_sat), 'LineWidth', 6, 'Color', cmap(1,:))
semilogx(10^(-3)*L, (A_mal2-A_sat2)/(A0-A_sat2), 'LineWidth', 4, 'Color', cmap(3,:))
semilogx(10^(-3)*L, (A_measp2-A_sat2)/(A0-A_sat2), 'LineWidth', 4, 'Color', cmap(5,:))
semilogx(10^(-3)*L, (A_mal3-A_sat3)/(A0-A_sat3), 'LineWidth', 2)
semilogx(10^(-3)*L, (A_measp3-A_sat3)/(A0-A_sat3), 'LineWidth', 2, 'Color', cmap(6,:))
errorbar(maltose_dose_response(:,1), maltose_dose_response(:,2), maltose_dose_response(:,3), ...
    'o', 'MarkerSize', 10, 'Color', 'black', 'MarkerFaceColor', cmap(4,:), 'MarkerEdgeColor', 'black')
errorbar(measp_dose_response(:,1), measp_dose_response(:,2), measp_dose_response(:,3), ...
    'o', 'MarkerSize', 10, 'Color', 'black', 'MarkerFaceColor', cmap(7,:), 'MarkerEdgeColor', 'black')
ylim([0 1])
xlabel('Attractant (mM)')
ylabel('Relative kinase activity')
grid on
set(findall(gcf,'-property','FontSize'),'FontSize',16)
legend(['indirect binding model: maltose'], ['indirect binding model: MeAsp'], ...
    ['TS model, Fit A: maltose'], ['TS model, Fit A: MeAsp'], ['TS model, Fit B: maltose'], ['TS model, Fit B: MeAsp'], ...
    ['TS model, Fit C: maltose'], ['TS model, Fit C: MeAsp'],['maltose'], ['MeAsp'], 'FontSize', 14)
hold off

figure()
semilogx(10^(-3)*L, (A0-A_mal_3_N)/(A0-A_sat_N), ':', 'LineWidth', 8, 'Color', 'black')
hold on
semilogx(10^(-3)*L, (A0-A_measp_3_N)/(A0-A_sat_N), ':', 'LineWidth', 8, 'Color', cmap(8,:))
semilogx(10^(-3)*L, (A0-A_mal_3)/(A0-A_sat), 'LineWidth', 6, 'Color', cmap(2,:))
semilogx(10^(-3)*L, (A0-A_measp_3)/(A0-A_sat), 'LineWidth', 6, 'Color', cmap(1,:))
semilogx(10^(-3)*L, (A0-A_mal_3_2)/(A0-A_sat2), 'LineWidth', 4, 'Color', cmap(3,:))
semilogx(10^(-3)*L, (A0-A_measp_3_2)/(A0-A_sat2), 'LineWidth', 4, 'Color', cmap(5,:))
semilogx(10^(-3)*L, (A0-A_mal_3_3)/(A0-A_sat3), 'LineWidth', 2)
semilogx(10^(-3)*L, (A0-A_measp_3_3)/(A0-A_sat3), 'LineWidth', 2, 'Color', cmap(6,:))
errorbar(maltose_range(:,1), maltose_range(:,2), maltose_range(:,3), ...
    'o', 'MarkerSize', 10, 'Color', 'black', 'MarkerFaceColor', cmap(4,:), 'MarkerEdgeColor', 'black')
errorbar(measp_range(:,1), measp_range(:,2), measp_range(:,3), ...
    'o', 'MarkerSize', 10, 'Color', 'black', 'MarkerFaceColor', cmap(7,:), 'MarkerEdgeColor', 'black')
ylim([0 1])
xlabel('Attractant (mM)')
ylabel('Relative response')
grid on
set(findall(gcf,'-property','FontSize'),'FontSize',16)
hold off


%%%%%%%%%

figure()
semilogx(10^(-3)*L, (A_mal-A_sat)/(A0-A_sat), 'LineWidth', 6, 'Color', cmap(2,:))
hold on
semilogx(10^(-3)*L, (A_mal3-A_sat3)/(A0-A_sat3), 'LineWidth', 4, 'Color', 'black')
ylim([0 1])
xlabel('Attractant (mM)')
ylabel('Relative kinase activity')
grid on
set(findall(gcf,'-property','FontSize'),'FontSize',16)
legend(['TS model, Fit A: maltose'], ['TS model, Fit A: maltose, [BP] \times 4'], 'FontSize', 14)
hold off

figure()
semilogx(10^(-3)*L, (A0-A_mal_3)/(A0-A_sat), 'LineWidth', 6, 'Color', cmap(2,:))
hold on
semilogx(10^(-3)*L, (A0-A_mal_3_3)/(A0-A_sat3), 'LineWidth', 4, 'Color', 'black')
ylim([0 1])
xlabel('Attractant (mM)')
ylabel('Relative response')
grid on
set(findall(gcf,'-property','FontSize'),'FontSize',16)
hold off












