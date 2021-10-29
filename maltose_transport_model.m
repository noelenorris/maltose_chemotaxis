
%% 
%
%  Porin-limited maltose transport 
%
%       This code generates plots using model of maltose transport,
%       demonstrating effects of porin-limited transport on the
%       concentrations of free and bound maltose in the periplasm.
%
%  Noele Norris
% 
%%



clear all
clc

syms Vc Kc Kbp BP Vp Kp Lext Lp

Lp_f = simplify(solve(Vc*(BP/(Kc+BP))*Lp/(Kc*Kbp/(Kc+BP)+Lp) == Vp*(Lext-Lp)/(Kp+Lext+Lp), Lp))

L_BP_f = simplify(BP*Lp_f(2)/(Kbp+Lp_f(2)))



%%
BP_v = 1000;
Kc_v = 100;
Kbp_v = 2;
Kp_v = 10^4;
Vp_v = 1;
Vc_v = Vp_v*0.00075;
% 
% BP_v = 1000;
% Kc_v = 100;
% Kbp_v = 2;
% Kp_v = 10^4;
% Vp_v = 1;
% Vc_v = Vp_v*1e-4;

Lext_it = 10.^[-2:.01:2];
%Lext_it = 1;

L_BP = [];
L_BP_approx = [];
L_P = [];
for i=1:1:length(Lext_it)
    Lext_v = Lext_it(i);
    L_P(i) = subs(Lp_f(end), {BP, Kc, Kbp, Kp, Vp, Vc, Lext}, {BP_v, Kc_v, Kbp_v, Kp_v, Vp_v, Vc_v, Lext_v});
    L_BP(i) = subs(L_BP_f(end), {BP, Kc, Kbp, Kp, Vp, Vc, Lext}, {BP_v, Kc_v, Kbp_v, Kp_v, Vp_v, Vc_v, Lext_v});
    L_BP_approx(i) = (Kc_v*Vp_v)*Lext_v/(Kp_v*Vc_v);
end

alpha = (Kc_v*Vp_v)/(Kp_v*Vc_v)

% Vc_v = Vp_v*1e-5;
% 
% 
% L_P2 = [];
% for i=1:1:length(Lext_it)
%     Lext_v = Lext_it(i);
%     L_P2(i) = subs(Lp_f(end), {BP, Kc, Kbp, Kp, Vp, Vc, Lext}, {BP_v, Kc_v, Kbp_v, Kp_v, Vp_v, Vc_v, Lext_v});
% end
% 
% 
% Vc_v = Vp_v*1e-3;
% 
% 
% L_P3 = [];
% for i=1:1:length(Lext_it)
%     Lext_v = Lext_it(i);
%     L_P3(i) = subs(Lp_f(end), {BP, Kc, Kbp, Kp, Vp, Vc, Lext}, {BP_v, Kc_v, Kbp_v, Kp_v, Vp_v, Vc_v, Lext_v});
% end

L_BP_np = BP_v*Lext_it./(Lext_it + Kbp_v);

%%



figure;
loglog(Lext_it, L_P,'LineWidth',5)
hold on
vline1 = xline(0.057)
vline1.LineWidth = 1;
set(vline1,'LineStyle','--')
set(vline1, 'Color', 'black')
vline2 = xline(14.3)
vline2.LineWidth = 1;
set(vline2,'LineStyle','--')
set(vline2, 'Color', 'black')
loglog(Lext_it, Lext_it, 'LineStyle', '-.', 'Color', 'black', 'LineWidth', 1)
xlabel('[Mal]_{ext} (\muM)')
ylabel('[Mal]_p (\muM)')
set(findall(gca,'-property','FontSize'),'FontSize',18)
grid on
hold off

figure;
loglog(Lext_it, L_BP,'LineWidth',5)
hold on
loglog(Lext_it, L_BP_approx,'--','LineWidth',5)
loglog(Lext_it, Lext_it, 'LineStyle', '-.', 'Color', 'black', 'LineWidth', 1)
loglog(Lext_it, L_BP_np, 'LineStyle', '-.', 'Color', 'black', 'LineWidth', 1)
vline1 = xline(0.057)
vline1.LineWidth = 1;
set(vline1,'LineStyle','--')
set(vline1, 'Color', 'black')
vline2 = xline(14.3)
vline2.LineWidth = 1;
set(vline2,'LineStyle','--')
set(vline2, 'Color', 'black')
xlabel('[Mal]_{ext} (\muM)')
ylabel('[Mal:BP] (\muM)')
set(findall(gca,'-property','FontSize'),'FontSize',18)
ylim([10^(-1) 10^4])
grid on
hold off
