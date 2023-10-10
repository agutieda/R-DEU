% This script:
% - Generate scatterplots of estimates of sigma_r and sigma_delta by subject

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0) Initialize and Import Estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rdeuTab = readtable('./input/rdeu_subject_mle.csv');
rumLTab = readtable('./input/luce_subject_mle.csv');
rumWTab = readtable('./input/wilcox_subject_mle.csv');
speRTab = readtable('./input/spe_subject_ra.csv','ReadRowNames',true);
speTTab = readtable('./input/spe_subject_dr.csv','ReadRowNames',true);

% Plot options
fill_color_main = "#0066E6";
edge_color_main = "#002E8A";
fill_color_alt1 = "#FF4D00";
edge_color_alt1 = "#D90000";
fill_color_alt2 = "#4D19FF";
edge_color_alt2 = "#3612B3";
line_color = [0.5, 0.5, 0.5];
font_text = 14;

% Limits
ra_LB_x  =  0.0;
ra_UB_x  =  2.5;
ra_LB_y  =  0.0;
ra_UB_y  =  2.5;

dr_LB_x  =  0.0;
dr_UB_x  =  0.7;
dr_LB_y  =  0.0;
dr_UB_y  =  0.7;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) Recover MLE and truncate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract variables
rdeu_ra = rdeuTab.sig_ra;
rdeu_dr = rdeuTab.sig_dr;
rumL_ra = rumLTab.sig_ra;
rumL_dr = rumLTab.sig_dr;
rumW_ra = rumWTab.sig_ra;
rumW_dr = rumWTab.sig_dr;
sp_ra   = spRTab.std_ra;
sp_dr   = spTTab.std_dr;
sp_ra_nS = spRTab.nNoSwitch;
sp_dr_nS = spTTab.nNoSwitch;

% Truncate
rdeu_ra = min(max(rdeu_ra,ra_LB_y),ra_UB_y);
rdeu_dr = min(max(rdeu_dr,dr_LB_y),dr_UB_y);
rumL_ra = min(max(rumL_ra,ra_LB_y),ra_UB_y);
rumL_dr = min(max(rumL_dr,dr_LB_y),dr_UB_y);
rumW_ra = min(max(rumW_ra,ra_LB_y),ra_UB_y);
rumW_dr = min(max(rumW_dr,dr_LB_y),dr_UB_y);
sp_ra   = min(max(sp_ra,ra_LB_x),ra_UB_x);
sp_dr   = min(max(sp_dr,dr_LB_x),dr_UB_x);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2) RDEU vs SP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Risk Aversion
selectedSample = sp_ra_nS>=1;

fig = figure;
ax = axes('Parent',fig);
hold(ax,'on');
scatter(sp_ra(~selectedSample), rdeu_ra(~selectedSample), 60,...
    'MarkerFaceColor',fill_color_main,...
    'MarkerEdgeColor',edge_color_main,...
    'Marker','o',...    
    'MarkerFaceAlpha', 0.6, ...
    'MarkerEdgeAlpha', 0.9, ...
    'LineWidth',1.0);
hold on;
scatter(sp_ra(selectedSample), rdeu_ra(selectedSample), 60,...
    'MarkerFaceColor',fill_color_alt1,...
    'MarkerEdgeColor',edge_color_alt1,...
    'Marker','o',...    
    'MarkerFaceAlpha', 0.6, ...
    'MarkerEdgeAlpha', 0.9, ...
    'LineWidth',1.0);
xlabel('SPE','Interpreter','latex', 'FontSize', font_text)
ylabel('RDEU','Interpreter','latex', 'FontSize', font_text)
xlim([ra_LB_x,ra_UB_x]);
ylim([ra_LB_y,ra_UB_y]);
regLine = refline(1,0);
set(regLine, 'Parent', ax, 'LineWidth',1.0,'LineStyle','--','Color',line_color);
set(ax,'XGrid','on','YGrid','on');
print(fig,'./output/dmpl_scatter_rdeu_sig_ra','-depsc2');

% Discount Rate
selectedSample = sp_dr_nS>=1;
fig = figure;
ax = axes('Parent',fig);
hold(ax,'on');
scatter( sp_dr(~selectedSample), rdeu_dr(~selectedSample), 60,...
    'MarkerFaceColor',fill_color_main,...
    'MarkerEdgeColor',edge_color_main,...
    'Marker','o',...    
    'MarkerFaceAlpha', 0.6, ...
    'MarkerEdgeAlpha', 0.9, ...
    'LineWidth',1.0);
hold on;
scatter(sp_dr(selectedSample), rdeu_dr(selectedSample), 60,...
    'MarkerFaceColor',fill_color_alt2,...
    'MarkerEdgeColor',edge_color_alt2,...
    'Marker','o',...    
    'MarkerFaceAlpha', 0.6, ...
    'MarkerEdgeAlpha', 0.9, ...
    'LineWidth',1.0);
xlabel('SPE','Interpreter','latex', 'FontSize', font_text)
ylabel('RDEU','Interpreter','latex', 'FontSize', font_text)
xlim([dr_LB_x,dr_UB_x]);
ylim([dr_LB_y,dr_UB_y]);
regLine = refline(1,0);
set(regLine, 'Parent', ax, 'LineWidth',1.0,'LineStyle','--','Color',line_color);
set(ax,'XGrid','on','YGrid','on');
print(fig,'./output/dmpl_scatter_rdeu_sig_dr','-depsc2');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3) RUM-LUCE vs SP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Risk Aversion
selectedSample = sp_ra_nS>=1;
fig = figure;
ax = axes('Parent',fig);
hold(ax,'on');
scatter(sp_ra(~selectedSample), rumL_ra(~selectedSample), 60,...
    'MarkerFaceColor',fill_color_main,...
    'MarkerEdgeColor',edge_color_main,...
    'Marker','o',...    
    'MarkerFaceAlpha', 0.6, ...
    'MarkerEdgeAlpha', 0.9, ...
    'LineWidth',1.0);
hold on;
scatter(sp_ra(selectedSample), rumL_ra(selectedSample), 60,...
    'MarkerFaceColor',fill_color_alt1,...
    'MarkerEdgeColor',edge_color_alt1,...
    'Marker','o',...    
    'MarkerFaceAlpha', 0.6, ...
    'MarkerEdgeAlpha', 0.9, ...
    'LineWidth',1.0);
xlabel('SPE','Interpreter','latex', 'FontSize', font_text)
ylabel('LUCE','Interpreter','latex', 'FontSize', font_text)
xlim([ra_LB_x,ra_UB_x]);
ylim([ra_LB_y,ra_UB_y]);
regLine = refline(1,0);
set(regLine, 'Parent', ax, 'LineWidth',1.0,'LineStyle','--','Color',line_color);
set(ax,'XGrid','on','YGrid','on');
print(fig,'./output/dmpl_scatter_rum_luce_sig_ra','-depsc2');

% Discount Rate
selectedSample = sp_dr_nS>=1;
fig = figure;
ax = axes('Parent',fig);
hold(ax,'on');
scatter( sp_dr(~selectedSample), rumL_dr(~selectedSample), 60,...
    'MarkerFaceColor',fill_color_main,...
    'MarkerEdgeColor',edge_color_main,...
    'Marker','o',...    
    'MarkerFaceAlpha', 0.6, ...
    'MarkerEdgeAlpha', 0.9, ...
    'LineWidth',1.0);
hold on;
scatter(sp_dr(selectedSample), rumL_dr(selectedSample), 60,...
    'MarkerFaceColor',fill_color_alt2,...
    'MarkerEdgeColor',edge_color_alt2,...
    'Marker','o',...    
    'MarkerFaceAlpha', 0.6, ...
    'MarkerEdgeAlpha', 0.9, ...
    'LineWidth',1.0);
xlabel('SPE','Interpreter','latex', 'FontSize', font_text)
ylabel('LUCE','Interpreter','latex', 'FontSize', font_text)
xlim([dr_LB_x,dr_UB_x]);
ylim([dr_LB_y,dr_UB_y]);
regLine = refline(1,0);
set(regLine, 'Parent', ax, 'LineWidth',1.0,'LineStyle','--','Color',line_color);
set(ax,'XGrid','on','YGrid','on');
print(fig,'./output/dmpl_scatter_rum_luce_sig_dr','-depsc2');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4) RUM-WILCOX vs SP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Risk Aversion
selectedSample = sp_ra_nS>=1;
fig = figure;
ax = axes('Parent',fig);
hold(ax,'on');
scatter(sp_ra(~selectedSample), rumW_ra(~selectedSample), 60,...
    'MarkerFaceColor',fill_color_main,...
    'MarkerEdgeColor',edge_color_main,...
    'Marker','o',...    
    'MarkerFaceAlpha', 0.6, ...
    'MarkerEdgeAlpha', 0.9, ...
    'LineWidth',1.0);
hold on;
scatter(sp_ra(selectedSample), rumW_ra(selectedSample), 60,...
    'MarkerFaceColor',fill_color_alt1,...
    'MarkerEdgeColor',edge_color_alt1,...
    'Marker','o',...    
    'MarkerFaceAlpha', 0.6, ...
    'MarkerEdgeAlpha', 0.9, ...
    'LineWidth',1.0);
xlabel('SPE','Interpreter','latex', 'FontSize', font_text)
ylabel('WILCOX','Interpreter','latex', 'FontSize', font_text)
xlim([ra_LB_x,ra_UB_x]);
ylim([ra_LB_y,ra_UB_y]);
regLine = refline(1,0);
set(regLine, 'Parent', ax, 'LineWidth',1.0,'LineStyle','--','Color',line_color);
set(ax,'XGrid','on','YGrid','on');
print(fig,'./output/dmpl_scatter_rum_wilcox_sig_ra','-depsc2');

% Discount Rate
selectedSample = sp_dr_nS>=1;
fig = figure;
ax = axes('Parent',fig);
hold(ax,'on');
scatter( sp_dr(~selectedSample), rumW_dr(~selectedSample), 60,...
    'MarkerFaceColor',fill_color_main,...
    'MarkerEdgeColor',edge_color_main,...
    'Marker','o',...    
    'MarkerFaceAlpha', 0.6, ...
    'MarkerEdgeAlpha', 0.9, ...
    'LineWidth',1.0);
hold on;
scatter(sp_dr(selectedSample), rumW_dr(selectedSample), 60,...
    'MarkerFaceColor',fill_color_alt2,...
    'MarkerEdgeColor',edge_color_alt2,...
    'Marker','o',...    
    'MarkerFaceAlpha', 0.6, ...
    'MarkerEdgeAlpha', 0.9, ...
    'LineWidth',1.0);
xlabel('SPE','Interpreter','latex', 'FontSize', font_text)
ylabel('WILCOX','Interpreter','latex', 'FontSize', font_text)
xlim([dr_LB_x,dr_UB_x]);
ylim([dr_LB_y,dr_UB_y]);
regLine = refline(1,0);
set(regLine, 'Parent', ax, 'LineWidth',1.0,'LineStyle','--','Color',line_color);
set(ax,'XGrid','on','YGrid','on');
print(fig,'./output/dmpl_scatter_rum_wilcox_sig_dr','-depsc2');
