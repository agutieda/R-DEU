% This script:
% - Generate scatterplots of estimates of sigma_r and sigma_delta by subject

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0) Initialize and Import Estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot options
fill_color_main = "#000000";
edge_color_main = "#000000";
fill_color_alt1 = "#FFFFFF";
edge_color_alt1 = "#000000";
fill_color_alt2 = "#FFFFFF";
edge_color_alt2 = "#000000";
line_color = [0.5, 0.5, 0.5];
font_text = 14;

marker_face_alpha = 0.20;
marker_edge_alpha = 1.00;
marker_face_alpha_alt = 0.90;
marker_edge_alpha_alt = 1.00;



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
scatter(sp_ra(~selectedSample), rdeu_ra(~selectedSample), 50,...
    'MarkerFaceColor',fill_color_main,...
    'MarkerEdgeColor',edge_color_main,...
    'Marker','o',...    
    'MarkerFaceAlpha', marker_face_alpha, ...
    'MarkerEdgeAlpha', marker_edge_alpha, ...
    'LineWidth',1.0);
hold on;
scatter(sp_ra(selectedSample), rdeu_ra(selectedSample), 50,...
    'MarkerFaceColor',fill_color_alt1,...
    'MarkerEdgeColor',edge_color_alt1,...
    'Marker','^',...    
    'MarkerFaceAlpha', marker_face_alpha_alt, ...
    'MarkerEdgeAlpha', marker_edge_alpha_alt, ...
    'LineWidth',1.0);
xlabel('SPE','Interpreter','latex', 'FontSize', font_text)
ylabel('R-DEU','Interpreter','latex', 'FontSize', font_text)
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
scatter( sp_dr(~selectedSample), rdeu_dr(~selectedSample), 50,...
    'MarkerFaceColor',fill_color_main,...
    'MarkerEdgeColor',edge_color_main,...
    'Marker','o',...    
    'MarkerFaceAlpha', marker_face_alpha, ...
    'MarkerEdgeAlpha', marker_edge_alpha, ...
    'LineWidth',1.0);
hold on;
scatter(sp_dr(selectedSample), rdeu_dr(selectedSample), 60,...
    'MarkerFaceColor',fill_color_alt2,...
    'MarkerEdgeColor',edge_color_alt2,...
    'Marker','square',...    
    'MarkerFaceAlpha', marker_face_alpha_alt, ...
    'MarkerEdgeAlpha', marker_edge_alpha_alt, ...
    'LineWidth',1.0);
xlabel('SPE','Interpreter','latex', 'FontSize', font_text)
ylabel('R-DEU','Interpreter','latex', 'FontSize', font_text)
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
scatter(sp_ra(~selectedSample), rumL_ra(~selectedSample), 50,...
    'MarkerFaceColor',fill_color_main,...
    'MarkerEdgeColor',edge_color_main,...
    'Marker','o',...    
    'MarkerFaceAlpha', marker_face_alpha, ...
    'MarkerEdgeAlpha', marker_edge_alpha, ...
    'LineWidth',1.0);
hold on;
scatter(sp_ra(selectedSample), rumL_ra(selectedSample), 50,...
    'MarkerFaceColor',fill_color_alt1,...
    'MarkerEdgeColor',edge_color_alt1,...
    'Marker','^',...    
    'MarkerFaceAlpha', marker_face_alpha_alt, ...
    'MarkerEdgeAlpha', marker_edge_alpha_alt, ...
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
scatter( sp_dr(~selectedSample), rumL_dr(~selectedSample), 50,...
    'MarkerFaceColor',fill_color_main,...
    'MarkerEdgeColor',edge_color_main,...
    'Marker','o',...    
    'MarkerFaceAlpha', marker_face_alpha, ...
    'MarkerEdgeAlpha', marker_edge_alpha, ...
    'LineWidth',1.0);
hold on;
scatter(sp_dr(selectedSample), rumL_dr(selectedSample), 60,...
    'MarkerFaceColor',fill_color_alt2,...
    'MarkerEdgeColor',edge_color_alt2,...
    'Marker','square',...    
    'MarkerFaceAlpha', marker_face_alpha_alt, ...
    'MarkerEdgeAlpha', marker_edge_alpha_alt, ...
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
scatter(sp_ra(~selectedSample), rumW_ra(~selectedSample), 50,...
    'MarkerFaceColor',fill_color_main,...
    'MarkerEdgeColor',edge_color_main,...
    'Marker','o',...    
    'MarkerFaceAlpha', marker_face_alpha, ...
    'MarkerEdgeAlpha', marker_edge_alpha, ...
    'LineWidth',1.0);
hold on;
scatter(sp_ra(selectedSample), rumW_ra(selectedSample), 50,...
    'MarkerFaceColor',fill_color_alt1,...
    'MarkerEdgeColor',edge_color_alt1,...
    'Marker','^',...    
    'MarkerFaceAlpha', marker_face_alpha_alt, ...
    'MarkerEdgeAlpha', marker_edge_alpha_alt, ...
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
scatter( sp_dr(~selectedSample), rumW_dr(~selectedSample), 50,...
    'MarkerFaceColor',fill_color_main,...
    'MarkerEdgeColor',edge_color_main,...
    'Marker','o',...    
    'MarkerFaceAlpha', marker_face_alpha, ...
    'MarkerEdgeAlpha', marker_edge_alpha, ...
    'LineWidth',1.0);
hold on;
scatter(sp_dr(selectedSample), rumW_dr(selectedSample), 60,...
    'MarkerFaceColor',fill_color_alt2,...
    'MarkerEdgeColor',edge_color_alt2,...
    'Marker','square',...    
    'MarkerFaceAlpha', marker_face_alpha_alt, ...
    'MarkerEdgeAlpha', marker_edge_alpha_alt, ...
    'LineWidth',1.0);
xlabel('SPE','Interpreter','latex', 'FontSize', font_text)
ylabel('WILCOX','Interpreter','latex', 'FontSize', font_text)
xlim([dr_LB_x,dr_UB_x]);
ylim([dr_LB_y,dr_UB_y]);
regLine = refline(1,0);
set(regLine, 'Parent', ax, 'LineWidth',1.0,'LineStyle','--','Color',line_color);
set(ax,'XGrid','on','YGrid','on');
print(fig,'./output/dmpl_scatter_rum_wilcox_sig_dr','-depsc2');
