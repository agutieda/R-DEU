% This script:
% - Predicted and observed choice distribution under RUM
% - Menu 71

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Predicted and Observed Choice Distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Tables with menu information
menuTab = readtable('./input/menuTab.csv');

% Tables with estimated probability of choice in each menu
rhoY_hat = readtable('./input/rum_rhoY_hat.csv');

% Tables with observed probability of choice in each menu
rhoY_obs = readtable('./input/rum_rhoY_obs.csv');

% Auxiliary objects
yList = (0:0.05:1)';
nY = length(yList);
nM = height(menuTab);

% chosen menu 
chosen_menu = 71;

% Predicted probability of choosing each Y unconditionally
predictedDist_Y = rhoY_hat{chosen_menu,:};

% Observed probability of choosing each Y unconditionally
observedDist_Y = rhoY_obs{chosen_menu,:};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot options
face_alpha_plot = 0.6;
edge_alpha_plot = 0.9;

bar_color_obs  = [0 0.4470 0.7410];
edge_color_obs = [0, 0, 0];

bar_color_hat = [0.8500, 0.3250, 0.0980];
edge_color_hat = [0, 0, 0];

line_color = [0.6350, 0.0780, 0.1840];
line_width = 1.5;
font_text = 14;

% Print individual
fig = figure;

bar(yList, observedDist_Y, ...
    'LineWidth', line_width, ...
    'FaceColor', bar_color_obs, ...
    'EdgeColor', edge_color_obs, ...
    'FaceAlpha', face_alpha_plot, ...
    'EdgeAlpha',edge_alpha_plot );

hold on;

bar(yList, predictedDist_Y, ...
    'LineWidth', line_width, ...
    'FaceColor', bar_color_hat, ...
    'EdgeColor', edge_color_hat, ...
    'FaceAlpha', face_alpha_plot, ...
    'EdgeAlpha',edge_alpha_plot );

hold off;

legend({'Observed', 'Predicted'},'Location','northwest','FontSize',font_text);
ylim([0,1.05]);
grid('on');

xlabel('Share of Tokens $a$','FontSize',font_text,'Interpreter','latex');
ylabel('Frequency','FontSize',font_text);

print(fig,'./output/cb_rum_predicted_choice_menu','-depsc2');
