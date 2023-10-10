% This script:
% - Predicted and observed choice distribution under RDEU
% - By risk task

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Predicted and Observed Choice Distribution by task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Tables with menu information
menuTab = readtable('./input/menuTab.csv');

% Tables with estimated probability of choice in each menu
rhoY_hat = readtable('./input/rdeu_rhoY_hat.csv');

% Tables with observed probability of choice in each menu
rhoY_obs = readtable('./input/rdeu_rhoY_obs.csv');

% Auxiliary objects
yList = (0:0.05:1)';
nY = length(yList);
nM = height(menuTab);


% Compute predicted and observed choice distribution by task set
predictedDist_Y  = cell(6,1);
observedDist_Y = cell(6,1);
for jTask = 1:6

    rhoY_hat_j  = rhoY_hat{menuTab.taskID == jTask,:} ;
    rhoY_obs_j  = rhoY_obs{menuTab.taskID == jTask,:} ;
    nM_j = sum(menuTab.taskID == jTask);

    % Predicted probability of choosing each Y unconditionally
    predictedDist_Y_j = sum( rhoY_hat_j./nM_j , 1 );

    % Observed probability of choosing each Y unconditionally
    observedDist_Y_j = sum(rhoY_obs_j./nM_j , 1 );

    % Store
    predictedDist_Y{jTask} = predictedDist_Y_j;
    observedDist_Y{jTask} = observedDist_Y_j;

end


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

for jTask = 1:6

    fig = figure;

    bar(yList, observedDist_Y{jTask}, ...
        'LineWidth', line_width, ...
        'FaceColor', bar_color_obs, ...
        'EdgeColor', edge_color_obs, ...
        'FaceAlpha', face_alpha_plot, ...
        'EdgeAlpha',edge_alpha_plot );

    hold on;

    bar(yList, predictedDist_Y{jTask}, ...
        'LineWidth', line_width, ...
        'FaceColor', bar_color_hat, ...
        'EdgeColor', edge_color_hat, ...
        'FaceAlpha', face_alpha_plot, ...
        'EdgeAlpha',edge_alpha_plot );

    hold off;

    legend({'Observed', 'Predicted'},'Location','northwest','FontSize',font_text);
    ylim([0,0.8]);
    grid('on');

    xlabel('Share of Tokens $a$','FontSize',font_text,'Interpreter','latex');
    ylabel('Frequency','FontSize',font_text);

    dir_route = ['./output/cb_rdeu_predicted_choice_task_',num2str(jTask)];
    print(fig,dir_route,'-depsc2');


end