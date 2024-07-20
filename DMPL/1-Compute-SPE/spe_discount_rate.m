% This file:
% - Find intervals at which each subject changes from early payoff to delayed
%   payoff in each task
% - Compute discount rate revealed by each task as mid-point of switching points
%   in the interval, given each draw of risk aversion found in the risk
%   eliciation tasks
% - Compute mean and standard deviation of risk aversion across tasks
% - Use baseline background consumption of 118 DKK

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0) Initialize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear;
clc;

% Tables with menu information and choices
menuTab = readtable('./input/menuTab.csv');
obsTab = readtable('./input/obsTab.csv');
raTab = readtable('./output/raTab.csv');

% Select observations and menus corresponding to risk only
menuTab = menuTab(strcmp(menuTab.taskType,'Time'),:);
obsTab = obsTab(strcmp(obsTab.taskType,'Time'),:);

%%% Level of background consumption
% recall that payoffs in the menuTab are expressed as thousands of DKK
omega = 118/1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) Elicit implied coefficient of risk aversion by each subject-task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Subject list
subjectList = unique(obsTab.subjectID);
nSubjects = length(subjectList);

% Preallocate the output
drTab = table();

for subjectIdx = 1:nSubjects

    subject = subjectList(subjectIdx);

    % Create table for this subject
    subjectTab = obsTab(obsTab.subjectID == subject, :);
    taskList = unique(subjectTab.taskID);
    nTasks = length(taskList);
    raList = raTab.ra(raTab.subjectID == subject);

    for iR = 1:length(raList)

        % For each menu, find the discount rate that makes a risk-averse DM with
        % CRRA utility and exponential discounting indifferent between the payoffs

        ra_i = raList(iR);
        K_At_r = compute_thresholds_time(menuTab, omega, ra_i);
        subjectTab.K_At_r = K_At_r(subjectTab.menuID - 40);

        % Iterate across tasks
        for taskIdx = 1:nTasks

            task = taskList(taskIdx);

            % Choose menus in this task
            taskTab = subjectTab(subjectTab.taskID == task, :);

            % Drop observations where indifference is reported
            taskTab(taskTab.Y == 0, :) = [];

            % Elicit the level of "dr" under exponential discounting and risk
            % neutrality consistent with the observed choice
            if all(taskTab.Y == 1)

                % If DM always chooses early payoff in this task, impute highest
                % discount rate in the task list
                dr_elicited = max(taskTab.K_At_r);
                noSwitch = 1;

            elseif all(taskTab.Y == 2)

                % If DM always chooses delayed payoff in this task, impute the
                % discount rate in the task list
                dr_elicited = min(taskTab.K_At_r);
                noSwitch = 2;

            else

                % In this case, there is at least one switching point
                noSwitch = 0;

                % Find position of switching points
                taskTab.diffY = [0; diff(taskTab.Y)];
                idxSwitch = find(taskTab.diffY == 1);

                % Impute r for each switching point
                nSwitch = length(idxSwitch);
                dr_elicited = zeros(nSwitch, 1);
                for iD = 1:nSwitch
                    dr_u = taskTab.K_At_r(idxSwitch(iD));
                    dr_l = taskTab.K_At_r(idxSwitch(iD) - 1);
                    dr_elicited(iD) = (dr_u + dr_l) / 2;
                end

            end

            % Add information to the list
            switchID = 0;
            for iD = 1:length(dr_elicited)
                switchID = switchID + 1;
                newDataRow = table( ...
                    subject, task, switchID, dr_elicited(iD), ra_i, noSwitch, ...
                    'VariableNames', ...
                    {'subjectID', 'taskSetID', 'switchID', 'dr', 'ra', 'flag'});
                drTab = [drTab; newDataRow];
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2) Create table with SPE by subject
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

drTabBySubj = table();

for iSubject = 1:nSubjects

    subjectID  = subjectList(iSubject);
    subjectTab = drTab(drTab.subjectID==subjectID,:);
    avg_dr     = mean(subjectTab.dr);
    std_dr     = std(subjectTab.dr);
    nAllEarly  = sum(subjectTab.flag == 1) / 4;
    nAllLater  = sum(subjectTab.flag == 2) / 4;
    nNoSwitch  = sum(subjectTab.flag == 1 | subjectTab.flag == 2) / 4;

    newDataRow = table( ...
        iSubject, subjectID, avg_dr, std_dr, ...
        nAllEarly, nAllLater, nNoSwitch, ...
        'VariableNames', ...
        {'subjectNum', 'subjectID', 'avg_dr', 'std_dr', ...
        'nAllEarly', 'nAllLater', 'nNoSwitch'} );

    drTabBySubj = [drTabBySubj; newDataRow];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3) Create table with SPE at population level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu_dr  = mean(drTab.dr);
sig_dr = std(drTab.dr);
rho    = corr(drTab.ra, drTab.dr);
drTab_EST  = table(mu_dr, sig_dr, rho, 'VariableNames', {'mu_dr', 'sig_dr', 'rho'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4) Export
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

writetable(drTab_EST,'./output/spe_pooled_dr.csv', ...
    'WriteRowNames',true);

writetable(drTabBySubj,'./output/spe_subject_dr.csv', ...
    'WriteRowNames',true);


