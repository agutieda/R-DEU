% This file:
% - Find intervals at which each subject changes from risky to safe lottery in
%   each task
% - Compute risk aversion revealed by each task as mid-point of switching points
%   in the interval
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

% Select observations and menus corresponding to risk only
menuTab = menuTab(strcmp(menuTab.taskType,'Risk'),:);
obsTab = obsTab(strcmp(obsTab.taskType,'Risk'),:);

%%% Level of background consumption
% recall that payoffs in the menuTab are expressed as thousands of DKK
omega = 118/1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) Compute indifference thresholds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For each menu, find the coefficient of risk aversion that makes a DM with
% CRRA utility indifferent between the two alternatives in a menu
menuTab.K_Ar = compute_thresholds_risk(menuTab, omega);
obsTab.K_Ar = menuTab.K_Ar(obsTab.menuID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2) Elicit implied coefficient of risk aversion by each subject-task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Subject list
subjectList = unique(obsTab.subjectID);
nSubjects = length(subjectList);

% Preallocate the output
raTab = table();

for subjectIdx = 1:nSubjects

    subject = subjectList(subjectIdx);

    % Create table for this subject
    subjectTab = obsTab(obsTab.subjectID == subject, :);
    taskList = unique(subjectTab.taskID);
    nTasks = length(taskList);

    % Iterate across tasks
    for taskIdx = 1:nTasks

        task = taskList(taskIdx);

        % Choose menus in this task
        taskTab = subjectTab(subjectTab.taskID == task, :);

        % Drop observations where indifference is reported
        taskTab(taskTab.Y == 0, :) = [];

        % Elicit the level of "r" under CRRA consistent with the observed choice
        if all(taskTab.Y == 1)

            % If DM always chooses safe lottery in this task, impute highest
            % level of risk aversion in the task list plus an epsilon equal
            % to the average difference between two consecutive levels of r
            ra_elicited = max(taskTab.K_Ar) + 0.38;
            noSwitch = 1;

        elseif all(taskTab.Y == 2)

            % If DM always chooses risky lottery in this task, impute lowest
            % level of risk aversion in the task list minus an epsilon equal
            % to the average difference between two consecutive levels of r
            ra_elicited = min(taskTab.K_Ar) - 0.38;
            noSwitch = 2;

        else

            % In this case, there is at least one switching point
            noSwitch = 0;

            % Find position of switching points
            taskTab.diffY = [0; diff(taskTab.Y)];
            idxSwitch = find(taskTab.diffY == 1);

            % Impute r for each switching point
            nSwitch = length(idxSwitch);
            ra_elicited = zeros(nSwitch, 1);
            for iR = 1:nSwitch
                r_u = taskTab.K_Ar(idxSwitch(iR));
                r_l = taskTab.K_Ar(idxSwitch(iR) - 1);
                ra_elicited(iR) = (r_u + r_l) / 2;
            end
        end

        % Add information to the list
        switchID = 0;
        for iR = 1:length(ra_elicited)
            switchID = switchID + 1;
            newDataRow = table( ...
                subject, task, switchID, ra_elicited(iR), noSwitch, ...
                'VariableNames', ...
                {'subjectID', 'taskSetID', 'switchID', 'ra', 'flag'});
            raTab = [raTab; newDataRow];
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2) Create table with elicited values and summary statistics by subject
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

raTabBySubj = table();

for iSubject = 1:nSubjects

    subjectID  = subjectList(iSubject);
    subjectTab = raTab(raTab.subjectID==subjectID,:);
    avg_ra     = mean(subjectTab.ra);
    std_ra     = std(subjectTab.ra);
    nAllSafe   = sum(subjectTab.flag == 1);
    nAllRisk   = sum(subjectTab.flag == 2);
    nNoSwitch  = sum(subjectTab.flag == 1 | subjectTab.flag == 2);

    newDataRow = table( ...
        iSubject, subjectID, avg_ra, std_ra, ...
        nAllSafe, nAllRisk, nNoSwitch, ...
        'VariableNames', ...
        {'subjectNum', 'subjectID', 'avg_ra', 'std_ra', ...
        'nAllSafe', 'nAllRisk', 'nNoSwitch'} );

    raTabBySubj = [raTabBySubj; newDataRow];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4) Create table with SPE at population level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu_ra  = mean(raTab.ra);
sig_ra = std(raTab.ra);
raTab_EST = table(mu_ra, sig_ra, 'VariableNames', {'mu_ra', 'sig_ra'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3) Export
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

writetable(raTab_EST,'./output/spe_pooled_ra.csv', ...
    'WriteRowNames',true);

writetable(raTab,'./output/raTab.csv', ...
    'WriteRowNames',true);

writetable(raTabBySubj,'./output/spe_subject_ra.csv', ...
    'WriteRowNames',true);
