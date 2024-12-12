% Created: October 16, 2024
% Last updated: December 11, 2024
% Kyra Veprek, kyraveprek24@gmail.com
%
% This code reads raw data from the csv files of Kyra Veprek's FYP for Experiment 5 and
% outputs a structure (dataSorted) which nicely organizes the trial data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;

N_SUBJECTS = 12; %Number of subjects
N_TRIALS = 122; %Number of Trials + Practice Trials
hZ = 75;

shift = -40;
theta = deg2rad(shift);
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];

% Read in the data from csv files 
M1 = readmatrix('filesSubj01.csv', 'Range', 1);
M2 = readmatrix('filesSubj02.csv', 'Range', 1);
M3 = readmatrix('filesSubj03.csv', 'Range', 1);
M4 = readmatrix('filesSubj04.csv', 'Range', 1);
M5 = readmatrix('filesSubj05.csv', 'Range', 1);
M6 = readmatrix('filesSubj06.csv', 'Range', 1);
M7 = readmatrix('filesSubj07.csv', 'Range', 1);
M8 = readmatrix('filesSubj08.csv', 'Range', 1);
M9 = readmatrix('filesSubj09.csv', 'Range', 1);
M10 = readmatrix('filesSubj10.csv', 'Range', 1);
M11 = readmatrix('filesSubj11.csv', 'Range', 1);
M12 = readmatrix('filesSubj12.csv', 'Range', 1);

%%
% Remove unnecessary timesteps before the obstacle has appeared
obstacleColumn = 9;
M1_new = remove_rampup_rows(M1, obstacleColumn);
M2_new = remove_rampup_rows(M2, obstacleColumn);
M3_new = remove_rampup_rows(M3, obstacleColumn);
M4_new = remove_rampup_rows(M4, obstacleColumn);
M5_new = remove_rampup_rows(M5, obstacleColumn);
M6_new = remove_rampup_rows(M6, obstacleColumn);
M7_new = remove_rampup_rows(M7, obstacleColumn);
M8_new = remove_rampup_rows(M8, obstacleColumn);
M9_new = remove_rampup_rows(M9, obstacleColumn);
M10_new = remove_rampup_rows(M10, obstacleColumn);
M11_new = remove_rampup_rows(M11, obstacleColumn);
M12_new = remove_rampup_rows(M12, obstacleColumn);

%Determine where each trial begins and ends by identifying the cell locations of the 'NaNs'
pt1 = find(isnan(M1_new(3:size(M1_new),1)));
pt2 = find(isnan(M2_new(3:size(M2_new),1)));
pt3 = find(isnan(M3_new(3:size(M3_new),1)));
pt4 = find(isnan(M4_new(3:size(M4_new),1)));
pt5 = find(isnan(M5_new(3:size(M5_new),1)));
pt6 = find(isnan(M6_new(3:size(M6_new),1)));
pt7 = find(isnan(M7_new(3:size(M7_new),1)));
pt8 = find(isnan(M8_new(3:size(M8_new),1)));
pt9 = find(isnan(M9_new(3:size(M9_new),1)));
pt10 = find(isnan(M10_new(3:size(M10_new),1)));
pt11 = find(isnan(M11_new(3:size(M11_new),1)));
pt12 = find(isnan(M12_new(3:size(M12_new),1)));

%% Subject 1

%dataSorted(N_SUBJECTS,N_TRIALS) = struct;

for iS = 1
    for iT = 1:N_TRIALS
        [start_row, end_row, info_row] = get_trial_rows(iT, pt1, height(M1_new));
        dataSorted(iS,iT).info.TNum = M1_new(info_row,2);
        dataSorted(iS,iT).info.Ang = M1_new(info_row,4);
        dataSorted(iS,iT).info.Speed = M1_new(info_row,6);
        dataSorted(iS,iT).info.Size = M1_new(info_row,8);
        dataSorted(iS,iT).Time = M1_new(start_row:end_row, 1);
        dataSorted(iS,iT).SubjectX = M1_new(start_row:end_row,3);
        dataSorted(iS,iT).SubjectZ = M1_new(start_row:end_row,5);
        dataSorted(iS,iT).SubjectYaw = M1_new(start_row:end_row,7);
        dataSorted(iS,iT).SubjectPitch = M1_new(start_row:end_row,6);
        dataSorted(iS,iT).SubjectRoll = M1_new(start_row:end_row,8);
        dataSorted(iS,iT).obstX = M1_new(start_row:end_row,9);
        dataSorted(iS,iT).obstY = M1_new(start_row:end_row,10);
        dataSorted(iS,iT).GoalX = M1_new(start_row:end_row,11);
        dataSorted(iS,iT).GoalY = M1_new(start_row:end_row,12);
        if dataSorted(iS,iT).info.Ang ~= 0
            dataSorted(iS,iT).Xfilt = filter_butter(hZ,dataSorted(iS,iT).SubjectX);
            dataSorted(iS,iT).Zfilt = filter_butter(hZ,dataSorted(iS,iT).SubjectZ);
            dataSorted(iS,iT).P(:,1) = dataSorted(iS,iT).Xfilt;
            dataSorted(iS,iT).P(:,2) = dataSorted(iS,iT).Zfilt;
            dataSorted(iS,iT).A1P(:,1) = dataSorted(iS,iT).obstX;
            dataSorted(iS,iT).A1P(:,2) = dataSorted(iS,iT).obstY;
            dataSorted(iS, iT).distance = pdist2(dataSorted(iS, iT).P, dataSorted(iS, iT).A1P, 'euclidean');
            for i = 1:length(dataSorted(iS,iT).P)
               dataSorted(iS,iT).distance(1,i) = pdist2(dataSorted(iS,iT).P(i,:), dataSorted(iS,iT).A1P(i,:), 'euclidean');
               dataSorted(iS,iT).theta(1,i) = 2*atan(dataSorted(iS,iT).info.Size/dataSorted(iS,iT).distance(1,i));
            end
        end
        fprintf('Completed subject: %d, trial: %d\n', iS, iT);
    end
end

%%
for iS = 2
    for iT = 1:N_TRIALS
        [start_row, end_row, info_row] = get_trial_rows(iT, pt2, height(M2_new));
        dataSorted(iS,iT).info.TNum = M2_new(info_row,2);
        dataSorted(iS,iT).info.Ang = M2_new(info_row,4);
        dataSorted(iS,iT).info.Speed = M2_new(info_row,6);
        dataSorted(iS,iT).info.Size = M2_new(info_row,8);
        dataSorted(iS,iT).Time = M2_new(start_row:end_row, 1);
        dataSorted(iS,iT).SubjectX = M2_new(start_row:end_row,3);
        dataSorted(iS,iT).SubjectZ = M2_new(start_row:end_row,5);
        dataSorted(iS,iT).SubjectYaw = M2_new(start_row:end_row,7);
        dataSorted(iS,iT).SubjectPitch = M2_new(start_row:end_row,6);
        dataSorted(iS,iT).SubjectRoll = M2_new(start_row:end_row,8);
        dataSorted(iS,iT).obstX = M2_new(start_row:end_row,9);
        dataSorted(iS,iT).obstY = M2_new(start_row:end_row,10);
        dataSorted(iS,iT).GoalX = M2_new(start_row:end_row,11);
        dataSorted(iS,iT).GoalY = M2_new(start_row:end_row,12);
        if dataSorted(iS,iT).info.Ang ~= 0
            dataSorted(iS,iT).Xfilt = filter_butter(hZ,dataSorted(iS,iT).SubjectX);
            dataSorted(iS,iT).Zfilt = filter_butter(hZ,dataSorted(iS,iT).SubjectZ);
            dataSorted(iS,iT).P(:,1) = dataSorted(iS,iT).Xfilt;
            dataSorted(iS,iT).P(:,2) = dataSorted(iS,iT).Zfilt;
            dataSorted(iS,iT).A1P(:,1) = dataSorted(iS,iT).obstX;
            dataSorted(iS,iT).A1P(:,2) = dataSorted(iS,iT).obstY;
            for i = 1:length(dataSorted(iS,iT).P)
                dataSorted(iS,iT).distance(1,i) = pdist2(dataSorted(iS,iT).P(i,:), dataSorted(iS,iT).A1P(i,:), 'euclidean');
                dataSorted(iS,iT).theta(1,i) = 2*atan(dataSorted(iS,iT).info.Size/dataSorted(iS,iT).distance(1,i));
            end
        end
        fprintf('Completed subject: %d, trial: %d\n', iS, iT);
    end
end

%%
for iS = 3
    for iT = 1:N_TRIALS
        [start_row, end_row, info_row] = get_trial_rows(iT, pt3, height(M3_new));
        dataSorted(iS,iT).info.TNum = M3_new(info_row,2);
        dataSorted(iS,iT).info.Ang = M3_new(info_row,4);
        dataSorted(iS,iT).info.Speed = M3_new(info_row,6);
        dataSorted(iS,iT).info.Size = M3_new(info_row,8);
        dataSorted(iS,iT).Time = M3_new(start_row:end_row, 1);
        dataSorted(iS,iT).SubjectX = M3_new(start_row:end_row,3);
        dataSorted(iS,iT).SubjectZ = M3_new(start_row:end_row,5);
        dataSorted(iS,iT).SubjectYaw = M3_new(start_row:end_row,7);
        dataSorted(iS,iT).SubjectPitch = M3_new(start_row:end_row,6);
        dataSorted(iS,iT).SubjectRoll = M3_new(start_row:end_row,8);
        dataSorted(iS,iT).obstX = M3_new(start_row:end_row,9);
        dataSorted(iS,iT).obstY = M3_new(start_row:end_row,10);
        dataSorted(iS,iT).GoalX = M3_new(start_row:end_row,11);
        dataSorted(iS,iT).GoalY = M3_new(start_row:end_row,12);
        if dataSorted(iS,iT).info.Ang ~= 0
            dataSorted(iS,iT).Xfilt = filter_butter(hZ,dataSorted(iS,iT).SubjectX);
            dataSorted(iS,iT).Zfilt = filter_butter(hZ,dataSorted(iS,iT).SubjectZ);
            dataSorted(iS,iT).P(:,1) = dataSorted(iS,iT).Xfilt;
            dataSorted(iS,iT).P(:,2) = dataSorted(iS,iT).Zfilt;
            dataSorted(iS,iT).A1P(:,1) = dataSorted(iS,iT).obstX;
            dataSorted(iS,iT).A1P(:,2) = dataSorted(iS,iT).obstY;
            for i = 1:length(dataSorted(iS,iT).P)
                dataSorted(iS,iT).distance(1,i) = pdist2(dataSorted(iS,iT).P(i,:), dataSorted(iS,iT).A1P(i,:), 'euclidean');
                dataSorted(iS,iT).theta(1,i) = 2*atan(dataSorted(iS,iT).info.Size/dataSorted(iS,iT).distance(1,i));
            end
        end
    end
end


for iS = 4
    for iT = 1:N_TRIALS
        [start_row, end_row, info_row] = get_trial_rows(iT, pt4, height(M4_new));
        dataSorted(iS,iT).info.TNum = M4_new(info_row,2);
        dataSorted(iS,iT).info.Ang = M4_new(info_row,4);
        dataSorted(iS,iT).info.Speed = M4_new(info_row,6);
        dataSorted(iS,iT).info.Size = M4_new(info_row,8);
        dataSorted(iS,iT).Time = M4_new(start_row:end_row, 1);
        dataSorted(iS,iT).SubjectX = M4_new(start_row:end_row,3);
        dataSorted(iS,iT).SubjectZ = M4_new(start_row:end_row,5);
        dataSorted(iS,iT).SubjectYaw = M4_new(start_row:end_row,7);
        dataSorted(iS,iT).SubjectPitch = M4_new(start_row:end_row,6);
        dataSorted(iS,iT).SubjectRoll = M4_new(start_row:end_row,8);
        dataSorted(iS,iT).obstX = M4_new(start_row:end_row,9);
        dataSorted(iS,iT).obstY = M4_new(start_row:end_row,10);
        dataSorted(iS,iT).GoalX = M4_new(start_row:end_row,11);
        dataSorted(iS,iT).GoalY = M4_new(start_row:end_row,12);
        if dataSorted(iS,iT).info.Ang ~= 0
            dataSorted(iS,iT).Xfilt = filter_butter(hZ,dataSorted(iS,iT).SubjectX);
            dataSorted(iS,iT).Zfilt = filter_butter(hZ,dataSorted(iS,iT).SubjectZ);
            dataSorted(iS,iT).P(:,1) = dataSorted(iS,iT).Xfilt;
            dataSorted(iS,iT).P(:,2) = dataSorted(iS,iT).Zfilt;
            dataSorted(iS,iT).A1P(:,1) = dataSorted(iS,iT).obstX;
            dataSorted(iS,iT).A1P(:,2) = dataSorted(iS,iT).obstY;
            for i = 1:length(dataSorted(iS,iT).P)
                dataSorted(iS,iT).distance(1,i) = pdist2(dataSorted(iS,iT).P(i,:), dataSorted(iS,iT).A1P(i,:), 'euclidean');
                dataSorted(iS,iT).theta(1,i) = 2*atan(dataSorted(iS,iT).info.Size/dataSorted(iS,iT).distance(1,i));
            end
        end
        fprintf('Completed subject: %d, trial: %d\n', iS, iT);
    end
end

for iS = 5
    for iT = 1:N_TRIALS
        [start_row, end_row, info_row] = get_trial_rows(iT, pt5, height(M5_new));
        dataSorted(iS,iT).info.TNum = M5_new(info_row,2);
        dataSorted(iS,iT).info.Ang = M5_new(info_row,4);
        dataSorted(iS,iT).info.Speed = M5_new(info_row,6);
        dataSorted(iS,iT).info.Size = M5_new(info_row,8);
        dataSorted(iS,iT).Time = M5_new(start_row:end_row, 1);
        dataSorted(iS,iT).SubjectX = M5_new(start_row:end_row,3);
        dataSorted(iS,iT).SubjectZ = M5_new(start_row:end_row,5);
        dataSorted(iS,iT).SubjectYaw = M5_new(start_row:end_row,7);
        dataSorted(iS,iT).SubjectPitch = M5_new(start_row:end_row,6);
        dataSorted(iS,iT).SubjectRoll = M5_new(start_row:end_row,8);
        dataSorted(iS,iT).obstX = M5_new(start_row:end_row,9);
        dataSorted(iS,iT).obstY = M5_new(start_row:end_row,10);
        dataSorted(iS,iT).GoalX = M5_new(start_row:end_row,11);
        dataSorted(iS,iT).GoalY = M5_new(start_row:end_row,12);
        if dataSorted(iS,iT).info.Ang ~= 0
            dataSorted(iS,iT).Xfilt = filter_butter(hZ,dataSorted(iS,iT).SubjectX);
            dataSorted(iS,iT).Zfilt = filter_butter(hZ,dataSorted(iS,iT).SubjectZ);
            dataSorted(iS,iT).P(:,1) = dataSorted(iS,iT).Xfilt;
            dataSorted(iS,iT).P(:,2) = dataSorted(iS,iT).Zfilt;
            dataSorted(iS,iT).A1P(:,1) = dataSorted(iS,iT).obstX;
            dataSorted(iS,iT).A1P(:,2) = dataSorted(iS,iT).obstY;
            for i = 1:length(dataSorted(iS,iT).P)
                dataSorted(iS,iT).distance(1,i) = pdist2(dataSorted(iS,iT).P(i,:), dataSorted(iS,iT).A1P(i,:), 'euclidean');
                dataSorted(iS,iT).theta(1,i) = 2*atan(dataSorted(iS,iT).info.Size/dataSorted(iS,iT).distance(1,i));
            end
        end
        fprintf('Completed subject: %d, trial: %d\n', iS, iT);
    end
end

for iS = 6
    for iT = 1:N_TRIALS
        [start_row, end_row, info_row] = get_trial_rows(iT, pt6, height(M6_new));
        dataSorted(iS,iT).info.TNum = M6_new(info_row,2);
        dataSorted(iS,iT).info.Ang = M6_new(info_row,4);
        dataSorted(iS,iT).info.Speed = M6_new(info_row,6);
        dataSorted(iS,iT).info.Size = M6_new(info_row,8);
        dataSorted(iS,iT).Time = M6_new(start_row:end_row, 1);
        dataSorted(iS,iT).SubjectX = M6_new(start_row:end_row,3);
        dataSorted(iS,iT).SubjectZ = M6_new(start_row:end_row,5);
        dataSorted(iS,iT).SubjectYaw = M6_new(start_row:end_row,7);
        dataSorted(iS,iT).SubjectPitch = M6_new(start_row:end_row,6);
        dataSorted(iS,iT).SubjectRoll = M6_new(start_row:end_row,8);
        dataSorted(iS,iT).obstX = M6_new(start_row:end_row,9);
        dataSorted(iS,iT).obstY = M6_new(start_row:end_row,10);
        dataSorted(iS,iT).GoalX = M6_new(start_row:end_row,11);
        dataSorted(iS,iT).GoalY = M6_new(start_row:end_row,12);
        if dataSorted(iS,iT).info.Ang ~= 0
            dataSorted(iS,iT).Xfilt = filter_butter(hZ,dataSorted(iS,iT).SubjectX);
            dataSorted(iS,iT).Zfilt = filter_butter(hZ,dataSorted(iS,iT).SubjectZ);
            dataSorted(iS,iT).P(:,1) = dataSorted(iS,iT).Xfilt;
            dataSorted(iS,iT).P(:,2) = dataSorted(iS,iT).Zfilt;
            dataSorted(iS,iT).A1P(:,1) = dataSorted(iS,iT).obstX;
            dataSorted(iS,iT).A1P(:,2) = dataSorted(iS,iT).obstY;
            for i = 1:length(dataSorted(iS,iT).P)
                dataSorted(iS,iT).distance(1,i) = pdist2(dataSorted(iS,iT).P(i,:), dataSorted(iS,iT).A1P(i,:), 'euclidean');
                dataSorted(iS,iT).theta(1,i) = 2*atan(dataSorted(iS,iT).info.Size/dataSorted(iS,iT).distance(1,i));
            end
        end
        fprintf('Completed subject: %d, trial: %d\n', iS, iT);
    end
end

for iS = 7
    for iT = 1:N_TRIALS
        [start_row, end_row, info_row] = get_trial_rows(iT, pt7, height(M7_new));
        dataSorted(iS,iT).info.TNum = M7_new(info_row,2);
        dataSorted(iS,iT).info.Ang = M7_new(info_row,4);
        dataSorted(iS,iT).info.Speed = M7_new(info_row,6);
        dataSorted(iS,iT).info.Size = M7_new(info_row,8);
        dataSorted(iS,iT).Time = M7_new(start_row:end_row, 1);
        dataSorted(iS,iT).SubjectX = M7_new(start_row:end_row,3);
        dataSorted(iS,iT).SubjectZ = M7_new(start_row:end_row,5);
        dataSorted(iS,iT).SubjectYaw = M7_new(start_row:end_row,7);
        dataSorted(iS,iT).SubjectPitch = M7_new(start_row:end_row,6);
        dataSorted(iS,iT).SubjectRoll = M7_new(start_row:end_row,8);
        dataSorted(iS,iT).obstX = M7_new(start_row:end_row,9);
        dataSorted(iS,iT).obstY = M7_new(start_row:end_row,10);
        dataSorted(iS,iT).GoalX = M7_new(start_row:end_row,11);
        dataSorted(iS,iT).GoalY = M7_new(start_row:end_row,12);
        if dataSorted(iS,iT).info.Ang ~= 0
            dataSorted(iS,iT).Xfilt = filter_butter(hZ,dataSorted(iS,iT).SubjectX);
            dataSorted(iS,iT).Zfilt = filter_butter(hZ,dataSorted(iS,iT).SubjectZ);
            dataSorted(iS,iT).P(:,1) = dataSorted(iS,iT).Xfilt;
            dataSorted(iS,iT).P(:,2) = dataSorted(iS,iT).Zfilt;
            dataSorted(iS,iT).A1P(:,1) = dataSorted(iS,iT).obstX;
            dataSorted(iS,iT).A1P(:,2) = dataSorted(iS,iT).obstY;
            for i = 1:length(dataSorted(iS,iT).P)
                dataSorted(iS,iT).distance(1,i) = pdist2(dataSorted(iS,iT).P(i,:), dataSorted(iS,iT).A1P(i,:), 'euclidean');
                dataSorted(iS,iT).theta(1,i) = 2*atan(dataSorted(iS,iT).info.Size/dataSorted(iS,iT).distance(1,i));
            end
        end
        fprintf('Completed subject: %d, trial: %d\n', iS, iT);
    end
end

for iS = 8
    for iT = 1:N_TRIALS
        [start_row, end_row, info_row] = get_trial_rows(iT, pt8, height(M8_new));
        dataSorted(iS,iT).info.TNum = M8_new(info_row,2);
        dataSorted(iS,iT).info.Ang = M8_new(info_row,4);
        dataSorted(iS,iT).info.Speed = M8_new(info_row,6);
        dataSorted(iS,iT).info.Size = M8_new(info_row,8);
        dataSorted(iS,iT).Time = M8_new(start_row:end_row, 1);
        dataSorted(iS,iT).SubjectX = M8_new(start_row:end_row,3);
        dataSorted(iS,iT).SubjectZ = M8_new(start_row:end_row,5);
        dataSorted(iS,iT).SubjectYaw = M8_new(start_row:end_row,7);
        dataSorted(iS,iT).SubjectPitch = M8_new(start_row:end_row,6);
        dataSorted(iS,iT).SubjectRoll = M8_new(start_row:end_row,8);
        dataSorted(iS,iT).obstX = M8_new(start_row:end_row,9);
        dataSorted(iS,iT).obstY = M8_new(start_row:end_row,10);
        dataSorted(iS,iT).GoalX = M8_new(start_row:end_row,11);
        dataSorted(iS,iT).GoalY = M8_new(start_row:end_row,12);
        if dataSorted(iS,iT).info.Ang ~= 0
            dataSorted(iS,iT).Xfilt = filter_butter(hZ,dataSorted(iS,iT).SubjectX);
            dataSorted(iS,iT).Zfilt = filter_butter(hZ,dataSorted(iS,iT).SubjectZ);
            dataSorted(iS,iT).P(:,1) = dataSorted(iS,iT).Xfilt;
            dataSorted(iS,iT).P(:,2) = dataSorted(iS,iT).Zfilt;
            dataSorted(iS,iT).A1P(:,1) = dataSorted(iS,iT).obstX;
            dataSorted(iS,iT).A1P(:,2) = dataSorted(iS,iT).obstY;
            for i = 1:length(dataSorted(iS,iT).P)
                dataSorted(iS,iT).distance(1,i) = pdist2(dataSorted(iS,iT).P(i,:), dataSorted(iS,iT).A1P(i,:), 'euclidean');
                dataSorted(iS,iT).theta(1,i) = 2*atan(dataSorted(iS,iT).info.Size/dataSorted(iS,iT).distance(1,i));
            end
        end
        fprintf('Completed subject: %d, trial: %d\n', iS, iT);
    end
end

for iS = 9
    for iT = 1:N_TRIALS
        [start_row, end_row, info_row] = get_trial_rows(iT, pt9, height(M9_new));
        dataSorted(iS,iT).info.TNum = M9_new(info_row,2);
        dataSorted(iS,iT).info.Ang = M9_new(info_row,4);
        dataSorted(iS,iT).info.Speed = M9_new(info_row,6);
        dataSorted(iS,iT).info.Size = M9_new(info_row,8);
        dataSorted(iS,iT).Time = M9_new(start_row:end_row, 1);
        dataSorted(iS,iT).SubjectX = M9_new(start_row:end_row,3);
        dataSorted(iS,iT).SubjectZ = M9_new(start_row:end_row,5);
        dataSorted(iS,iT).SubjectYaw = M9_new(start_row:end_row,7);
        dataSorted(iS,iT).SubjectPitch = M9_new(start_row:end_row,6);
        dataSorted(iS,iT).SubjectRoll = M9_new(start_row:end_row,8);
        dataSorted(iS,iT).obstX = M9_new(start_row:end_row,9);
        dataSorted(iS,iT).obstY = M9_new(start_row:end_row,10);
        dataSorted(iS,iT).GoalX = M9_new(start_row:end_row,11);
        dataSorted(iS,iT).GoalY = M9_new(start_row:end_row,12);
        if dataSorted(iS,iT).info.Ang ~= 0
            dataSorted(iS,iT).Xfilt = filter_butter(hZ,dataSorted(iS,iT).SubjectX);
            dataSorted(iS,iT).Zfilt = filter_butter(hZ,dataSorted(iS,iT).SubjectZ);
            dataSorted(iS,iT).P(:,1) = dataSorted(iS,iT).Xfilt;
            dataSorted(iS,iT).P(:,2) = dataSorted(iS,iT).Zfilt;
            dataSorted(iS,iT).A1P(:,1) = dataSorted(iS,iT).obstX;
            dataSorted(iS,iT).A1P(:,2) = dataSorted(iS,iT).obstY;
            for i = 1:length(dataSorted(iS,iT).P)
                dataSorted(iS,iT).distance(1,i) = pdist2(dataSorted(iS,iT).P(i,:), dataSorted(iS,iT).A1P(i,:), 'euclidean');
                dataSorted(iS,iT).theta(1,i) = 2*atan(dataSorted(iS,iT).info.Size/dataSorted(iS,iT).distance(1,i));
            end
        end
        fprintf('Completed subject: %d, trial: %d\n', iS, iT);
    end
end

for iS = 10
    for iT = 1:N_TRIALS
        [start_row, end_row, info_row] = get_trial_rows(iT, pt10, height(M10_new));
        dataSorted(iS,iT).info.TNum = M10_new(info_row,2);
        dataSorted(iS,iT).info.Ang = M10_new(info_row,4);
        dataSorted(iS,iT).info.Speed = M10_new(info_row,6);
        dataSorted(iS,iT).info.Size = M10_new(info_row,8);
        dataSorted(iS,iT).Time = M10_new(start_row:end_row, 1);
        dataSorted(iS,iT).SubjectX = M10_new(start_row:end_row,3);
        dataSorted(iS,iT).SubjectZ = M10_new(start_row:end_row,5);
        dataSorted(iS,iT).SubjectYaw = M10_new(start_row:end_row,7);
        dataSorted(iS,iT).SubjectPitch = M10_new(start_row:end_row,6);
        dataSorted(iS,iT).SubjectRoll = M10_new(start_row:end_row,8);
        dataSorted(iS,iT).obstX = M10_new(start_row:end_row,9);
        dataSorted(iS,iT).obstY = M10_new(start_row:end_row,10);
        dataSorted(iS,iT).GoalX = M10_new(start_row:end_row,11);
        dataSorted(iS,iT).GoalY = M10_new(start_row:end_row,12);
        if dataSorted(iS,iT).info.Ang ~= 0
            dataSorted(iS,iT).Xfilt = filter_butter(hZ,dataSorted(iS,iT).SubjectX);
            dataSorted(iS,iT).Zfilt = filter_butter(hZ,dataSorted(iS,iT).SubjectZ);
            dataSorted(iS,iT).P(:,1) = dataSorted(iS,iT).Xfilt;
            dataSorted(iS,iT).P(:,2) = dataSorted(iS,iT).Zfilt;
            dataSorted(iS,iT).A1P(:,1) = dataSorted(iS,iT).obstX;
            dataSorted(iS,iT).A1P(:,2) = dataSorted(iS,iT).obstY;
            for i = 1:length(dataSorted(iS,iT).P)
                dataSorted(iS,iT).distance(1,i) = pdist2(dataSorted(iS,iT).P(i,:), dataSorted(iS,iT).A1P(i,:), 'euclidean');
                dataSorted(iS,iT).theta(1,i) = 2*atan(dataSorted(iS,iT).info.Size/dataSorted(iS,iT).distance(1,i));
            end
        end
        fprintf('Completed subject: %d, trial: %d\n', iS, iT);
    end
end

for iS = 11
    for iT = 1:N_TRIALS
        [start_row, end_row, info_row] = get_trial_rows(iT, pt11, height(M11_new));
        dataSorted(iS,iT).info.TNum = M11_new(info_row,2);
        dataSorted(iS,iT).info.Ang = M11_new(info_row,4);
        dataSorted(iS,iT).info.Speed = M11_new(info_row,6);
        dataSorted(iS,iT).info.Size = M11_new(info_row,8);
        dataSorted(iS,iT).Time = M11_new(start_row:end_row, 1);
        dataSorted(iS,iT).SubjectX = M11_new(start_row:end_row,3);
        dataSorted(iS,iT).SubjectZ = M11_new(start_row:end_row,5);
        dataSorted(iS,iT).SubjectYaw = M11_new(start_row:end_row,7);
        dataSorted(iS,iT).SubjectPitch = M11_new(start_row:end_row,6);
        dataSorted(iS,iT).SubjectRoll = M11_new(start_row:end_row,8);
        dataSorted(iS,iT).obstX = M11_new(start_row:end_row,9);
        dataSorted(iS,iT).obstY = M11_new(start_row:end_row,10);
        dataSorted(iS,iT).GoalX = M11_new(start_row:end_row,11);
        dataSorted(iS,iT).GoalY = M11_new(start_row:end_row,12);
        if dataSorted(iS,iT).info.Ang ~= 0
            dataSorted(iS,iT).Xfilt = filter_butter(hZ,dataSorted(iS,iT).SubjectX);
            dataSorted(iS,iT).Zfilt = filter_butter(hZ,dataSorted(iS,iT).SubjectZ);
            dataSorted(iS,iT).P(:,1) = dataSorted(iS,iT).Xfilt;
            dataSorted(iS,iT).P(:,2) = dataSorted(iS,iT).Zfilt;
            dataSorted(iS,iT).A1P(:,1) = dataSorted(iS,iT).obstX;
            dataSorted(iS,iT).A1P(:,2) = dataSorted(iS,iT).obstY;
            for i = 1:length(dataSorted(iS,iT).P)
                dataSorted(iS,iT).distance(1,i) = pdist2(dataSorted(iS,iT).P(i,:), dataSorted(iS,iT).A1P(i,:), 'euclidean');
                dataSorted(iS,iT).theta(1,i) = 2*atan(dataSorted(iS,iT).info.Size/dataSorted(iS,iT).distance(1,i));
            end
        end
        fprintf('Completed subject: %d, trial: %d\n', iS, iT);
    end
end

for iS = 12
    for iT = 1:N_TRIALS
        [start_row, end_row, info_row] = get_trial_rows(iT, pt12, height(M12_new));
        dataSorted(iS,iT).info.TNum = M12_new(info_row,2);
        dataSorted(iS,iT).info.Ang = M12_new(info_row,4);
        dataSorted(iS,iT).info.Speed = M12_new(info_row,6);
        dataSorted(iS,iT).info.Size = M12_new(info_row,8);
        dataSorted(iS,iT).Time = M12_new(start_row:end_row, 1);
        dataSorted(iS,iT).SubjectX = M12_new(start_row:end_row,3);
        dataSorted(iS,iT).SubjectZ = M12_new(start_row:end_row,5);
        dataSorted(iS,iT).SubjectYaw = M12_new(start_row:end_row,7);
        dataSorted(iS,iT).SubjectPitch = M12_new(start_row:end_row,6);
        dataSorted(iS,iT).SubjectRoll = M12_new(start_row:end_row,8);
        dataSorted(iS,iT).obstX = M12_new(start_row:end_row,9);
        dataSorted(iS,iT).obstY = M12_new(start_row:end_row,10);
        dataSorted(iS,iT).GoalX = M12_new(start_row:end_row,11);
        dataSorted(iS,iT).GoalY = M12_new(start_row:end_row,12);
        if dataSorted(iS,iT).info.Ang ~= 0
            dataSorted(iS,iT).Xfilt = filter_butter(hZ,dataSorted(iS,iT).SubjectX);
            dataSorted(iS,iT).Zfilt = filter_butter(hZ,dataSorted(iS,iT).SubjectZ);
            dataSorted(iS,iT).P(:,1) = dataSorted(iS,iT).Xfilt;
            dataSorted(iS,iT).P(:,2) = dataSorted(iS,iT).Zfilt;
            dataSorted(iS,iT).A1P(:,1) = dataSorted(iS,iT).obstX;
            dataSorted(iS,iT).A1P(:,2) = dataSorted(iS,iT).obstY;
            for i = 1:length(dataSorted(iS,iT).P)
                dataSorted(iS,iT).distance(1,i) = pdist2(dataSorted(iS,iT).P(i,:), dataSorted(iS,iT).A1P(i,:), 'euclidean');
                dataSorted(iS,iT).theta(1,i) = 2*atan(dataSorted(iS,iT).info.Size/dataSorted(iS,iT).distance(1,i));
            end
        end
        fprintf('Completed subject: %d, trial: %d\n', iS, iT);
    end
end

%% Rotate trials
for iS = 1:N_SUBJECTS
    for iT = 1:N_TRIALS
        dataSorted(iS,iT).info.NSubject = iS;
        if dataSorted(iS,iT).info.Ang ~= 0
            t_leng = length(dataSorted(iS,iT).Time);
            center = repmat(0, t_leng, 2);
            if dataSorted(iS,iT).GoalX == -9
                sign = -1;
            else
                sign = 1;
            end
            sPart = dataSorted(iS,iT).P - center;
            sA1 = dataSorted(iS,iT).A1P - center;
            soPart = sPart*R;
            soA1 = sA1*R;
            dataSorted(iS,iT).rotP = soPart + center;
            dataSorted(iS,iT).rotA1P = soA1 + center;
            dataSorted(iS, iT).pFilt(:,1) = dataSorted(iS, iT).rotP(:,1);
            dataSorted(iS, iT).pFilt(:,2) = dataSorted(iS, iT).rotP(:,2);
            dataSorted(iS, iT).pFilt(60:end,1) = filter_butter(hZ, dataSorted(iS, iT).rotP(60:end,1));
            dataSorted(iS, iT).pFilt(60:end,2) = filter_butter(hZ, dataSorted(iS, iT).rotP(60:end,2));
            %saving obstacle position into single variable in dataSorted
            f0 = dataSorted(iS,iT).pFilt(1,:);
            %.P is the final rotated positional data!!
            %dataSorted(iS,iT).P(:,2) = sign*(dataSorted(iS,iT).P(:,2)- f0(2));
            dataSorted(iS,iT).A1P(:,1) = sign*(dataSorted(iS,iT).rotA1P(:,1)- f0(1));
            dataSorted(iS,iT).A1P(:,2) = sign*(dataSorted(iS,iT).rotA1P(:,2)- f0(2));
            dataSorted(iS,iT).pFilt(:,1) = sign*(dataSorted(iS,iT).pFilt(:,1)- f0(1));
            dataSorted(iS,iT).pFilt(:,2) = sign*(dataSorted(iS,iT).pFilt(:,2)- f0(2));
            %The below three rows were added to try to handle the data loss
            %jump at the beginning
            dataSorted(iS, iT).Time = dataSorted(iS, iT).Time(1:end-1);
            dataSorted(iS,iT).pFinal(:,1) = dataSorted(iS,iT).pFilt(2:end,1);
            dataSorted(iS,iT).pFinal(:,2) = dataSorted(iS,iT).pFilt(2:end,2);
            [spd, hdn] = get_speed_heading(dataSorted(iS,iT).pFinal(:,1), dataSorted(iS,iT).pFinal(:,2), 90);
            dataSorted(iS,iT).hdn = hdn;
            dataSorted(iS,iT).spd = spd;
        end
    end
    fprintf('Completed subject: %d, trial: %d\n', iS, iT);
end

%%
clearvars -except dataSorted

save dataSorted_Exp5a_ObstAvoidSize
