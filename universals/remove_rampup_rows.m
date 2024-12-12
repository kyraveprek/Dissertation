function [outputMatrix] = remove_rampup_rows(inputMatrix, inputColumn)

% created: December 11, 2024
% Kyra Veprek, kyraveprek24@gmail.com
%
% This function reads in the matrices from Kyra's experiments and removes
% the rows before the obejct has appeared by looking for the places where
% it remains stationary

% Find indices where the value in inputColumn changes between rows
changeIndices = [true; diff(inputMatrix(:, inputColumn)) ~= 0];

% Extract the rows corresponding to these changes
outputMatrix = inputMatrix(changeIndices, :);
end
% row = 1;
% for iR = 1:(size(inputMatrix)-1)
%     if inputMatrix(iR,inputColumn) ~= inputMatrix(iR+1,inputColumn)
%         outputMatrix(row,:) = inputMatrix(iR,:);
%         row = row+1;
%     end
% end
