% -------------------------------------------------------------------------
% Script: Compatibility_check_adjacency_matrix.m
%
% Description:
%   This script performs compatibility check and constructs an interaction
%   adjacency matrix (Figure 3c of the maintext) for the binding interaction graph (Figure 3d of the maintext). 
%   For each possible pair of combinations, two distances are computed using a pre-computed
%   pharmacophore distance matrix.


%  Geometric Compatibility Check of Binding Pairs, which generates The adjacency matrix J and J_complement:
%  The geometric compatibility of pharmacophore point(PP) pairs is assessed based on their Euclidean distances. 
%  Specifically, for two selected pharmacophore binding pairs, (L1, P1) and (L2, P2), 
%  where D1 represents the distance between PPs L1 and L2 on the ligand, and D2 represents 
%  the distance between PPs P1 and P2 on the receptor protein, the binding pairs are consider geometrically 
%  compatible if and only if the absolute difference between D1 and D2 satisfies the condition: 
%  |D1 – D2| ≤ f + 2e. If the absolute difference between the corresponding pair of distances is below the
%  threshold (f + 2e),the corrsponding element is assigned 1; otherwise, 0.
%  Here, f denotes the flexibility constant, and e signifies the interaction 
%  distance tolerance. In this study, f and e is set to 0.1 Å and 2.8 Å, respectively.

%
% Workflow:
%   1. Generate all ligand–protein pharmacophore string combinations.
%   2. Convert each combination to numeric form based on label mapping.
%   3. Compute two distances per pair using the distance matrix.
%   4. Apply the interaction threshold to form the adjacency matrix J.
%   5. Output J and its complement J_complement.
%
% Notes:
%   - This script is designed to be used after running Extract_coordinates_pairwise_distance.m
%     where the DistanceMatrix is generated, which we have also reproduced in this script.
%   - replacePatterns and compareAndConcatenate are local helper functions.
%
% Expected runtime: ~1 second (MATLAB, CPU)
% --------------------------------------------------------------------------------------------------------------------------------------- 



%% Generate all combinations of ligand-side and protein-side pharmacophore binding pairs
clc,clear;
tic;

% Define ligand-side pharmacophore point labels (lowercase)
String1_values = {'hp1', 'hp2', 'hp3','hp4','hp5','ha1'}; % Ligand pharmacophore labels, 6 PPs on the ligand
% Define protein-side harmacophore point labels (uppercase)
String2_values = {'HP1', 'HP2', 'HP3','HP4', 'HP5','HD1','HD2'}; % Protein pharmacophore labels, 7 PPs on the protein

% Precompute sizes
size1 = size(String1_values, 2); % Count of ligand pharmacophore labels
size2 = size(String2_values, 2); % Count of protein pharmacophore labels

% Preallocate cell array to store all combinations
combinations = cell(size1*size2, 1); % Each entry becomes e.g., 'hp1HP1'
index = 1;  % Index for filling combinations

% Generate all ligand–protein concatenated combinations
for i = 1:length(String1_values)    % Loop over ligand labels
    for j = 1:length(String2_values) % Loop over protein labels
        combinations{index} = [String1_values{i}, String2_values{j}];  % Concatenate strings
        index = index + 1;  % Increment index
    end
end

% The combination of PP pairs will be (corrsponds to nodes 1~42, respectively, shown as Figure 3b of the maintext): 
% 'hp1HP1','hp1HP2','hp1HP3','hp1HP4','hp1HP5','hp1HD1','hp1HD2','hp2HP1','hp2HP2','hp2HP3','hp2HP4','hp2HP5','hp2HD1','hp2HD2',
% 'hp3HP1','hp3HP2','hp3HP3','hp3HP4','hp3HP5','hp3HD1','hp3HD2','hp4HP1','hp4HP2','hp4HP3','hp4HP4','hp4HP5','hp4HD1','hp4HD2',
% 'hp5HP1','hp5HP2','hp5HP3','hp5HP4','hp5HP5','hp5HD1','hp5HD2','ha1HP1','ha1HP2','ha1HP3','ha1HP4','ha1HP5','ha1HD1','ha1HD2'.





%% ------------------------------------------------------------------------
% Load precomputed PP pairwise distance matrix
% (from the output of Extract_coordinates_pairwise_distance.m)
% ------------------------------------------------------------------------
DistanceMatrix=  [      0   20.0722   22.6127   10.4383   14.8783   10.9496   14.8931   15.5257   13.3347   13.2472    3.7965    9.6504   14.2628
                  20.0722         0   11.9910   21.3126   32.5873   22.7246   24.4162   18.2274   21.7904   22.5376   20.9224   20.5442   30.7876
                  22.6127   11.9910         0   18.0202   36.6403   28.4966   32.5445   22.9000   29.1186   30.5464   23.6533   18.6840   35.9702
                  10.4383   21.3126   18.0202         0   22.1937   20.5552   24.8449   21.2104   23.0820   23.3756   10.9018    3.8837   22.9735
                  14.8783   32.5873   36.6403   22.1937         0   17.7437   14.1598   27.3442   19.3098   15.2958   13.0030   20.2503    3.9308
                  10.9496   22.7246   28.4966   20.5552   17.7437         0   11.7288   10.5580    3.9552    8.8999   13.4970   20.4276   15.1274
                  14.8931   24.4162   32.5445   24.8449   14.1598   11.7288         0   20.2843   10.0137    3.1419   14.3646   22.8275   10.4808
                  15.5257   18.2274   22.9000   21.2104   27.3442   10.5580   20.2843         0   11.1227   17.1529   18.9048   22.2450   25.1167
                  13.3347   21.7904   29.1186   23.0820   19.3098    3.9552   10.0137   11.1227         0    6.9919   15.3609   22.5558   16.1500
                  13.2472   22.5376   30.5464   23.3756   15.2958    8.8999    3.1419   17.1529    6.9919         0   13.4647   21.7315   11.7497
                   3.7965   20.9224   23.6533   10.9018   13.0030   13.4970   14.3646   18.9048   15.3609   13.4647         0    8.8203   12.7977
                   9.6504   20.5442   18.6840    3.8837   20.2503   20.4276   22.8275   22.2450   22.5558   21.7315    8.8203         0   21.0327
                  14.2628   30.7876   35.9702   22.9735    3.9308   15.1274   10.4808   25.1167   16.1500   11.7497   12.7977   21.0327         0];

distance_matrix = round(DistanceMatrix, 3);  % Round to three decimals

% Initialize a size1*size2 x size1*size2 matrix to store the results
D1 = zeros(size1*size2, size1*size2);  % Preallocate first distance matrix D1
D2 = zeros(size1*size2, size1*size2);  % Preallocate second distance matrix D2

% Compute the matrix values
for i = 1:length(combinations)
    for j = 1:length(combinations)

        % Convert the combinations to number strings
        newString1 = replacePatterns(combinations{i});
        newString2 = replacePatterns(combinations{j});

        % Compute two distances via helper function
        [distance1, distance2] = compareAndConcatenate(newString1, newString2, distance_matrix);
        D1(i, j) = distance1; 
        D2(i, j) = distance2; 
    end
end

%% Perform Compatibility check and generate adjacency matrix J
e = 2.8; % the interaction distance tolerance parameter
f = 0.1; % the flexibility constant

% Initialize adjacency matrix
J = zeros(size1*size2, size1*size2);

% Fill adjacency matrix based on compatibility condition
for i = 1:size1*size2
    for j = 1:size1*size2
        if abs(D1(i, j) - D2(i, j)) <= (2*e + f)
            J(i, j) = 1;   % Mark as compatible
        end
    end
end

% Zero out diagonal (no self-interaction)
for i = 1:size1*size2
    J(i, i) = 0;
end


figure
heatmap(J); %shown as Figure 3c of the manuscript


% Compute complement matrix
J_complement=1-J;
fprintf("J_complement\n");
J_complement(logical(eye(size(J_complement))))=0
% heatmap(J_complement);

toc;

%% ------------------------------------------------------------------------
% Helper function 1: replace pharmacophore labels with numeric codes
% ------------------------------------------------------------------------
function convertedStr = replacePatterns(inputStr)
% Mapping between pharmacophore labels and numeric codes
   patternMap = {'HP2', '01';'HP5', '02';'HD2', '03';'HP1', '04'; 'HP4', '05';'HP3', '06'; 'HD1', '07';'hp5', '08';'hp3', '09';'ha1', '10';'hp2', '11';'hp1', '12';'hp4', '13'};
    
    % Initialize the converted string with the input string
    convertedStr = inputStr;
    
    % Loop over the map and replace each pattern
    for i = 1:size(patternMap, 1)
        oldPattern = patternMap{i, 1};
        newPattern = patternMap{i, 2};
        convertedStr = strrep(convertedStr, oldPattern, newPattern);
    end
end

%% ------------------------------------------------------------------------
% Helper function 2: compare the first and second two characters of str1 and str2, to detemined the distance
% ------------------------------------------------------------------------
function [result1, result2] = compareAndConcatenate(str1, str2, b)
    % Compare the first two characters of str1 and str2
    if strcmp(str1(1:2), str2(1:2))
        result1 = 0; % If they match, return 1
    else
        % Convert the first two characters from each string to numbers
        x1 = str2double(str1(1:2));
        y1 = str2double(str2(1:2));
        % Return the value at the corresponding location in matrix b
        result1 = b(x1, y1);
    end

        % Compare the second two characters of str1 and str2
    if strcmp(str1(3:4), str2(3:4))
        result2 = 0; % If they match, return 1
    else
        % Convert the first two characters from each string to numbers
        x2 = str2double(str1(3:4));
        y2 = str2double(str2(3:4));
        % Return the value at the corresponding location in matrix b
        result2 = b(x2, y2);
    end
    
end