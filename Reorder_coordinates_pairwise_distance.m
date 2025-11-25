% --------------------------------------------------------------------------------------------------------------------------------------- 
% Script: Reorder_pharmacophore_coordinates_and_pairwise_distance.m
% Description:
%   The results of this script provides SOURCE DATA for Supplementary Figure 2. 
%   This section reorders the extracted pharmacophore point coordinates according
%   to a specified canonical atom-label sequence. The original label list
%   corresponds to the order in which coordinates were extracted from the
%   PDB file. The targetOrder list defines a desired standardized ordering
%   (e.g., lower-case ligand labels first, followed by upper-case protein
%   labels). A new coordinate matrix is constructed accordingly, and a new
%   pairwise Euclidean distance matrix is computed from the reordered list.
%
%    
%
% Workflow:
%   1. Define original and target label sequences.
%   2. Reorder coordinate rows by matching labels to targetOrder.
%   3. Construct reordered coordinate matrix.
%   4. Compute full pairwise Euclidean distance matrix for the new ordering.
%
% Expected runtime: ~0.1 second (MATLAB, CPU)
% --------------------------------------------------------------------------------------------------------------------------------------- 
clc;clear;
tic;
% Original Coordinates matrix
Coordinates=[ 100.3890  104.4800  107.2090
              111.2000  108.9630   90.9020
              115.0430   97.6160   91.4140
              106.6270   96.1110  107.2770
               96.3530  111.2670  119.8190
               91.8160  109.6240  102.7440
               96.7790  118.9120  107.9080
               93.5240  105.2170   93.3030
               92.3590  113.1420  101.0200
               96.0710  116.9000  105.6010
              102.7290  105.9750  109.7980
              108.3470   99.2520  108.7800
               95.3630  113.6810  116.8790];

% Original labels
labels = {'HP2','HP5', 'HD2', 'HP1', 'HP4', 'HP3',  'HD1', 'hp5', 'hp3', 'ha1', 'hp2', 'hp1', 'hp4'}; % Original label order

% Target ordering to standardize pharmacophore index sequence
targetOrder = {'hp1', 'hp2', 'hp3','hp4','hp5','ha1', 'HP1', 'HP2', 'HP3','HP4', 'HP5','HD1','HD2'}; % Desired label order

% Initialize reordered coordinate matrix (same size as original)
reordered_Coordinates = zeros(size(Coordinates));           % Preallocate matrix

% Perform reordering by matching each target label to its original index
for i = 1:length(targetOrder)                               % Loop through target order list
    index = find(strcmp(labels, targetOrder{i}));           % Locate matching label in original list
    reordered_Coordinates(i, :) = Coordinates(index, :);    % Copy matching coordinate row
end

% Display reordered coordinate matrix
fprintf("Reordered coordinates of pharmacophore points\n");
disp(reordered_Coordinates);

% The Coordinates matrix will be displayed as:  
%   108.3470   99.2520  108.7800 - hp1
%   102.7290  105.9750  109.7980 - hp2
%    92.3590  113.1420  101.0200 - hp3
%    95.3630  113.6810  116.8790 - hp4
%    93.5240  105.2170   93.3030 - hp5
%    96.0710  116.9000  105.6010 - ha1
%   106.6270   96.1110  107.2770 - HP1
%   100.3890  104.4800  107.2090 - HP2
%    91.8160  109.6240  102.7440 - HP3
%    96.3530  111.2670  119.8190 - HP4
%   111.2000  108.9630   90.9020 - HP5
%    96.7790  118.9120  107.9080 - HD1
%   115.0430   97.6160   91.4140 - HD2
% --------------------------------------------------------------------------------------------------------------------------------------- 

% Compute Pairwise Distance Matrix Based on Reordered Coordinates
% Initialize distance matrix for the new ordering
sizeOfDisMatrix = size(Coordinates, 2);                                 % Number of pharmacophore points
reordered_DistanceMatrix = zeros(sizeOfDisMatrix, sizeOfDisMatrix); % Preallocate N×N matrix

% Compute Euclidean distances between all coordinate pairs
for i = 1:size(reordered_Coordinates, 1)
    for j = 1:size(reordered_Coordinates, 1)
        if i ~= j
            x1 = reordered_Coordinates(i, 1);
            y1 = reordered_Coordinates(i, 2);
            z1 = reordered_Coordinates(i, 3);
            x2 = reordered_Coordinates(j, 1);
            y2 = reordered_Coordinates(j, 2);
            z2 = reordered_Coordinates(j, 3);

            % Euclidean distance calculation
            reordered_DistanceMatrix(i, j) = sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2);
        end
    end
end

% Display reordered distance matrix
fprintf("Reordered pairwise distance matrix of pharmacophore points\n");
disp(reordered_DistanceMatrix);                                     % Print full pairwise distance table
 
% The reordered_DistanceMatrix will be displayed as (SOURCE DATA FOR Supplementary Figure 2):  
%         hp1       hp2       hp3       hp4       hp5       ha1       HP1       HP2       HP3       HP4       HP5       HD1       HD2
% hp1         0    8.8203   22.5558   21.0327   22.2450   21.7315    3.8837    9.6504   20.4276   20.2503   20.5442   22.8275   18.6840
% hp2    8.8203         0   15.3609   12.7977   18.9048   13.4647   10.9018    3.7965   13.4970   13.0030   20.9224   14.3646   23.6533
% hp3   22.5558   15.3609         0   16.1500   11.1227    6.9919   23.0820   13.3347    3.9552   19.3098   21.7904   10.0137   29.1186
% hp4   21.0327   12.7977   16.1500         0   25.1167   11.7497   22.9735   14.2628   15.1274    3.9308   30.7876   10.4808   35.9702
% hp5   22.2450   18.9048   11.1227   25.1167         0   17.1529   21.2104   15.5257   10.5580   27.3442   18.2274   20.2843   22.9000
% ha1   21.7315   13.4647    6.9919   11.7497   17.1529         0   23.3756   13.2472    8.8999   15.2958   22.5376    3.1419   30.5464
% HP1    3.8837   10.9018   23.0820   22.9735   21.2104   23.3756         0   10.4383   20.5552   22.1937   21.3126   24.8449   18.0202
% HP2    9.6504    3.7965   13.3347   14.2628   15.5257   13.2472   10.4383         0   10.9496   14.8783   20.0722   14.8931   22.6127
% HP3   20.4276   13.4970    3.9552   15.1274   10.5580    8.8999   20.5552   10.9496         0   17.7437   22.7246   11.7288   28.4966
% HP4   20.2503   13.0030   19.3098    3.9308   27.3442   15.2958   22.1937   14.8783   17.7437         0   32.5873   14.1598   36.6403
% HP5   20.5442   20.9224   21.7904   30.7876   18.2274   22.5376   21.3126   20.0722   22.7246   32.5873         0   24.4162   11.9910
% HD1   22.8275   14.3646   10.0137   10.4808   20.2843    3.1419   24.8449   14.8931   11.7288   14.1598   24.4162         0   32.5445
% HD2   18.6840   23.6533   29.1186   35.9702   22.9000   30.5464   18.0202   22.6127   28.4966   36.6403   11.9910   32.5445         0
% --------------------------------------------------------------------------------------------------------------------------------------- 
toc;