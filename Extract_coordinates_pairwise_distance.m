% --------------------------------------------------------------------------------------------------------------------------------------- 
% Script: Extract_coordinates_pairwise_distance.m
% Description:
%   This script extracts the 3D coordinates of user-specified pharmacophore
%   points from a PDB file (7arm.pdb) and computes the full pairwise Euclidean 
%   distance matrix among these points. The script identifies atoms based on 
%   their PDB record strings (e.g., 'ATOM   2057', 'HETATM10922'), parses their
%   Cartesian coordinates using standard fixed-column PDB formatting, and
%   stores the resulting coordinates in an ordered matrix. A square
%   distance matrix is then generated, where each entry (i,j) represents
%   the Euclidean distance between pharmacophore points i and j.
%
% Workflow:
%   1. File Input:
%       - Opens the specified PDB file for reading.
%       - Validates successful file access.
%       - Defines a list of atom records corresponding to the pharmacophore
%         points of interest.
%
%   2. Coordinate Extraction:
%       - Reads the PDB file line by line.
%       - Compares each line against the target atom record identifiers.
%       - For matching lines, extracts the x, y, z coordinates from the
%         appropriate PDB columns and appends them to the coordinate matrix.
%
%   3. Distance Computation:
%       - Initializes an N×N matrix, where N is the number of selected atoms.
%       - Computes Euclidean distances between all distinct atom pairs.
%       - Outputs the complete distance matrix.
%
% Notes:
%   - Coordinate extraction follows the standard PDB fixed-column format:
%         X: columns 31–38
%         Y: columns 39–46
%         Z: columns 47–54
%   - Euclidean distance is computed as:
%         d = sqrt((x2 − x1)^2 + (y2 − y1)^2 + (z2 − z1)^2)
%
% Expected runtime: ~0.2 second (MATLAB, CPU)
% --------------------------------------------------------------------------------------------------------------------------------------- 


clc;clear;
tic;
% Open the pdb file
fileID = fopen('7arm.pdb', 'r');              % Open the PDB file for reading

% Check if the file is opened successfully
if fileID == -1
    error('File cannot be opened');           % Throw an error if the file failed to open
end

% --------------------------------------------------------------------------------------------------------------------------------------- 
% Select the pharmacophore points of interest.
% Select order is: first, pharmacophore points (in the PDB file, this refers to "ATOM") on the protein side, 
%                         then arranged in ascending order of their numerical values; 
%                  second, pharmacophore points (in the PDB file, this refers to "HETATM") on the ligand side, and  
%                         then arranged in ascending order of their numerical values; 
% Notes: Spacing must strictly match the original PDB file. 
% For example, 'ATOM    358' and 'ATOM    389' have four spaces in the middle; 
%              'ATOM   2057' has three spaces in the middle; and 'HETATM10932' has no spaces in the middle.

% Define target pharmacophore points:
 records = [ "ATOM   2057","ATOM   2758","ATOM   2805","ATOM   3214" ,"ATOM   3313","ATOM   5034","ATOM   9421","HETATM10922","HETATM10932","HETATM10938","HETATM10952","HETATM10959","HETATM10966"]; 
% Mapping between pharmacophore labels and index order    
 patternMap = {'HP2', '01';'HP5', '02';'HD2', '03';'HP1', '04'; 'HP4', '05';'HP3', '06'; 'HD1', '07';'hp5', '08';'hp3', '09';'ha1', '10';'hp2', '11';'hp1', '12';'hp4', '13'};
% --------------------------------------------------------------------------------------------------------------------------------------- 
      

% Initialize the matrix used to store coordinates
Coordinates = [];   % Empty array to store x,y,z coordinates

% Read the PDB file line by line 
while ~feof(fileID)                          % Loop until end of file
    line = fgets(fileID);                    % Read a line from the PDB file
    
    % Check if it is the record of pharmacophore points of interest
    for record = records                     % Loop through selected atom identifiers
        if contains(line, record)            % Match the line to the target record
            % Extract x, y, z coordinates (columns follow standard PDB fixed format)
            x = str2double(line(31:38));     % Extract X coordinate
            y = str2double(line(39:46));     % Extract Y coordinate
            z = str2double(line(47:54));     % Extract Z coordinate
             % Add to coordinate matrix 
            Coordinates = [Coordinates; x, y, z];   % Append coordinates as a new row
        end
    end
end

% Close file
fclose(fileID);                              % Close the PDB file after reading

% Display the results matrix
fprintf("Coordinates of pharmacophore points\n");
disp(Coordinates);                           % Print all extracted coordinates

% The Coordinates matrix will be displayed as:  
%   100.3890  104.4800  107.2090 - HP2
%   111.2000  108.9630   90.9020 - HP5
%   115.0430   97.6160   91.4140 - HD2
%   106.6270   96.1110  107.2770 - HP1
%    96.3530  111.2670  119.8190 - HP4 
%    91.8160  109.6240  102.7440 - HP3
%    96.7790  118.9120  107.9080 - HD1
%    93.5240  105.2170   93.3030 - hp5
%    92.3590  113.1420  101.0200 - hp3
%    96.0710  116.9000  105.6010 - ha1
%   102.7290  105.9750  109.7980 - hp2
%   108.3470   99.2520  108.7800 - hp1
%    95.3630  113.6810  116.8790 - hp4
% --------------------------------------------------------------------------------------------------------------------------------------- 


% Initialize distance matrix
sizeOfDisMatrix = size(records, 2);          % Number of extracted pharmacophore points 
DistanceMatrix = zeros(sizeOfDisMatrix, sizeOfDisMatrix); % Preallocate distance matrix

% Calculate the distances between all coordinate pairs
for i = 1:size(Coordinates, 1)
    for j = 1:size(Coordinates, 1)
        if i ~= j
            x1 = Coordinates(i, 1);
            y1 = Coordinates(i, 2);
            z1 = Coordinates(i, 3);
            x2 = Coordinates(j, 1);
            y2 = Coordinates(j, 2);
            z2 = Coordinates(j, 3);

            % Euclidean distance calculation
            DistanceMatrix(i, j) = sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2);  % Fill symmetric distance entry
        end
    end
end

% Display distance matrix
fprintf("Pairwise distance matrix of pharmacophore points\n");
disp(DistanceMatrix);

% The DistanceMatrix will be displayed as:  
%          HP2       HP5       HD2       HP1       HP4       HP3       HD1       hp5       hp3       ha1       hp2       hp1       hp4
% HP2          0   20.0722   22.6127   10.4383   14.8783   10.9496   14.8931   15.5257   13.3347   13.2472    3.7965    9.6504   14.2628
% HP5    20.0722         0   11.9910   21.3126   32.5873   22.7246   24.4162   18.2274   21.7904   22.5376   20.9224   20.5442   30.7876
% HD2    22.6127   11.9910         0   18.0202   36.6403   28.4966   32.5445   22.9000   29.1186   30.5464   23.6533   18.6840   35.9702
% HP1    10.4383   21.3126   18.0202         0   22.1937   20.5552   24.8449   21.2104   23.0820   23.3756   10.9018    3.8837   22.9735
% HP4    14.8783   32.5873   36.6403   22.1937         0   17.7437   14.1598   27.3442   19.3098   15.2958   13.0030   20.2503    3.9308
% HP3    10.9496   22.7246   28.4966   20.5552   17.7437         0   11.7288   10.5580    3.9552    8.8999   13.4970   20.4276   15.1274
% HD1    14.8931   24.4162   32.5445   24.8449   14.1598   11.7288         0   20.2843   10.0137    3.1419   14.3646   22.8275   10.4808
% hp5    15.5257   18.2274   22.9000   21.2104   27.3442   10.5580   20.2843         0   11.1227   17.1529   18.9048   22.2450   25.1167
% hp3    13.3347   21.7904   29.1186   23.0820   19.3098    3.9552   10.0137   11.1227         0    6.9919   15.3609   22.5558   16.1500
% ha1    13.2472   22.5376   30.5464   23.3756   15.2958    8.8999    3.1419   17.1529    6.9919         0   13.4647   21.7315   11.7497
% hp2     3.7965   20.9224   23.6533   10.9018   13.0030   13.4970   14.3646   18.9048   15.3609   13.4647         0    8.8203   12.7977
% hp1     9.6504   20.5442   18.6840    3.8837   20.2503   20.4276   22.8275   22.2450   22.5558   21.7315    8.8203         0   21.0327
% hp4    14.2628   30.7876   35.9702   22.9735    3.9308   15.1274   10.4808   25.1167   16.1500   11.7497   12.7977   21.0327         0
% --------------------------------------------------------------------------------------------------------------------------------------- 

toc;
