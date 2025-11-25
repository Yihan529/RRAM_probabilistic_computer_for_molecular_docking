% ============================================================================================================================
% Script: MWCP_solver_42nodes_pcomputing.m
%
% Description:
%   This script implements an illustrative probabilistic computing solver for the
%   42-node Maximum Weighted Clique Problem (MWCP) derived from a molecular docking
%   instance (LolA–LolCDE complex, PDB ID: 7arm). The MWCP is reformulated into a QUBO
%   instance using the adjacency complement matrix and pharmacophore-based bias
%   potentials. The solver operates using our Gaussian RNG-driven p-bit update rule,
%   executing sequential asynchronous updates over a predefined number of iterations.
%
%   The purpose of this script is to provide a reference implementation that clarifies 
%   the workflow used in our hardware solver, ensuring transparency and reproducibility 
%   of the computational method.
%
%   Expected runtime: ~90 seconds (MATLAB, CPU)
%
% Workflow:
%   1. Load the QUBO interaction matrix (J_complement) for the 42-node problem.
%   2. Load the pharmacophore-derived bias vector h.
%   3. Execute probabilistic p-bit updates for multiple independent runs.
%   4. Convert each sampled state into a unique decimal code for frequency analysis.
%   5. Identify the dominant clique from each run and compare with the optimal solution: maximum weighted clique (MWC).
%   6. Report the success rate across all runs.
%
% Notes:
%   - Gaussian RNG and update logic follow the analytical derivation in Methods and Supplementary Note 8.
%   - The solver is simplied and omits the dynamic slope annealing schedule, which was executed in hardware.
% ============================================================================================================================

clc;clear;
tic

% ============================================================================================================================
% Define the QUBO problem matrix
% ============================================================================================================================
% Load J_complement matrix for the 42-node docking problem
% (from the output of Compatibility_check_adjacency_matrix.m)
% Each entry J(i,j) represents the interaction strength between node (or p-bit) i and node (or p-bit) j.
P=9; % Penalty factor P shown as Equation (4) of the main text 
J_complement_penalty=-P*  [0	1	1	1	1	1	1	1	0	1	1	1	1	1	1	1	0	0	0	0	0	1	1	0	0	0	0	0	1	1	0	0	0	0	0	1	1	0	0	0	0	0
        1	0	1	1	1	1	1	0	1	0	1	1	1	1	1	1	1	1	0	1	0	1	1	1	1	0	1	0	1	1	1	1	0	1	0	1	1	1	1	0	1	0
        1	1	0	1	1	1	1	1	0	1	1	1	0	1	0	1	1	0	0	1	1	0	1	1	0	0	1	1	0	1	1	0	0	1	1	0	1	1	0	0	1	1
        1	1	1	0	1	1	1	1	1	1	1	1	0	1	0	1	0	1	1	1	1	0	1	0	1	1	1	1	0	1	0	1	1	1	1	0	1	0	1	1	1	1
        1	1	1	1	0	1	1	1	1	1	1	1	1	0	0	0	0	1	1	0	1	0	0	0	1	1	0	1	0	0	0	1	1	0	1	0	0	0	1	1	0	1
        1	1	1	1	1	0	1	1	1	0	0	1	1	1	0	1	1	1	0	1	1	0	1	1	1	0	1	1	0	1	1	1	0	1	1	0	1	1	1	0	1	1
        1	1	1	1	1	1	0	1	1	1	1	0	1	1	0	0	1	1	1	1	1	0	0	1	1	1	1	1	0	0	1	1	1	1	1	0	0	1	1	1	1	1
        1	0	1	1	1	1	1	0	1	1	1	1	1	1	1	0	0	1	1	1	0	1	0	1	1	1	1	0	1	1	0	0	0	1	0	1	0	1	1	1	1	0
        0	1	0	1	1	1	1	1	0	1	1	1	1	1	0	1	0	0	0	0	1	0	1	0	0	1	0	1	1	1	1	0	0	0	0	0	1	0	0	1	0	1
        1	0	1	1	1	0	1	1	1	0	1	1	1	1	0	0	1	0	1	0	1	1	0	1	0	1	0	1	0	1	1	0	0	1	1	1	0	1	0	1	0	1
        1	1	1	1	1	0	1	1	1	1	0	1	1	1	1	0	0	1	1	0	1	1	0	0	1	1	0	1	0	0	0	1	1	0	1	1	0	0	1	1	0	1
        1	1	1	1	1	1	0	1	1	1	1	0	1	1	1	0	1	1	1	1	0	1	1	1	1	1	1	0	0	0	0	1	1	0	1	1	1	1	1	1	1	0
        1	1	0	0	1	1	1	1	1	1	1	1	0	1	1	0	0	0	1	1	1	1	0	0	0	1	1	1	1	0	1	0	0	1	1	1	0	0	0	1	1	1
        1	1	1	1	0	1	1	1	1	1	1	1	1	0	0	1	1	1	0	1	1	0	1	1	1	0	1	1	0	0	1	1	1	1	1	0	1	1	1	0	1	1
        1	1	0	0	0	0	0	1	0	0	1	1	1	0	0	1	1	1	1	1	1	1	1	0	1	0	1	0	1	0	1	1	1	1	1	1	0	1	1	1	1	1
        1	1	1	1	0	1	0	0	1	0	0	0	0	1	1	0	1	1	1	1	1	1	1	0	0	0	0	1	0	1	0	0	1	0	1	0	1	0	1	1	1	1
        0	1	1	0	0	1	1	0	0	1	0	1	0	1	1	1	0	1	1	1	1	0	0	1	0	1	0	1	1	0	1	1	1	0	1	1	0	1	1	1	0	1
        0	1	0	1	1	1	1	1	0	0	1	1	0	1	1	1	1	0	1	1	1	1	0	0	1	1	0	1	1	0	1	1	1	0	1	1	1	1	1	1	1	1
        0	0	0	1	1	0	1	1	0	1	1	1	1	0	1	1	1	1	0	1	1	0	0	1	1	1	1	0	1	1	1	1	1	1	0	1	1	1	1	1	1	0
        0	1	1	1	0	1	1	1	0	0	0	1	1	1	1	1	1	1	1	0	1	1	0	0	0	1	1	1	1	0	0	0	1	1	1	1	1	0	1	1	1	1
        0	0	1	1	1	1	1	0	1	1	1	0	1	1	1	1	1	1	1	1	0	0	1	1	1	0	1	1	1	1	1	1	0	1	1	1	1	1	1	0	1	1
        1	1	0	0	0	0	0	1	0	1	1	1	1	0	1	1	0	1	0	1	0	0	1	1	1	1	1	1	1	1	0	0	0	0	1	1	0	1	1	1	1	1
        1	1	1	1	0	1	0	0	1	0	0	1	0	1	1	1	0	0	0	0	1	1	0	1	1	1	1	1	1	1	1	1	0	1	0	0	1	0	0	1	0	1
        0	1	1	0	0	1	1	1	0	1	0	1	0	1	0	0	1	0	1	0	1	1	1	0	1	1	1	1	0	1	1	1	0	1	0	1	0	1	1	1	0	1
        0	1	0	1	1	1	1	1	0	0	1	1	0	1	1	0	0	1	1	0	1	1	1	1	0	1	1	1	0	1	1	1	1	1	1	1	0	1	1	1	0	1
        0	0	0	1	1	0	1	1	1	1	1	1	1	0	0	0	1	1	1	1	0	1	1	1	1	0	1	1	0	0	0	1	1	0	1	1	1	1	1	1	1	0
        0	1	1	1	0	1	1	1	0	0	0	1	1	1	1	0	0	0	1	1	1	1	1	1	1	1	0	1	0	1	1	1	0	1	1	1	0	0	0	1	1	1
        0	0	1	1	1	1	1	0	1	1	1	0	1	1	0	1	1	1	0	1	1	1	1	1	1	1	1	0	1	0	0	1	1	1	1	1	1	1	1	0	1	1
        1	1	0	0	0	0	0	1	1	0	0	0	1	0	1	0	1	1	1	1	1	1	1	0	0	0	0	1	0	1	1	1	1	1	1	1	1	0	0	0	1	0
        1	1	1	1	0	1	0	1	1	1	0	0	0	0	0	1	0	0	1	0	1	1	1	1	1	0	1	0	1	0	1	1	1	1	1	1	1	1	0	0	0	0
        0	1	1	0	0	1	1	0	1	1	0	0	1	1	1	0	1	1	1	0	1	0	1	1	1	0	1	0	1	1	0	1	1	1	1	0	1	1	0	0	0	1
        0	1	0	1	1	1	1	0	0	0	1	1	0	1	1	0	1	1	1	0	1	0	1	1	1	1	1	1	1	1	1	0	1	1	1	0	0	0	1	1	0	1
        0	0	0	1	1	0	1	0	0	0	1	1	0	1	1	1	1	1	1	1	0	0	0	0	1	1	0	1	1	1	1	1	0	1	1	0	0	0	1	1	1	0
        0	1	1	1	0	1	1	1	0	1	0	0	1	1	1	0	0	0	1	1	1	0	1	1	1	0	1	1	1	1	1	1	1	0	1	1	0	0	0	1	1	1
        0	0	1	1	1	1	1	0	0	1	1	1	1	1	1	1	1	1	0	1	1	1	0	0	1	1	1	1	1	1	1	1	1	1	0	0	0	1	1	0	1	1
        1	1	0	0	0	0	0	1	0	1	1	1	1	0	1	0	1	1	1	1	1	1	0	1	1	1	1	1	1	1	0	0	0	1	0	0	1	1	1	1	1	1
        1	1	1	1	0	1	0	0	1	0	0	1	0	1	0	1	0	1	1	1	1	0	1	0	0	1	0	1	1	1	1	0	0	0	0	1	0	1	1	1	1	1
        0	1	1	0	0	1	1	1	0	1	0	1	0	1	1	0	1	1	1	0	1	1	0	1	1	1	0	1	0	1	1	0	0	0	1	1	1	0	1	1	1	1
        0	1	0	1	1	1	1	1	0	0	1	1	0	1	1	1	1	1	1	1	1	1	0	1	1	1	0	1	0	0	0	1	1	0	1	1	1	1	0	1	1	1
        0	0	0	1	1	0	1	1	1	1	1	1	1	0	1	1	1	1	1	1	0	1	1	1	1	1	1	0	0	0	0	1	1	1	0	1	1	1	1	0	1	1
        0	1	1	1	0	1	1	1	0	0	0	1	1	1	1	1	0	1	1	1	1	1	0	0	0	1	1	1	1	0	0	0	1	1	1	1	1	1	1	1	0	1
        0	0	1	1	1	1	1	0	1	1	1	0	1	1	1	1	1	1	0	1	1	1	1	1	1	0	1	1	0	0	1	1	0	1	1	1	1	1	1	1	1	0];  % 42×42 QUBO J_complement matrix with penalty for the MWCP


% ============================================================================================================================
% Define bias vector (knowledge-based pharmacophore potential derived from the PDBbind dataset, table 1 of the main text)
% The elements of h are sequentially assigned to the following features (same with results in 'combinations'):
% 'hp1HP1','hp1HP2','hp1HP3','hp1HP4','hp1HP5','hp1HD1','hp1HD2','hp2HP1','hp2HP2','hp2HP3','hp2HP4','hp2HP5','hp2HD1','hp2HD2',
% 'hp3HP1','hp3HP2','hp3HP3','hp3HP4','hp3HP5','hp3HD1','hp3HD2','hp4HP1','hp4HP2','hp4HP3','hp4HP4','hp4HP5','hp4HD1','hp4HD2',
% 'hp5HP1','hp5HP2','hp5HP3','hp5HP4','hp5HP5','hp5HD1','hp5HD2','ha1HP1','ha1HP2','ha1HP3','ha1HP4','ha1HP5','ha1HD1','ha1HD2'.
% ============================================================================================================================
h=10*[0.0504,0.0504,0.0504,0.0504,0.0504,0.1453,0.1453,0.0504,0.0504,0.0504,0.0504,0.0504,0.1453,0.1453,...
      0.0504,0.0504,0.0504,0.0504,0.0504,0.1453,0.1453,0.0504,0.0504,0.0504,0.0504,0.0504,0.1453,0.1453,...
      0.0504,0.0504,0.0504,0.0504,0.0504,0.1453,0.1453,0.2317,0.2317,0.2317,0.2317,0.2317,0.6686,0.6686]; % Size = 1×42


Vdd=1;                   % Global output voltage of each node
num_nodes=42;            % Total number of nodes (p-bits)
time_step = 10000;       % Number of iterations per run/trial
G0=0.4;                  % Global scaling parameter for interaction strength
tool=0;                  % Intermediate variable to connect ∑Jij*xj with hi

% ============================================================================================================================
% Simulation loop
% ============================================================================================================================
num_runs = 100;   % Total number of independent runs
target_solution = [1,9,17,25,41];  % Expected optimal solution, i.e., the maximum weighted clique (MWC) (1-based indexing)
success_count = 0;  % Counter for successful runs

for run_no=1:1:num_runs % Total number of independent runs    % Number of runs
         Pbit = randi([0, 1], num_nodes, 1); % Initialize node states (p-bits) randomly with 0 or 1
    % Iterate through all time steps
    for i = 1:time_step-1   
        % Update each node sequentially
        for row=1:length(h)%
                    tool1=0;
                    for colomn=1:length(h) % The following loop calculates ∑Jij*xj for each node
                        if colomn>=row
                            tool1=(J_complement_penalty(row,colomn)*Pbit(colomn,1))+tool1;
                        elseif colomn<row
                            tool1=(J_complement_penalty(row,colomn)*Pbit(colomn,1))+tool1;
                        end
                        tool=tool+J_complement_penalty(row,colomn);
                    end
                    tool1=tool1+h(row); % After the calculation of ∑Jij*xj, add hi (local bias) to each node
                     Pbit(row,1)=Gaussian(tool1*G0);
        end

% % Convert the binary state vector into a unique decimal representation.This allows us to track the most frequently visited state configuration.
          Y(1, i) = sum(Pbit' .* 2.^(41:-1:0));
    end

%============================================================================================================================
% Analyze the most frequently occurring solution
%============================================================================================================================
[maxValue, maxIndex] = mode(Y); % Find the state with highest visit frequency
% Decode the binary indices of the clique (positions of '1's in the string)
maxProbability_index3{run_no,1} = find(dec2bin(maxValue, 42) == '1');
% Compute the occurrence probability (%) of the most frequent solution
maxProbability_index3{run_no,2} = maxIndex / time_step * 100;

% Check if the identified clique matches the expected solution
    if isequal(maxProbability_index3{run_no,1}, target_solution)
        success_count = success_count + 1;
    end

end

success_rate = success_count / num_runs * 100;
fprintf('Success rate of finding the optimal solution MWC: [1,9,17,25,41] over %d runs: %.2f%%\n', num_runs, success_rate);
toc;


%%
%============================================================================================================================
% Function: Gaussian p-bit update function (sigmoid-like transfer function)
% Detailed analytical derivation and fitting results of Gaussian Cumulative distribution function 
% with the sigmoid function can be found in Methods section of the maintext and Supplementary Note 8
%============================================================================================================================
% Input:
%   Iin  - Input current (interaction + bias term)
% Output:
%   Vout - Node state (0 or Vdd=1)
% Description:
%   Generates a Gaussian random number u, then compares it with Iin.
%   If Iin ≥ u, the output is 1; otherwise, the output is 0.
% This models the stochastic switching behavior of a Gaussian RNG-driven p-bit.

function Vout = Gaussian(Iin)
VDD = 1;
mu = 0;    % Mean of Gaussian RNG
sigma_u = 0.5;  % Standard deviation
u = mu + sigma_u * randn;  % Gaussian-distributed random number
% Stochastic decision rule
Vout(Iin >= u) = VDD;
Vout(Iin < u) = 0;
end

