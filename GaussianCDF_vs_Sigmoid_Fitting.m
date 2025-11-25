% ============================================================================================================================
% Script: GaussianCDF_vs_Sigmoid_Fitting.m
%
% Description:
%   This script visualizes the relationship between Gaussian cumulative distribution
%   functions (CDFs) with varying standard deviation σ, and their corresponding
%   logistic-sigmoid approximations parameterized by β. The Gaussian CDFs model the
%   statistical behavior of the Gaussian RNG-based p-bit input dependence, while the
%   sigmoid functions approximate these CDFs for analytical modeling.
%
% The script supports two coloring modes:
%   - PAIR_SAME_COLOR = true  → Each Gaussian–sigmoid pair uses the same color.
%   - PAIR_SAME_COLOR = false → Each curve is assigned a unique color (10 distinct colors in total).
% The figure provides visual confirmation of the quality of the Gaussian→sigmoid
%
%   approximation across various σ values relevant to the stochastic p-bit device.
% Expected runtime: ~1.5 seconds (MATLAB, CPU)
% ============================================================================================================================
clc;clear
tic;
figure('Color','w'); hold on; grid on;  % Create figure with white background, enable grid

% Define standard deviations σ and their corresponding fitted β parameters
sigmas = [0.2 0.5 1.0 1.5 2.0];   
betas  = [8.50 3.35 1.65 1.15 0.85];

% Define input range for evaluation (x-axis)
xgrid  = linspace(-6,6,2000);   % High resolution grid for smooth curves

% ------------------------------------------------------------------------
% Coloring setup
% ------------------------------------------------------------------------
PAIR_SAME_COLOR = false;               % If true → same color for each Gaussian–sigmoid pair
colsPair = lines(numel(sigmas));       % Color palette with 5 colors (for pairs)
cols10   = lines(numel(sigmas)*2);     % Color palette with 10 distinct colors (for all curves)

% Pre-allocate storage for plot handles and legend labels
hPair   = gobjects(numel(sigmas), 2);  % Handles for plotting each pair (Gaussian, Sigmoid)
labPair = strings(numel(sigmas), 2);   % Labels for legend entries

% ============================================================================================================================
% Loop over each σ–β pair and plot Gaussian CDF and corresponding sigmoid fit
% ============================================================================================================================
for k = 1:numel(sigmas)
    s    = sigmas(k);  % Current σ value
    beta = betas(k);   % Corresponding fitted β value

    % Assign colors depending on mode
    if PAIR_SAME_COLOR
        c1 = colsPair(k,:); c2 = c1;               % Use same color for both Gaussian and sigmoid
    else
        c1 = cols10(2*k-1,:); c2 = cols10(2*k,:);  % Use two distinct colors for Gaussian and sigmoid
    end

    % -----------------------------
    % Gaussian CDF (solid line)
    % -----------------------------
    cdf_gauss   = 0.5 * (1 + erf(xgrid / (sqrt(2) * s))); % Closed-form Gaussian CDF
    hPair(k,1)  = plot(xgrid, cdf_gauss, '-', 'LineWidth', 2, 'Color', c1);
    labPair(k,1)= sprintf('Gaussian CDF ($\\sigma=%.1f$)', s);  % Legend entry with LaTeX

    % -----------------------------
    % Sigmoid approximation (dashed line)
    % -----------------------------
    cdf_sig   = 1 ./ (1 + exp(-beta * xgrid));     % Logistic sigmoid
    hPair(k,2)= plot(xgrid, cdf_sig, '--', 'LineWidth', 2, 'Color', c2);
    labPair(k,2)= sprintf('Sigmoid ($\\beta=%.2f$)', beta); % Legend entry with LaTeX
end

% ============================================================================================================================
% Figure annotation
% ============================================================================================================================
xlabel('$s_{in}$','Interpreter','latex','FontSize',16,'FontAngle','italic');  % X-axis label
ylabel('$P(x_i^{new}=1)$','Interpreter','latex','FontSize',16,'FontAngle','italic'); % Y-axis label
title('Gaussian CDFs (solid) and Sigmoid fits (dashed)', ...
      'Interpreter','latex','FontSize',18); % Figure title

set(gca,'LineWidth',1.2);  % Set axis line width for better visibility
xlim([-6 6]);              % Fix x-axis range

% ============================================================================================================================
% Legend setup
% ============================================================================================================================
% Flatten the handle and label arrays in row-major order, then display legend
% with two columns for compact layout.
legend(hPair(:), cellstr(labPair(:)), ...
       'Interpreter','latex', 'NumColumns', 2, ...
       'Location','best', 'Box','on', 'FontSize',12);
toc;