% optimizeDSS.m
% Dynamic Spectrum Sharing Parameter Optimization Script
% This script optimizes the four spectrum sharing parameters (initialNRSplit, upperBound, lowerBound, allocStep)
% for the DSS simulation using a coordinate descent approach with binary search for one-dimensional optimization.
% The optimization is performed sequentially for each parameter, fixing the others.
% For each parameter, it performs a binary search within a range defined by initial step to find the value that maximizes the objective function.
%
% Revised Objective Function (based on literature review for balanced LTE-NR coexistence):
%   Score = mean( weighted_throughput × avg_rat_fairness / (1 + blocking_prob) ) + mean(cell_edge_rate)
% Where:
%   - weighted_throughput = trafficMix × NR_throughput + (1 - trafficMix) × LTE_throughput
%     (load-proportional aggregation to favor NR in high-mix scenarios)
%   - avg_rat_fairness = existing global Jain_fairness (proxy for average per-RAT fairness;
%     for precise per-RAT Jain (0.5 × Jain_LTE + 0.5 × Jain_NR), modify DSS.m to compute and store
%     jainLTE and jainNR in results struct using user-specific throughputs)
% This balances total performance, fairness, blocking, and edge metrics while addressing NR bias.
%
% Assumptions and Setup:
% - Parameter bounds: initialNRSplit [0.1, 0.9], upperBound [0.6, 1.0], lowerBound [0.0, 0.4], allocStep [0.01, 0.2]
%   (Adjusted to ensure feasibility: upper > lower + 0.1, etc.)
% - Initial values: initialNRSplit=0.5, upperBound=0.9, lowerBound=0.5, allocStep=0.05
% - Initial step size: 0.1 for splits/bounds, 0.02 for step
% - Binary search precision: epsilon = 0.001
% - Max iterations per parameter: 20 (to prevent excessive computation)
% - Requires DSS.m in the path.
%
% Usage:
%   Run this script directly. It will output the optimized parameters and final objective score.
%   Results are stored in 'optimizedResults.mat' for further analysis.
%
% Tuning Notes:
% - Adjust initial values, steps, bounds, or objective function as needed.
% - Increase max iterations for finer optimization (but longer runtime).
% - The script assumes DSS returns a struct array with fields: trafficMix, mobilitySpeed, lteThroughput, etc.

%% Initialization
clear; clc; close all;

% Initial parameter values
params.initialNRSplit = 0.5;
params.upperBound = 0.8;
params.lowerBound = 0.4;
params.allocStep = 0.2;

% Parameter names and initial step sizes (for binary search range initialization)
paramNames = {'initialNRSplit', 'upperBound', 'lowerBound', 'allocStep'};
initialSteps = [0.3, 0.2, 0.2, 0.1];  % Initial range half-width for binary search

% Bounds for each parameter [min, max]
paramBounds = [0.1, 0.9;   % initialNRSplit
               0.6, 1.0;   % upperBound
               0.0, 0.4;   % lowerBound
               0.001, 0.4]; % allocStep

% Optimization settings
epsilon = 0.0001;  % Binary search stopping criterion
maxIter = 50;     % Max binary search iterations per parameter
numOptCycles = 20; % Number of full cycles over all parameters (for multi-pass refinement)

% Display initial parameters
fprintf('Starting Optimization with Initial Parameters:\n');
disp(params);
fprintf('\n');

%% Optimization Loop
bestScore = -Inf;
optimizedParams = params;

for cycle = 1:numOptCycles
    fprintf('Optimization Cycle %d/%d\n', cycle, numOptCycles);
    
    for pIdx = 1:length(paramNames)
        paramName = paramNames{pIdx};
        initStep = initialSteps(pIdx);
        lb = paramBounds(pIdx, 1);
        ub = paramBounds(pIdx, 2);
        
        % Initial test: current +/- initStep (clamped to bounds)
        currentVal = params.(paramName);
        testLow = max(lb, currentVal - initStep);
        testHigh = min(ub, currentVal + initStep);
        
        % Evaluate at testLow and testHigh
        tempParamsLow = params;
        tempParamsLow.(paramName) = testLow;
        scoreLow = evaluateObjective(tempParamsLow);
        
        tempParamsHigh = params;
        tempParamsHigh.(paramName) = testHigh;
        scoreHigh = evaluateObjective(tempParamsHigh);
        
        % Select better as new current (and set search range)
        if scoreLow > scoreHigh
            currentVal = testLow;
            searchHigh = testHigh;
            searchLow = max(lb, testLow - initStep);  % Extend low if possible
        else
            currentVal = testHigh;
            searchLow = testLow;
            searchHigh = min(ub, testHigh + initStep);  % Extend high if possible
        end
        
        % Enforce constraint: upperBound > lowerBound + 0.1
        if strcmp(paramName, 'upperBound')
            currentVal = max(currentVal, params.lowerBound + 0.1);
            searchLow = max(searchLow, params.lowerBound + 0.1);
        elseif strcmp(paramName, 'lowerBound')
            currentVal = min(currentVal, params.upperBound - 0.1);
            searchHigh = min(searchHigh, params.upperBound - 0.1);
        end
        
        params.(paramName) = currentVal;
        
        % Binary search refinement
        iter = 0;
        while (searchHigh - searchLow) > epsilon && iter < maxIter
            midVal = (searchLow + searchHigh) / 2;
            
            % Temporary params for mid
            tempParamsMid = params;
            tempParamsMid.(paramName) = midVal;
            
            % Enforce constraints for mid
            if strcmp(paramName, 'upperBound')
                midVal = max(midVal, params.lowerBound + 0.1);
            elseif strcmp(paramName, 'lowerBound')
                midVal = min(midVal, params.upperBound - 0.1);
            end
            tempParamsMid.(paramName) = midVal;
            
            scoreMid = evaluateObjective(tempParamsMid);
            scoreCurrent = evaluateObjective(params);
            
            if scoreMid > scoreCurrent
                % Move towards mid
                if midVal > currentVal
                    searchLow = currentVal;
                else
                    searchHigh = currentVal;
                end
                currentVal = midVal;
                params.(paramName) = currentVal;
            else
                % Shrink towards current
                if midVal > currentVal
                    searchHigh = midVal;
                else
                    searchLow = midVal;
                end
            end
            
            iter = iter + 1;
        end
        
        fprintf('  Optimized %s: %.4f (Score: %.4f)\n', paramName, params.(paramName), evaluateObjective(params));
    end
    
    cycleScore = evaluateObjective(params);
    if cycleScore > bestScore
        bestScore = cycleScore;
        optimizedParams = params;
    end
end

%% Final Output
fprintf('\nOptimization Complete.\n');
fprintf('Best Parameters:\n');
disp(optimizedParams);
fprintf('Best Objective Score: %.4f\n', bestScore);

% Run final simulation with best params and save
finalResults = DSS(optimizedParams.initialNRSplit, optimizedParams.upperBound, ...
                   optimizedParams.lowerBound, optimizedParams.allocStep);
save('optimizedResults.mat', 'optimizedParams', 'finalResults', 'bestScore');

fprintf('\nFinal simulation results saved to optimizedResults.mat\n');

%% Objective Evaluation Function
function score = evaluateObjective(params)
    % Run DSS simulation
    results = DSS(params.initialNRSplit, params.upperBound, ...
                  params.lowerBound, params.allocStep);
    
    % Extract metrics across scenarios
    trafficMixes = [results.trafficMix];
    lteTP = [results.lteThroughput];
    nrTP = [results.nrThroughput];
    fairness = [results.jainFairness];  % Global Jain as proxy for avg_rat_fairness
    edgeRate = [results.cellEdgeRate];
    blocking = [results.blockingProb];
    
    % Compute weighted throughput (load-proportional)
    weightedTP = trafficMixes .* nrTP + (1 - trafficMixes) .* lteTP;
    
    % Revised composite score
    throughputFairnessTerm = weightedTP .* fairness ./ (1 + blocking);
    score = mean(throughputFairnessTerm) + mean(edgeRate);
    
    % Note: For precise avg_rat_fairness = 0.5 * Jain_LTE + 0.5 * Jain_NR,
    % modify DSS.m to compute and store per-RAT Jains in results (using LTE/NR user throughputs).
end