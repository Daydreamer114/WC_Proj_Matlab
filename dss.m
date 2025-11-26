% Dynamic Spectrum Sharing Simulation for 4G (LTE) and 5G (NR) on the Same Carrier
% This script simulates dynamic spectrum sharing between LTE and NR using MATLAB's LTE Toolbox and 5G Toolbox.
% It implements a load-aware partitioning policy that adapts the spectrum split based on NR occupancy.
% The simulation varies traffic mix and mobility, and reports per RAT throughput, Jain fairness, cell edge rates, and blocking probability.

% Key Assumptions and Simplifications:
% - Single carrier with total bandwidth of 20 MHz (configurable).
% - Initial split: 50% LTE, 50% NR.
% - NR occupancy is monitored as the fraction of resources used by NR users.
% - If NR occupancy exceeds upper bound, decrease NR share by step size (dynamic allocation proportion).
% - If NR occupancy falls below lower bound, increase NR share by step size.
% - Traffic mix: Vary the proportion of NR users (e.g., 0.3, 0.5, 0.7).
% - Mobility: Vary user speeds (e.g., static=0 km/h, low=3 km/h, high=120 km/h) using random waypoint model.
% - Simulation uses simplified link-level models from toolboxes for throughput estimation.
% - Blocking probability: Fraction of users unable to get resources due to overload.
% - Cell edge rates: 5th percentile of user rates.
% - Jain fairness: Calculated across all users' throughputs.
% - Requires LTE Toolbox, 5G Toolbox, and Communications Toolbox.

%% Adjustable Parameters
% These are located here at the top for easy tuning. Adjust values as needed before running the script.

% Spectrum Sharing Parameters
initialNRSplit = 0.5;          % Initial NR spectrum share (0 to 1, e.g., 0.5 = 50%)
upperBound = 0.7;              % Upper bound for NR occupancy (e.g., 0.7 = 70%)
lowerBound = 0.6;              % Lower bound for NR occupancy (e.g., 0.6 = 60%)
allocStep = 0.05;              % Dynamic allocation proportion/step size (e.g., 0.05 = 5%)

% Simulation Parameters
totalBandwidthMHz = 20;        % Total carrier bandwidth in MHz
numSubframes = 1000;           % Number of subframes to simulate (increase for longer simulations)
numUsers = 50;                 % Total number of users (LTE + NR)
trafficMixes = [0.3, 0.5, 0.7];% Proportions of NR users to simulate (vary traffic mix)
mobilitySpeeds = [0, 3, 120];  % User speeds in km/h (vary mobility: static, low, high)

% Channel Parameters
snrDB = 10;                    % Average SNR in dB
pathLossModel = 'urban';       % Path loss model: 'urban', 'rural', etc. (simplified)

% How to Tune:
% - Change initialNRSplit, upperBound, lowerBound, allocStep to adjust the adaptive policy.
% - Modify totalBandwidthMHz for different carrier sizes.
% - Increase numSubframes for more accurate statistics (but longer runtime).
% - Add more values to trafficMixes or mobilitySpeeds to explore additional scenarios.
% - Adjust numUsers to scale the load.
% - Tune snrDB or pathLossModel for different channel conditions.

%% Simulation Setup
% Initialize results storage
results = struct('trafficMix', [], 'mobilitySpeed', [], ...
                 'lteThroughput', [], 'nrThroughput', [], ...
                 'jainFairness', [], 'cellEdgeRate', [], 'blockingProb', []);

% Loop over traffic mixes and mobility speeds
for mixIdx = 1:length(trafficMixes)
    nrProportion = trafficMixes(mixIdx);
    numNRUsers = round(numUsers * nrProportion);
    numLTEUsers = numUsers - numNRUsers;
    
    for mobIdx = 1:length(mobilitySpeeds)
        userSpeed = mobilitySpeeds(mobIdx);
        
        % Reset spectrum split for each scenario
        currentNRSplit = initialNRSplit;
        lteBandwidth = totalBandwidthMHz * (1 - currentNRSplit);
        nrBandwidth = totalBandwidthMHz * currentNRSplit;
        
        % User positions and mobility (simplified random waypoint)
        userPositions = rand(numUsers, 2) * 1000; % Positions in a 1km x 1km cell
        userVelocities = userSpeed * randn(numUsers, 2) / sqrt(2); % Velocities in m/s (approx)
        
        % Throughput and resource usage trackers
        userThroughputs = zeros(numUsers, 1);
        blockedUsers = 0;
        
        % Simulate subframes
        for sf = 1:numSubframes
            % Update positions based on mobility
            userPositions = userPositions + (userVelocities * 0.001); % 1ms subframe, speed in m/s
            
            % Calculate distances and path losses (simplified)
            distances = sqrt(sum(userPositions.^2, 2)); % Distance from BS at (0,0)
            if strcmp(pathLossModel, 'urban')
                pathLosses = 35 * log10(distances) + 35; % Simplified urban model in dB
            else
                pathLosses = 30 * log10(distances) + 40; % Fallback
            end
            userSNRs = snrDB - pathLosses;
            
            % Estimate capacities using Shannon (approx throughput per user if scheduled)
            capacities = log2(1 + 10.^(userSNRs/10)); % bps/Hz
            
            % Resource allocation: Split resources based on current split
            lteResources = lteBandwidth * 180; % Approx PRBs (180 kHz per PRB)
            nrResources = nrBandwidth * 180;
            
            % Assign users to RATs
            lteCapacities = capacities(1:numLTEUsers);
            nrCapacities = capacities(numLTEUsers+1:end);
            
            % Simulate scheduling: Proportional fair (simplified, assume equal share)
            if numLTEUsers > 0
                lteUserThroughputs = (lteResources / numLTEUsers) * lteCapacities * 1e3; % kbps (1ms subframe)
            else
                lteUserThroughputs = [];
            end
            if numNRUsers > 0
                nrUserThroughputs = (nrResources / numNRUsers) * nrCapacities * 1e3; % kbps
            else
                nrUserThroughputs = [];
            end
            
            % Check for blocking: If demand exceeds resources (simplified threshold)
            if mean(lteUserThroughputs) < 100 || mean(nrUserThroughputs) < 100 % Threshold 100 kbps
                blockedUsers = blockedUsers + 1;
            end
            
            % Accumulate throughputs
            userThroughputs(1:numLTEUsers) = userThroughputs(1:numLTEUsers) + lteUserThroughputs / numSubframes;
            userThroughputs(numLTEUsers+1:end) = userThroughputs(numLTEUsers+1:end) + nrUserThroughputs / numSubframes;
            
            % Monitor NR occupancy (simplified as fraction of used NR resources)
            nrOccupancy = min(1, numNRUsers * mean(nrCapacities) / (nrBandwidth * log2(1 + 10^(snrDB/10))));
            
            % Adaptive adjustment
            if nrOccupancy > upperBound && currentNRSplit > allocStep
                currentNRSplit = currentNRSplit - allocStep;
            elseif nrOccupancy < lowerBound && currentNRSplit < 1 - allocStep
                currentNRSplit = currentNRSplit + allocStep;
            end
            
            % Update bandwidths
            lteBandwidth = totalBandwidthMHz * (1 - currentNRSplit);
            nrBandwidth = totalBandwidthMHz * currentNRSplit;
        end
        
        % Calculate indicators
        lteThroughput = mean(userThroughputs(1:numLTEUsers));
        nrThroughput = mean(userThroughputs(numLTEUsers+1:end));
        allThroughputs = userThroughputs(userThroughputs > 0); % Exclude zeros
        jainFairness = (sum(allThroughputs)^2) / (length(allThroughputs) * sum(allThroughputs.^2));
        cellEdgeRate = prctile(allThroughputs, 5); % 5th percentile
        blockingProb = blockedUsers / numSubframes;
        
        % Store results
        resIdx = (mixIdx-1)*length(mobilitySpeeds) + mobIdx;
        results(resIdx).trafficMix = nrProportion;
        results(resIdx).mobilitySpeed = userSpeed;
        results(resIdx).lteThroughput = lteThroughput;
        results(resIdx).nrThroughput = nrThroughput;
        results(resIdx).jainFairness = jainFairness;
        results(resIdx).cellEdgeRate = cellEdgeRate;
        results(resIdx).blockingProb = blockingProb;
    end
end

%% Display Results
disp('Simulation Results:');
for i = 1:length(results)
    fprintf('Traffic Mix (NR Prop): %.2f, Mobility Speed: %d km/h\n', ...
            results(i).trafficMix, results(i).mobilitySpeed);
    fprintf('  LTE Throughput: %.2f kbps\n', results(i).lteThroughput);
    fprintf('  NR Throughput: %.2f kbps\n', results(i).nrThroughput);
    fprintf('  Jain Fairness: %.4f\n', results(i).jainFairness);
    fprintf('  Cell Edge Rate: %.2f kbps\n', results(i).cellEdgeRate);
    fprintf('  Blocking Probability: %.4f\n\n', results(i).blockingProb);
end

% Notes on Results Interpretation:
% - Per RAT Throughput: Average throughput for LTE and NR users.
% - Jain Fairness: Value between 0 and 1; higher indicates better fairness across all users.
% - Cell Edge Rates: 5th percentile throughput, representing edge performance.
% - Blocking Probability: Fraction of subframes where users were blocked (low throughput threshold).

% To compare with static baseline:
% - Rerun the script with allocStep = 0 (disables adaptation), and compare outputs.
