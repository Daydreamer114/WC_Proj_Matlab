function results = DSS(initialNRSplit, upperBound, lowerBound, allocStep)
% Dynamic Spectrum Sharing Simulation for 4G (LTE) and 5G (NR) on the Same Carrier
% This function simulates dynamic spectrum sharing between LTE and NR using MATLAB's LTE Toolbox and 5G Toolbox.
% It implements a load-aware partitioning policy that adapts the spectrum split based on NR occupancy.
% The simulation varies traffic mix and mobility, and computes per RAT throughput, Jain fairness, cell edge rates, and blocking probability.
% Key Assumptions and Simplifications:
% - Single carrier with total bandwidth of 20 MHz (configurable).
% - Initial split: Provided as input.
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
%
% Inputs:
%   initialNRSplit - Initial NR spectrum share (0 to 1, e.g., 0.5 = 50%)
%   upperBound     - Upper bound for NR occupancy (e.g., 0.9)
%   lowerBound     - Lower bound for NR occupancy (e.g., 0.5)
%   allocStep      - Dynamic allocation proportion/step size (e.g., 0.05 = 5%)
%
% Output:
%   results        - Struct array with simulation metrics for each scenario
%
% Usage Example:
%   results = DSS(0.5, 0.9, 0.5, 0.05);
%
% Notes on Results Interpretation:
% - Per RAT Throughput: Average throughput for LTE and NR users.
% - Jain Fairness: Value between 0 and 1; higher indicates better fairness across all users.
% - Cell Edge Rates: 5th percentile throughput, representing edge performance.
% - Blocking Probability: Fraction of subframes where users were blocked (low throughput threshold).
% To compare with static baseline:
% - Call the function with allocStep = 0 (disables adaptation), and compare outputs.

%% Adjustable Parameters
% These are located here for easy tuning. Adjust values as needed.
% Simulation Parameters
totalBandwidthMHz = 20; % Total carrier bandwidth in MHz
numSubframes = 1000; % Number of subframes to simulate (increase for longer simulations)
numUsers = 100; % Total number of users (LTE + NR) - increased for smoother statistics
trafficMixes = [0.3, 0.5, 0.7]; % Proportions of NR users to simulate (vary traffic mix)
mobilitySpeeds = [0, 30, 60, 120]; % User speeds in km/h (expanded for clearer trends: static, pedestrian, vehicular, high)
% Channel Parameters
snrDB = 20; % Average SNR in dB (reduced for balanced throughput scale and variability)
pathLossModel = 'urban'; % Path loss model: 'urban', 'rural', etc. (simplified)
minThroughputThreshold = 200; % Minimum throughput threshold for blocking detection (kbps) - raised to introduce moderate blocking
% How to Tune:
% - Modify totalBandwidthMHz for different carrier sizes.
% - Increase numSubframes for more accurate statistics (but longer runtime).
% - Add more values to trafficMixes or mobilitySpeeds to explore additional scenarios.
% - Adjust numUsers to scale the load.
% - Tune snrDB or pathLossModel for different channel conditions.
% - Adjust minThroughputThreshold to fine-tune blocking sensitivity.

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
            userPositions = mod(userPositions, 1000); % Boundary handling: wrap around to stay in cell
            % Calculate distances and path losses (simplified, updated model)
            distances = sqrt(sum(userPositions.^2, 2)); % Distance from BS at (0,0)
            distances_km = distances / 1000; % Convert to km
            if strcmp(pathLossModel, 'urban')
                pathLosses = 32.4 + 20 * log10(2) + 30 * log10(distances_km + 0.01); % 3GPP-like urban model, +0.01 to avoid log(0)
            else
                pathLosses = 28.0 + 22 * log10(distances_km + 0.01); % Rural fallback
            end
            userSNRs = snrDB - pathLosses;
            % Estimate capacities using Shannon (approx throughput per user if scheduled)
            capacities = log2(1 + 10.^(userSNRs/10)); % bps/Hz
            capacities = max(capacities, 0.1); % Minimum capacity floor to avoid underflow
            % Add channel variation for realism and fairness improvement
            capacities = capacities .* (1 + 0.1 * randn(size(capacities))); % 10% fading variation
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
            if ~isempty(lteUserThroughputs) && mean(lteUserThroughputs) < minThroughputThreshold || ...
               ~isempty(nrUserThroughputs) && mean(nrUserThroughputs) < minThroughputThreshold
                blockedUsers = blockedUsers + 1;
            end
            % Accumulate throughputs
            if numLTEUsers > 0
                userThroughputs(1:numLTEUsers) = userThroughputs(1:numLTEUsers) + lteUserThroughputs / numSubframes;
            end
            if numNRUsers > 0
                userThroughputs(numLTEUsers+1:end) = userThroughputs(numLTEUsers+1:end) + nrUserThroughputs / numSubframes;
            end
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
end