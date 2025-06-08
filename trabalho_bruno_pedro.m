%% ==========================================================================
%  PEDRO FERREIRA 2222035 | BRUNO VICENTE 2210709
%  MOBILE COMMUNICATIONS 24/25
%
%  This script simulates and analyzes mobile network deployments in two scenarios:
%  1) 5G Network coverage in Leiria city
%  2) Raytracing simulation for signal propagation at ESTG campus
%% ==========================================================================

% Clear workspace, command window, and close all figures
clc; clear all; close all;

% Start timer to measure execution time
tic

% Display program header
fprintf(['\n\n\n',...
    '============================================================\n', ...
    ' PEDRO FERREIRA 2222035 | BRUNO VICENTE 2210709\n', ...
    ' COMUNICACOES MOVEIS 24/25\n', ...
    '============================================================\n\n']);

% Present scenario options to user
disp('Choose the scenario for analysis:');
disp('1 - 5G Network in Leiria city');
disp('2 - Raytracing Simulation ESTG');
choice = 0;
% Input validation loop to ensure user selects a valid option
while ~ismember(choice, [1, 2])
    try
        choice = input('Enter 1 or 2: ');
        if isempty(choice)
            disp('No input provided. Please enter 1 or 2.');
            continue;
        elseif ~isnumeric(choice) || ~ismember(choice, [1, 2])
            disp('Invalid option! Please enter 1 or 2.');
        end
    catch
        disp('Error in input. Please enter 1 or 2.');
    end
end

switch choice
    case 1
        %% ====================== 5G NETWORK SIMULATION ======================
        fprintf('You chose: 5G Network in Leiria city\n');

        % Define frequency constants
        RURAL_FREQ = 0.7e9;  % 700 MHz for rural areas
        URBAN_FREQ = 3.6e9;  % 3.6 GHz for urban areas

        % Define base station locations in Leiria
        % ZR - Rural Zone (700 MHz)
        % ZU - Urban Zone (3600 MHz)
        baseStations = [
            39.731823, -8.798338 % ZR 1 - Rural station
            39.736940, -8.799721 % ZU 2 - Urban station
            39.739715, -8.802933 % ZU 3 - Urban station
            39.737284, -8.810495 % ZU 4 - Urban station
            39.739443, -8.819369 % ZU 5 - Urban station
            39.732543, -8.824665 % ZR 6 - Rural station
            39.725831, -8.830228 % ZU 7 - Urban station
            39.744488, -8.817963 % ZR 8 - Rural station
            39.747375, -8.815152 % ZU 9 - Urban station
            39.746381, -8.804517 % ZU 10 - Urban station

            % Additional locations (commented out)
            %39.747029, -8.809635 %  Castelo Leiria
            %39.739583, -8.871774 %  ESTG
            %39.751326, -8.781058 %  Andrinos
            %39.712350, -8.849653
            ];

        numStations = size(baseStations, 1); % Calculate total number of stations

        %% Create patch microstrip antenna array for base stations
        % Define frequency and calculate wavelength (using rural frequency for example)
        % Later in the code, appropriate frequency will be selected based on station type
        fq = RURAL_FREQ; % Using rural frequency (700 MHz) for this antenna design example

        % Create patch element
        txElement = design(patchMicrostrip, fq);
        txElement.Tilt = 90; % Set tilt angle
        txElement.TiltAxis = [0 1 0]; % Set tilt axis

        % Define array dimensions (2x2)
        ntxrow = 2;
        ntxcol = 2;
        drow = lambda/2; % Half-wavelength spacing in rows
        dcol = lambda/2; % Half-wavelength spacing in columns

        % Create Uniform Rectangular Array (URA)
        txArray = phased.URA("Size", [ntxrow ntxcol], ...
            "Element", txElement, ...
            "ElementSpacing", [drow dcol]);

        % Calculate antenna patterns for visualization
        az = -180:1:180; % Azimuth range
        el = 0; % Elevation fixed at 0°
        az_el = 0; % Azimuth for elevation pattern

        % Calculate radiation patterns
        pat_az = pattern(txArray, fq, az, el); % Azimuth pattern
        pat_el = pattern(txArray, fq, az_el, -90:1:90); % Elevation pattern

        % Create figure for antenna visualization
        figure('Position',[100 100 1200 900]);

        % 1. Display array element positions
        subplot(2,2,1);
        pos = txArray.getElementPosition;
        scatter3(pos(1,:), pos(2,:), pos(3,:), 50,'r', 'filled');
        grid on; axis equal;
        xlabel('X (m)');
        ylabel('Y (m)');
        zlabel('Z (m)');
        title('Position of URA Array Elements');
        view(3);

        % 2. Display azimuth pattern
        subplot(2,2,2);
        patternAzimuth(txArray, fq);

        % 3. Display elevation pattern
        subplot(2,2,3);
        patternElevation(txArray, fq);

        % 4. Display 3D pattern
        subplot(2,2,4);
        pattern(txArray, fq, -180:180, -90:90, 'Type', 'directivity');
        title('3D Pattern of URA Array');

        %% Define network parameters
        % Specify which stations are in urban zones (using 3.6GHz)
        zu_indices = [2,3,4,5,7,9,10];

        % Define which station will use sectorization (3 sectors)
        sectorized_index = 5;

        % Azimuth angles for the three sectors (120° apart)
        azimutes = [0, 120, 240];

        % Set heights for base stations (default=28m, station 3=40m)
        stationHeights = 28 * ones(1, size(baseStations, 1));
        stationHeights(3) = 40;  % Set station #3 to be higher (40m)

        numStations = size(baseStations,1);

        % Initialize transmitter sites array
        txSites = txsite.empty;
        siteCount = 1;  % Global site counter

        %% Let user choose antenna type
        disp('Choose the antenna type:');
        disp('1 - Isotropic');
        disp('2 - Sectorized patch array');
        antennaChoice = 0;
        % Input validation for antenna choice
        while ~ismember(antennaChoice, [1, 2])
            antennaChoice = input('Enter 1 or 2: ');
            if ~ismember(antennaChoice, [1, 2])
                disp('Invalid option! Please enter 1 or 2.');
            end
        end

        % Initialize txSites array
        txSites = txsite.empty;
        siteCount = 1;

        %% Create base stations with appropriate configurations
        for i = 1:numStations
            % Define frequency and power based on urban/rural zone
            if ismember(i, zu_indices)
                freq = 3.6e9;  % 3.6 GHz for urban areas
                power = 5;     % Lower power for urban areas (more dense)
            else
                freq = 0.7e9;  % 700 MHz for rural areas (better propagation)
                power = 20;    % Higher power for rural coverage
            end

            % Create base station with appropriate configuration
            if i == sectorized_index
                % Create 3 sectors for the sectorized station
                for sector = 1:3
                    if antennaChoice == 1
                        % Create with isotropic antenna
                        txSites(siteCount) = txsite( ...
                            'Name', sprintf('Base%d_Sector%d', i, sector), ...
                            'Latitude', baseStations(i,1), ...
                            'Longitude', baseStations(i,2), ...
                            'AntennaHeight', stationHeights(i), ...
                            'TransmitterFrequency', freq, ...
                            'TransmitterPower', power, ...
                            'AntennaAngle', azimutes(sector));
                    else
                        % Create with patch array antenna
                        txSites(siteCount) = txsite( ...
                            'Name', sprintf('Base%d_Sector%d', i, sector), ...
                            'Latitude', baseStations(i,1), ...
                            'Longitude', baseStations(i,2), ...
                            'AntennaHeight', stationHeights(i), ...
                            'TransmitterFrequency', freq, ...
                            'TransmitterPower', power, ...
                            'AntennaAngle', azimutes(sector), ...
                            'Antenna', txArray);
                    end
                    siteCount = siteCount + 1;
                end
            else
                % Create regular non-sectorized station
                if antennaChoice == 1
                    % Create with isotropic antenna
                    txSites(siteCount) = txsite( ...
                        'Name', sprintf('BaseStation%d', i), ...
                        'Latitude', baseStations(i,1), ...
                        'Longitude', baseStations(i,2), ...
                        'AntennaHeight', stationHeights(i), ...
                        'TransmitterFrequency', freq, ...
                        'TransmitterPower', power);
                else
                    % Create with patch array antenna
                    txSites(siteCount) = txsite( ...
                        'Name', sprintf('BaseStation%d', i), ...
                        'Latitude', baseStations(i,1), ...
                        'Longitude', baseStations(i,2), ...
                        'AntennaHeight', stationHeights(i), ...
                        'TransmitterFrequency', freq, ...
                        'TransmitterPower', power, ...
                        'Antenna', txArray);
                end
                siteCount = siteCount + 1;
            end
        end

        %% Select propagation model
        modelList = {'close-in', 'longley-rice', 'freespace'};
        fprintf('Available propagation models:\n');
        for idx = 1:numel(modelList)
            fprintf('  %d - %s\n', idx, modelList{idx});
        end
        selection = input('Enter the number of the propagation model: ');

        % Apply selected propagation model
        if isnumeric(selection) && selection >= 1 && selection <= numel(modelList)
            modelName = modelList{selection};
            try
                pm = propagationModel(modelName);
                disp(['Using model: ', modelName]);

                % Calculate and display coverage
                coverage(txSites, ...
                    'SignalStrengths', -110:1:-10, ...  % Signal strength range (dBm)
                    'MaxRange', 5000, ...               % Maximum range in meters
                    'PropagationModel', pm, ...         % Selected propagation model
                    'Resolution', 20);                  % Resolution in meters
            catch ME
                warning("Error applying model '%s': %s", modelName, ME.message);
            end
        else
            disp('No model selected. Script terminated.');
        end

        % Display elapsed time
        toc

        %% Generate grid of points for detailed analysis
        % Set latitude/longitude limits based on base station locations
        latLim = [min([txSites.Latitude])-0.005, max([txSites.Latitude])+0.005];
        lonLim = [min([txSites.Longitude])-0.005, max([txSites.Longitude])+0.005];

        % Create grid with reduced resolution for faster processing
        gridResolution = 50;
        [latGrid, lonGrid] = meshgrid(linspace(latLim(1), latLim(2), gridResolution), ...
            linspace(lonLim(1), lonLim(2), gridResolution));

        fprintf('Calculating RSSI for grid %dx%d...\n', gridResolution, gridResolution);

        %% Calculate RSSI for all transmitter sites
        rssiMatrix = zeros([size(latGrid), numel(txSites)]);
        for i = 1:numel(txSites)
            fprintf('Calculating RSSI for base station %d of %d...\n', i, numel(txSites));
            % Create receiver sites at all grid points
            rx = rxsite('Latitude', latGrid(:), 'Longitude', lonGrid(:));

            % Calculate signal strength at each point
            sig = sigstrength(rx, txSites(i), pm);

            % Reshape result back to grid format
            rssiMatrix(:,:,i) = reshape(sig, size(latGrid));
        end

        %% Best Server Analysis and Co-Channel Interference Study
        % Find the strongest server at each point
        [bestRSSI, bestServerIndex] = max(rssiMatrix, [], 3);

        % Calculate statistics for each base station
        fprintf('Best Server Statistics:\n');
        for i = 1:numel(txSites)
            coverage_area = sum(bestServerIndex(:) == i) / numel(bestServerIndex) * 100;
            fprintf('Station %d: %.2f%% of the area\n', i, coverage_area);
        end

        % Calculate dominance margin (difference between strongest and second strongest)
        sortedRSSI = sort(rssiMatrix, 3, 'descend');
        margin = sortedRSSI(:,:,1) - sortedRSSI(:,:,2);

        % Define threshold for border zones
        threshold = 3; % 3 dB threshold for potential handover zones

        %% Visualization 1: Best Server Map
        figure('Name', 'Best Server Map', 'Position', [100 100 700 600]);
        pcolor(lonGrid, latGrid, bestServerIndex);
        shading flat;
        colormap(parula(numel(txSites)));  % Use distinct colors for each station
        colorbar('Ticks', 1:numel(txSites), 'TickLabels', {txSites.Name});
        axis xy;
        hold on;

        % Plot base station locations with labels
        for i = 1:numel(txSites)
            plot(txSites(i).Longitude, txSites(i).Latitude, 'r^', 'MarkerSize', 10, 'LineWidth', 2);
            text(txSites(i).Longitude, txSites(i).Latitude, [' BS', num2str(i)], 'Color', 'black', 'FontWeight', 'bold');
        end
        title('Best Server Map', 'FontWeight', 'bold');
        xlabel('Longitude'); ylabel('Latitude');
        hold off;

        %% Visualization 2: Dominance Margin Map
        figure('Name', 'Dominance Margin', 'Position', [150 150 700 600]);
        pcolor(lonGrid, latGrid, margin);
        shading flat;
        colormap(hot);
        colorbar;
        title('Dominance Margin (RSSI_1 - RSSI_2) [dB]', 'FontWeight', 'bold');
        xlabel('Longitude'); ylabel('Latitude');
        hold on;

        % Highlight areas with low margin (potential handover zones)
        borderZone = margin < threshold;
        contour(lonGrid, latGrid, borderZone, [1 1], 'LineColor', 'b', 'LineWidth', 1.5);
        legend('Margin < 3 dB (Instability Zone)', 'Location', 'northeast');
        hold off;

        %% Visualization 3: Estimated SINR Map
        figure('Name', 'Estimated SINR Map', 'Position', [200 200 700 600]);

        % Calculate SINR: Signal / (Interference + Noise)
        signal = sortedRSSI(:,:,1);

        % Convert RSSI to linear power, sum, then back to dB
        total_power = sum(10.^(rssiMatrix/10), 3);

        % Calculate interference (total power minus signal)
        interference = 10*log10(total_power - 10.^(signal/10) + 1e-10);

        % SINR = Signal - Interference in dB domain
        SINR = signal - interference;

        % Plot SINR map
        pcolor(lonGrid, latGrid, SINR);
        shading flat;
        colormap(parula);
        colorbar;
        title('Estimated SINR Map [dB]', 'FontWeight', 'bold');
        xlabel('Longitude'); ylabel('Latitude');

        %% Visualization 4: Best Server Distribution
        figure('Name', 'Best Server Distribution', 'Position', [250 250 700 600]);
        counts = histcounts(bestServerIndex, 1:(numel(txSites)+1));
        b = bar(counts);
        b.FaceColor = 'flat';  % Enable individual bar coloring

        % Apply colors from the map to match bars with coverage areas
        barColors = parula(numel(txSites));
        for i = 1:numel(txSites)
            b.CData(i,:) = barColors(i,:);
        end

        xticks(1:numel(txSites));
        xticklabels({txSites.Name});
        xlabel('Base Station');
        ylabel('Number of points as Best Server');
        title('Best Server Distribution', 'FontWeight', 'bold');
        grid on;

        %% Handover Analysis along a Predefined Path
        fprintf('\n--- Handover Study ---\n');
        fprintf('Chosen path type: Predefined Path \n');

        pathName = 'Predefined Path';

        % Define path coordinates that traverse the coverage area
        latitudes = [39.747029, 39.750000, 39.755000, 39.760000, 39.765000, ...
            39.770000, 39.775000, 39.778857, 39.776000, 39.774000, ...
            39.772000, 39.770000, 39.765000, 39.760000, 39.755000, ...
            39.751326, 39.748000, 39.745000, 39.740000, 39.735000, ...
            39.730000, 39.725000, 39.720000, 39.715000, 39.712350];

        longitudes = [-8.809635, -8.808000, -8.806000, -8.804000, -8.802000, ...
            -8.800000, -8.798000, -8.780254, -8.782000, -8.784000, ...
            -8.786000, -8.788000, -8.790000, -8.792000, -8.794000, ...
            -8.781058, -8.785000, -8.790000, -8.800000, -8.810000, ...
            -8.820000, -8.830000, -8.840000, -8.845000, -8.849653];

        numPoints = length(latitudes);
        fprintf('Path defined with %d points\n', numPoints);

        % Calculate RSSI along the path for all stations
        rssiMatrix = zeros(numPoints, numel(txSites));
        fprintf('Calculating RSSI along the path...\n');
        for i = 1:numel(txSites)
            rx = rxsite('Latitude', latitudes, 'Longitude', longitudes);
            sig = sigstrength(rx, txSites(i), pm);
            sig = sig(:);
            rssiMatrix(:, i) = sig;
        end

        % Determine best server at each point
        [bestRSSI, bestServerIndex] = max(rssiMatrix, [], 2);

        %% Visualization of Handover Study
        figure('Name','Handover Study - Path, RSSI and Best Server', 'Position',[100 100 1200 800]);

        % Plot 1: Map showing path and base stations
        subplot(3,1,1);
        geoplot(latitudes, longitudes, 'b-', 'LineWidth', 2); hold on;

        % Generate distinct colors for each base station
        colors = lines(numel(txSites));

        % Plot each base station with unique color
        for i = 1:numel(txSites)
            geoplot(txSites(i).Latitude, txSites(i).Longitude, ...
                '^', 'MarkerSize', 10, ...
                'MarkerEdgeColor', 'k', ...
                'MarkerFaceColor', colors(i,:), ...
                'DisplayName', txSites(i).Name);
        end

        title(['[Map] RX Path and Base Stations - ', pathName]);
        legend('RX Path', txSites.Name, 'Location','best');
        grid on;

        % Plot 2: Active base station along the path
        subplot(3,1,2);
        plot(1:numPoints, bestServerIndex, 'k', 'LineWidth', 1.5);
        xlabel('Point on path');
        ylabel('Active Base Station');
        yticks(1:numel(txSites));
        yticklabels({txSites.Name});
        title('[Chart] Serving Station (Best Server) Along the Path');
        grid on;

        % Plot 3: RSSI values along the path
        subplot(3,1,3);
        plot(1:numPoints, rssiMatrix, 'LineWidth', 1.2);
        xlabel('Point on path');
        ylabel('RSSI (dBm)');
        title('[Chart] RSSI Received by Base Station');
        legend({txSites.Name}, 'Location','best');
        grid on;

        % Identify handover points (where best server changes)
        handoverPoints = find(diff(bestServerIndex) ~= 0);
        fprintf('\n--- Handovers Detected ---\n');
        for i = 1:length(handoverPoints)
            idx = handoverPoints(i);
            fprintf('Handover at point %d: %s → %s\n', ...
                idx+1, txSites(bestServerIndex(idx)).Name, txSites(bestServerIndex(idx+1)).Name);
        end

        %% Capacity Analysis with User Distribution
        fprintf('\n--- Capacity Analysis with 11940 Users ---\n');

        % Total number of users to distribute
        nUsers = 128642;

        % Get base station positions for user distribution
        txPositions = zeros(numel(txSites), 2);
        for i = 1:numel(txSites)
            txPositions(i, 1) = txSites(i).Latitude;
            txPositions(i, 2) = txSites(i).Longitude;
        end

        % Initialize user coordinates
        userLat = zeros(nUsers, 1);
        userLon = zeros(nUsers, 1);

        % Define user distribution: 75% in urban areas, 25% in rural areas
        urbanUsers = round(0.75 * nUsers);
        ruralUsers = nUsers - urbanUsers;

        % Define influence radii (smaller for urban areas, larger for rural)
        urbanRadius = 0.015; % Urban radius in lat/lon degrees
        ruralRadius = 0.06;  % Rural radius in lat/lon degrees

        % Define which antennas are in urban areas
        urbanAntennas = [1, 2];  % Urban base station indices
        ruralAntennas = setdiff(1:numel(txSites), urbanAntennas);  % All others are rural

        % Distribution weights for urban antennas
        urbanWeights = [0.4, 0.6]; % 40% to antenna 1, 60% to antenna 2
        urbanDistribution = round(urbanUsers * urbanWeights);

        % Ensure exact total by adjusting last value
        urbanDistribution(end) = urbanUsers - sum(urbanDistribution(1:end-1));

        currentUserIndex = 1;

        % Distribute urban users around urban base stations
        for i = 1:length(urbanAntennas)
            antennaIndex = urbanAntennas(i);
            centerLat = txPositions(antennaIndex, 1);
            centerLon = txPositions(antennaIndex, 2);

            usersForThisAntenna = urbanDistribution(i);

            % Generate uniform distribution in circular area
            r = urbanRadius * sqrt(rand(usersForThisAntenna, 1));
            theta = 2 * pi * rand(usersForThisAntenna, 1);

            % Calculate positions using polar coordinates
            userLat(currentUserIndex:currentUserIndex+usersForThisAntenna-1) = centerLat + r .* cos(theta);
            userLon(currentUserIndex:currentUserIndex+usersForThisAntenna-1) = centerLon + r .* sin(theta);

            currentUserIndex = currentUserIndex + usersForThisAntenna;
        end

        % Distribute rural users around rural base stations
        ruralUsersPerAntenna = round(ruralUsers / length(ruralAntennas));

        for i = 1:length(ruralAntennas)
            antennaIndex = ruralAntennas(i);
            centerLat = txPositions(antennaIndex, 1);
            centerLon = txPositions(antennaIndex, 2);

            % For last antenna, ensure exact count
            usersForThisAntenna = ruralUsersPerAntenna;
            if i == length(ruralAntennas)
                usersForThisAntenna = ruralUsers - (i-1) * ruralUsersPerAntenna;
            end

            % Generate uniform distribution in circular area
            r = ruralRadius * sqrt(rand(usersForThisAntenna, 1));
            theta = 2 * pi * rand(usersForThisAntenna, 1);

            userLat(currentUserIndex:currentUserIndex+usersForThisAntenna-1) = centerLat + r .* cos(theta);
            userLon(currentUserIndex:currentUserIndex+usersForThisAntenna-1) = centerLon + r .* sin(theta);

            currentUserIndex = currentUserIndex + usersForThisAntenna;
        end

        % Ensure all coordinates are within map boundaries
        userLat = max(min(userLat, latLim(2)), latLim(1));
        userLon = max(min(userLon, lonLim(2)), lonLim(1));

        % Create receiver objects for all users
        rx = rxsite('Latitude', userLat, 'Longitude', userLon);

        % Calculate RSSI for all users and all base stations
        userRSSI = zeros(nUsers, numel(txSites));
        for i = 1:numel(txSites)
            userRSSI(:, i) = sigstrength(rx, txSites(i), pm);
        end

        % Define bandwidth for capacity calculation
        bandwidths = 100e6 * ones(1, numel(txSites));  % 100 MHz for all stations
        bandwidths([1, 2]) = 100e6;  % Same bandwidth for stations 1 and 2

        % Determine best server for each user based on RSSI
        [bestRSSI, userAssignedBS] = max(userRSSI, [], 2);

        % Define noise parameters
        noise_dBm = -100;  % Noise level in dBm
        noise_mW = 10^(noise_dBm / 10);  % Convert to mW

        % Convert RSSI to linear power (mW)
        userPower_mW = 10.^(userRSSI / 10);

        % Calculate SINR for all users at all base stations
        userSINR_all = zeros(nUsers, numel(txSites));
        for u = 1:nUsers
            for bs = 1:numel(txSites)
                signal = userPower_mW(u, bs);
                interference = sum(userPower_mW(u, :)) - signal;  % All other signals are interference
                userSINR_all(u, bs) = signal / (interference + noise_mW);  % SINR calculation
            end
        end

        % Get SINR with assigned base station for each user
        userSINR = zeros(nUsers, 1);
        for u = 1:nUsers
            userSINR(u) = userSINR_all(u, userAssignedBS(u));
        end

        % Calculate Shannon capacity for each user (bits per second)
        capacityPerUser_shannon = bandwidths(userAssignedBS)' .* log2(1 + userSINR);

        % Calculate total capacity per base station
        capacityPerStation = zeros(1, numel(txSites));
        for i = 1:numel(txSites)
            usersInBS = (userAssignedBS == i);
            capacityPerStation(i) = sum(capacityPerUser_shannon(usersInBS));
        end

        % Calculate average capacity per user per base station
        capacityPerUser = zeros(1, numel(txSites));
        for i = 1:numel(txSites)
            usersInBS = (userAssignedBS == i);
            if sum(usersInBS) > 0
                capacityPerUser(i) = capacityPerStation(i) / sum(usersInBS);
            end
        end

        % Count users per base station
        capacityCount = histcounts(userAssignedBS, 1:(numel(txSites)+1));

        %% Capacity Visualization
        figure('Name', 'Capacity with 11940 Users');

        % Average capacity per user by station (Mbps)
        subplot(1, 2, 2);
        bar(capacityPerUser / 1e6);  % Convert to Mbps
        xticks(1:numel(txSites));
        xticklabels({txSites.Name});
        title('Average Capacity per User');
        xlabel('Base Station');
        ylabel('Capacity (Mbps)');

        % User distribution by station
        figure('Name', 'User Distribution by Base Station');
        bar(capacityCount);
        xticks(1:numel(txSites));
        xticklabels({txSites.Name});
        title('User Distribution by Base Station');
        xlabel('Base Station');
        ylabel('Number of Users');
        hold on;

        % Highlight urban stations
        urbanAntennas = [1, 2];
        bar(urbanAntennas, capacityCount(urbanAntennas), 'r');
        legend('Rural Stations', 'Urban Stations');

        % Print capacity statistics
        fprintf('Capacity Statistics with %d users:\n', nUsers);
        totalUrbanUsers = 0;
        totalRuralUsers = 0;

        for i = 1:numel(txSites)
            % Identify station type (urban or rural)
            if ismember(i, urbanAntennas)
                stationType = 'Urban';
                totalUrbanUsers = totalUrbanUsers + capacityCount(i);
            else
                stationType = 'Rural';
                totalRuralUsers = totalRuralUsers + capacityCount(i);
            end

            % Print statistics for each station
            fprintf('Station %d (%s): %d users (%.1f%%) - Type: %s\n', ...
                i, txSites(i).Name, capacityCount(i), 100*capacityCount(i)/nUsers, stationType);

            if capacityCount(i) > 0
                fprintf('  Total capacity: %.2f Mbps\n', capacityPerStation(i)/1e6);
                fprintf('  Average capacity per user: %.2f Mbps\n', capacityPerUser(i)/1e6);
            else
                fprintf('  [No assigned users]\n');
            end
        end

        % Print final summary
        fprintf('\nFinal Summary:\n');
        fprintf('Total users: %d\n', nUsers);
        fprintf('Users in urban areas: %d (%.1f%%)\n', totalUrbanUsers, 100*totalUrbanUsers/nUsers);

    case 2
        %% ================ RAYTRACING SIMULATION FOR ESTG ================
        fprintf("=== STARTING SIMULATION ===\n");

        %% Initial configuration
        fprintf("Setting up base station and transmitter parameters...\n");

        % Define base station coordinates
        latitude = 39.735250;
        longitude = -8.820638;

        % Set transmitter parameters
        txHeight = 3;                   % Antenna height in meters
        Ptx_dBm = 0;                    % Transmission power in dBm
        Ptx_W = 10^((Ptx_dBm - 30)/10); % Convert dBm to Watts
        frequency = 3e9;                % Frequency: 3 GHz

        % Create transmitter site object
        tx = txsite("Name", "Base Station", ...
            "Latitude", latitude, ...
            "Longitude", longitude, ...
            "AntennaHeight", txHeight, ...
            "TransmitterPower", Ptx_W, ...
            "TransmitterFrequency", frequency);

        %% Configure ray tracing propagation model
        fprintf("Creating ray tracing propagation model...\n");

        % Create ray tracing model with specific parameters
        pm = propagationModel("raytracing", ...
            "Method", "sbr",                 % Shooting and Bouncing Rays method
        "MaxNumReflections", 2,          % Allow up to 2 reflections
        "SurfaceMaterial", "concrete");  % Use concrete material properties

        % Open site viewer with buildings and terrain data
        viewer1 = siteviewer("Name", "Ray Tracing Viewer", ...
            Buildings="estg.osm",      % Load OpenStreetMap building data
        Terrain="gmted2010");      % Load terrain elevation data

        fprintf("Plotting simulated coverage using ray tracing...\n");

        % Calculate coverage using ray tracing model
        coverage(tx, ...
            "PropagationModel", pm, ...
            "MaxRange", 2000,          % Maximum range in meters
        "Resolution", 10,           % Resolution in meters
        "SignalStrengths", [-120 -100 -90 -80 -70 -60 -50]);  % Signal levels to plot

        %% Import measured data
        fprintf("\nImporting measured data from file...\n");

        % Open and read data file
        fid = fopen('gupo1.txt');
        C = textscan(fid, '%f%f%f','Delimiter',',');
        fclose(fid);

        % Extract measurement data
        Latitude = C{1};
        Longitude = C{2};
        signalMeasured = C{3};

        fprintf("Import completed. %d points loaded.\n", length(Latitude));

        %% Visualize measured data
        fprintf("Plotting measured data on map...\n");

        viewer2 = siteviewer("Name", "Measured Data Viewer");

        % Create table and propagation data object for visualization
        tbl_measured = table(Latitude, Longitude, signalMeasured, ...
            'VariableNames', {'Latitude', 'Longitude', 'SignalStrength'});
        pd_measured = propagationData(tbl_measured);

        % Plot measured data on map
        plot(pd_measured, "Colormap", "turbo", "LegendTitle", "Measured RSSI (dBm)");
        contour(pd_measured, "Colormap", "turbo", "LegendTitle", "Measured RSSI (dBm)");

        %% Compute simulated RSSI at measurement points
        fprintf("\nComputing simulated RSSI for each measurement point...\n");

        n = length(Latitude);
        simulatedRSSI = zeros(size(Latitude));

        % Calculate simulated RSSI for each measured point
        for i = 1:n
            % Create receiver at measured location
            rx = rxsite("Latitude", Latitude(i), ...
                "Longitude", Longitude(i), ...
                "AntennaHeight", 1.5);  % Typical height for mobile device

            % Calculate signal strength
            simulatedRSSI(i) = sigstrength(rx, tx, pm);

            % Display progress
            if mod(i, ceil(n/10)) == 0 || i == n
                fprintf("Progress: %d/%d points (%.1f%%)\n", i, n, 100*i/n);
            end
        end

        fprintf("Simulated RSSI computation completed.\n");

        %% Visualize simulated RSSI
        fprintf("Plotting simulated data on map...\n");

        viewer3 = siteviewer("Name", "Simulated Data Viewer");

        % Create table and propagation data object for visualization
        tbl_simulated = table(Latitude, Longitude, simulatedRSSI, ...
            'VariableNames', {'Latitude', 'Longitude', 'SignalStrength'});
        pd_simulated = propagationData(tbl_simulated);

        % Plot simulated data on map
        plot(pd_simulated, "Colormap", "parula", "LegendTitle", "Simulated RSSI (dBm)");
        contour(pd_simulated, "Colormap", "parula", "LegendTitle", "Simulated RSSI (dBm)");

        %% Compute error between measured and simulated values
        fprintf("Computing error between measured and simulated RSSI...\n");

        errorRSSI = signalMeasured - simulatedRSSI;  % Error in dB
        absError = abs(errorRSSI);                   % Absolute error

        % Calculate error statistics
        meanError = mean(errorRSSI);                 % Mean error (bias)
        rmse = sqrt(mean(errorRSSI.^2));             % Root Mean Square Error
        mae = mean(absError);                        % Mean Absolute Error

        % Display error statistics
        fprintf("\n=== ERROR STATISTICS ===\n");
        fprintf("Mean Error (Measured - Simulated): %.2f dB\n", meanError);
        fprintf("RMSE: %.2f dB\n", rmse);
        fprintf("MAE: %.2f dB\n", mae);

        %% Visualize percentage error
        fprintf("Plotting percentage error heatmap...\n");

        % Filter out invalid RSSI values to avoid division issues
        minValidRSSI = -100;
        validIdx = abs(signalMeasured) > 1 & signalMeasured < 0;

        % Calculate percentage error
        percentageError = zeros(size(signalMeasured));
        percentageError(validIdx) = abs(errorRSSI(validIdx) ./ signalMeasured(validIdx)) * 100;

        % Limit maximum error for better visualization
        percentageError = min(percentageError, 100);  % Cap at 100%

        viewer5 = siteviewer("Name", "Percentage Error Viewer");

        % Create table and propagation data object for visualization
        tbl_percent = table(Latitude, Longitude, percentageError, ...
            'VariableNames', {'Latitude', 'Longitude', 'SignalStrength'});
        pd_percent = propagationData(tbl_percent);

        % Plot percentage error on map
        plot(pd_percent, "Colormap", flipud(hot), "LegendTitle", "RSSI % Error");
        contour(pd_percent, "Colormap", flipud(hot), "LegendTitle", "RSSI % Error");

        %% Error histogram
        fprintf("Plotting histogram of RSSI error...\n");

        figure;
        histogram(errorRSSI, 20);  % Create histogram with 20 bins
        title("RSSI Error Distribution (Measured - Simulated)");
        xlabel("Error (dB)");
        ylabel("Frequency");
        grid on;

        fprintf("\n=== SIMULATION COMPLETE ===\n");

    otherwise
        % Handle invalid choice (should not occur due to input validation)
        error('Invalid option! Please run again and choose 1 or 2.');
end
