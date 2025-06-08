% ============================================================
%  PEDRO FERREIRA 2222035 | BRUNO VICENTE 2210709
%  COMUNICACOES MOVEIS 24/25
% ===========================================================

clc; clear all; close all;

tic

fprintf(['\n\n\n',...
    '============================================================\n', ...
    ' PEDRO FERREIRA 2222035 | BRUNO VICENTE 2210709\n', ...
    ' COMUNICACOES MOVEIS 24/25\n', ...
    '============================================================\n\n']);

disp('Choose the scenario for analysis:');
disp('1 - 5G Network in Leiria city');
disp('2 - Raytracing Simulation ESTG');
choice = 0;
while ~ismember(choice, [1, 2])
    choice = input('Enter 1 or 2: ');
    if ~ismember(choice, [1, 2])
        disp('Invalid option! Please enter 1 or 2.');
    end
end

switch choice
    case 1
        fprintf('You choose: 5G Network in Leiria city\n');
        % Coordenadas das Base Stations

        % ZU - Zona Urbana (3600 MHz)
        % ZR - Zona Rural  ( 700 MHz)

        baseStations = [
            39.731823, -8.798338 %ZR 1
            39.736940, -8.799721 %ZU 2 -
            39.739715, -8.802933 %ZU 3 -
            39.737284, -8.810495 %ZU 4 -
            39.739443, -8.819369 %ZU 5 -
            39.732543, -8.824665 %ZR 6
            39.725831, -8.830228 %ZU 7 -
            39.744488, -8.817963 %ZR 8
            39.747375, -8.815152 %ZU 9 -
            39.746381, -8.804517 %ZU 10 -

            %39.747029, -8.809635 %  Castelo Leiria
            %39.739583, -8.871774 %  ESTG
            %39.751326, -8.781058 %  Andrinos
            %39.712350, -8.849653
            ];

        numStations = size(baseStations, 1);


        % Criação da antena nº9 com sectorização e sendo esta de patchs
        fq=0.7e9;
        lambda = physconst('lightspeed') / fq;
        txElement = design(patchMicrostrip, fq);
        txElement.Tilt = 90;
        txElement.TiltAxis = [0 1 0];
        ntxrow = 2;
        ntxcol = 2;
        drow = lambda/2;
        dcol = lambda/2;

        txArray = phased.URA("Size", [ntxrow ntxcol], ...
            "Element", txElement, ...
            "ElementSpacing", [drow dcol]);

        %Plots de diagrams da antena
        az = -180:1:180;
        el = 0;
        az_el = 0;

        % Padrões cortes azimute e elevação
        pat_az = pattern(txArray, fq, az, el);
        pat_el = pattern(txArray, fq, az_el, -90:1:90);

        figure('Position',[100 100 1200 900]);

        % 1. Visualizar layout do array
        subplot(2,2,1);
        pos = txArray.getElementPosition;
        scatter3(pos(1,:), pos(2,:), pos(3,:), 50,'r', 'filled');
        grid on; axis equal;
        xlabel('X (m)');
        ylabel('Y (m)');
        zlabel('Z (m)');
        title('Posição dos Elementos do Array URA');
        view(3);

        % 3. Corte em azimute (elevação = 0°)
        subplot(2,2,2);
        patternAzimuth(txArray, fq);

        % 4. Corte em elevação (azimute = 0°)
        subplot(2,2,3);
        patternElevation(txArray, fq);

        % 2. Padrão 3D do array
        subplot(2,2,4);
        pattern(txArray, fq, -180:180, -90:90, 'Type', 'directivity');
        title('Padrão 3D do Array URA');

        zu_indices = [2,3,4,5,7,9,10];  % Change from [1, 2, 9] to use only valid indices
        sectorized_index = 5;    % Apenas a 5 será setorizada
        azimutes = [0, 120, 240];

        % Define individual heights for each base station (replace single height variable)
        stationHeights = 28 * ones(1, size(baseStations, 1));
        stationHeights(3) = 40;  % Set base station 2 to be higher (40m)

        numStations = size(baseStations,1);

        % Inicializa vetor vazio de txSites
        txSites = txsite.empty;

        siteCount = 1;  % contador global de sites

        % Display antenna choices
        disp('Choose the antenna type:');
        disp('1 - Isotropic');
        disp('2 - Sectorized patch array');
        antennaChoice = 0;
        while ~ismember(antennaChoice, [1, 2])
            antennaChoice = input('Enter 1 or 2: ');
            if ~ismember(antennaChoice, [1, 2])
                disp('Invalid option! Please enter 1 or 2.');
            end
        end

        % Initialize txSites array
        txSites = txsite.empty;
        siteCount = 1;  % contador global de sites

        for i = 1:numStations
            % Define parâmetros de frequência e potência
            if ismember(i, zu_indices)
                freq = 3.6e9;
                power = 5;
            else
                freq = 0.7e9;
                power = 20;
            end

            % Create base station with the appropriate configuration
            if i == sectorized_index
                for sector = 1:3
                    if antennaChoice == 1
                        % Isotropic antenna (don't include Antenna parameter)
                        txSites(siteCount) = txsite( ...
                            'Name', sprintf('Base%d_Sector%d', i, sector), ...
                            'Latitude', baseStations(i,1), ...
                            'Longitude', baseStations(i,2), ...
                            'AntennaHeight', stationHeights(i), ...
                            'TransmitterFrequency', freq, ...
                            'TransmitterPower', power, ...
                            'AntennaAngle', azimutes(sector));
                    else
                        % Sectorized patch array
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
                if antennaChoice == 1
                    % Isotropic antenna (don't include Antenna parameter)
                    txSites(siteCount) = txsite( ...
                        'Name', sprintf('BaseStation%d', i), ...
                        'Latitude', baseStations(i,1), ...
                        'Longitude', baseStations(i,2), ...
                        'AntennaHeight', stationHeights(i), ...
                        'TransmitterFrequency', freq, ...
                        'TransmitterPower', power);
                else
                    % Sectorized patch array
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

        % Modelo de propagação
        modelList = {'close-in', 'longley-rice', 'freespace'};
        fprintf('Available propagation models:\n');
        for idx = 1:numel(modelList)
            fprintf('  %d - %s\n', idx, modelList{idx});
        end
        selection = input('Enter the number of the propagation model: ');
        if isnumeric(selection) && selection >= 1 && selection <= numel(modelList)
            modelName = modelList{selection};
            try
                pm = propagationModel(modelName);
                disp(['Using model: ', modelName]);
                coverage(txSites, ...
                    'SignalStrengths', -110:1:-10, ...
                    'MaxRange', 5000, ...
                    'PropagationModel', pm, ...
                    'Resolution', 20);
            catch ME
                warning("Error applying model '%s': %s", modelName, ME.message);
            end
        else
            disp('No model selected. Script terminated.');
        end

        toc

        %% Geração da grelha de pontos
        latLim = [min([txSites.Latitude])-0.005, max([txSites.Latitude])+0.005];
        lonLim = [min([txSites.Longitude])-0.005, max([txSites.Longitude])+0.005];

        % Resolução da grelha (reduzida para processamento mais rápido)
        gridResolution = 50;
        [latGrid, lonGrid] = meshgrid(linspace(latLim(1), latLim(2), gridResolution), ...
            linspace(lonLim(1), lonLim(2), gridResolution));

        fprintf('Calculating RSSI for grid %dx%d...\n', gridResolution, gridResolution);

        %% Calcular RSSI para todos os txSites
        rssiMatrix = zeros([size(latGrid), numel(txSites)]);
        for i = 1:numel(txSites)
            fprintf('Calculating RSSI for base station %d of %d...\n', i, numel(txSites));
            rx = rxsite('Latitude', latGrid(:), 'Longitude', lonGrid(:));
            sig = sigstrength(rx, txSites(i), pm);
            rssiMatrix(:,:,i) = reshape(sig, size(latGrid));
        end

        %% Estudo de Best Server e Interferencia Co-Canal

        % Cálculo do melhor servidor
        [bestRSSI, bestServerIndex] = max(rssiMatrix, [], 3);

        % Estatísticas por estação
        fprintf('Best Server Statistics:\n');
        for i = 1:numel(txSites)
            coverage_area = sum(bestServerIndex(:) == i) / numel(bestServerIndex) * 100;
            fprintf('Station %d: %.2f%% of the area\n', i, coverage_area);
        end

        % Margem de dominância (RSSI1 - RSSI2)
        sortedRSSI = sort(rssiMatrix, 3, 'descend');
        margin = sortedRSSI(:,:,1) - sortedRSSI(:,:,2);
        % Limite para zona de fronteira (margem pequena)
        threshold = 3; % dB
        % Create separate figure windows for each plot

        % Figure 1 – Mapa Best Server
        figure('Name', 'Mapa Best Server', 'Position', [100 100 700 600]);
        pcolor(lonGrid, latGrid, bestServerIndex);
        shading flat;
        colormap(parula(numel(txSites)));  % Changed from jet to parula for better color distinction
        colorbar('Ticks', 1:numel(txSites), 'TickLabels', {txSites.Name});
        axis xy;
        hold on;

        for i = 1:numel(txSites)
            plot(txSites(i).Longitude, txSites(i).Latitude, 'r^', 'MarkerSize', 10, 'LineWidth', 2);
            text(txSites(i).Longitude, txSites(i).Latitude, [' BS', num2str(i)], 'Color', 'black', 'FontWeight', 'bold');
        end
        title('Mapa Best Server', 'FontWeight', 'bold');
        xlabel('Longitude'); ylabel('Latitude');
        hold off;

        % Figure 2 – Margem de Dominância
        figure('Name', 'Margem de Dominância', 'Position', [150 150 700 600]);
        pcolor(lonGrid, latGrid, margin);
        shading flat;
        colormap(hot);
        colorbar;
        title('Margem de Domínio (RSSI_1 - RSSI_2) [dB]', 'FontWeight', 'bold');
        xlabel('Longitude'); ylabel('Latitude');
        hold on;

        borderZone = margin < threshold;
        contour(lonGrid, latGrid, borderZone, [1 1], 'LineColor', 'b', 'LineWidth', 1.5);
        legend('Margem < 3 dB (Zona de Instabilidade)', 'Location', 'northeast');
        hold off;

        % Figure 3 – Mapa SINR Estimado
        figure('Name', 'Mapa SINR Estimado', 'Position', [200 200 700 600]);
        signal = sortedRSSI(:,:,1);
        total_power = sum(10.^(rssiMatrix/10), 3);
        interference = 10*log10(total_power - 10.^(signal/10) + 1e-10);
        SINR = signal - interference;

        pcolor(lonGrid, latGrid, SINR);
        shading flat;
        colormap(parula);
        colorbar;
        title('Mapa de SINR Estimado [dB]', 'FontWeight', 'bold');
        xlabel('Longitude'); ylabel('Latitude');

        % Figure 4 – Distribuição de Best Server
        figure('Name', 'Distribuição de Best Server', 'Position', [250 250 700 600]);
        counts = histcounts(bestServerIndex, 1:(numel(txSites)+1));
        b = bar(counts);
        b.FaceColor = 'flat';  % Enable coloring each bar individually

        % Create a colormap for the bars matching the best server colors
        barColors = parula(numel(txSites));
        for i = 1:numel(txSites)
            b.CData(i,:) = barColors(i,:);
        end

        xticks(1:numel(txSites));
        xticklabels({txSites.Name});
        xlabel('Estação Base');
        ylabel('N.º de pontos como Best Server');
        title('Distribuição de Best Server', 'FontWeight', 'bold');

        % Add a grid to the bar chart for better readability
        grid on;


        %% Estudo do Handover
        fprintf('\n--- Handover Study ---\n');
        fprintf('Chosen path type: Random Path \n');
        % TRAJETO ALEATÓRIO
        pathName = 'Random Path';
        numPoints = 50;

        minLat = 39.78490706138443; maxLat = 39.72129535847994;
        minLon = -8.859864711773678; maxLon = -8.780557161553446;

        latitudes = minLat + (maxLat - minLat) * rand(1, numPoints);
        longitudes = minLon + (maxLon - minLon) * rand(1, numPoints);

        % CÁLCULO DO RSSI
        rssiMatrix = zeros(numPoints, numel(txSites));
        fprintf('Calculating RSSI along the path...\n');
        for i = 1:numel(txSites)
            rx = rxsite('Latitude', latitudes, 'Longitude', longitudes);
            sig = sigstrength(rx, txSites(i), pm);
            sig = sig(:);
            rssiMatrix(:, i) = sig;
        end

        [bestRSSI, bestServerIndex] = max(rssiMatrix, [], 2);

        % Figura única com 3 subplots
        figure('Name','Estudo de Handover - Trajeto, RSSI e Best Server', 'Position',[100 100 1200 800]);

        % Subplot 1 – Mapa com o trajeto e estações base (com cores únicas)
        subplot(3,1,1);
        geoplot(latitudes, longitudes, 'b-', 'LineWidth', 2); hold on;

        colors = lines(numel(txSites));  % Gera cores distintas automaticamente

        for i = 1:numel(txSites)
            geoplot(txSites(i).Latitude, txSites(i).Longitude, ...
                '^', 'MarkerSize', 10, ...
                'MarkerEdgeColor', 'k', ...
                'MarkerFaceColor', colors(i,:), ...
                'DisplayName', txSites(i).Name);
        end

        title(['[Mapa] Trajeto RX e Estações Base - ', pathName]);
        legend('Trajeto RX', txSites.Name, 'Location','best');
        grid on;

        % Subplot 3 – Estação ativa (Best Server)
        subplot(3,1,2);
        plot(1:numPoints, bestServerIndex, 'k', 'LineWidth', 1.5);
        xlabel('Ponto no trajeto');
        ylabel('Estação Base Ativa');
        yticks(1:numel(txSites));
        yticklabels({txSites.Name});
        title('[Gráfico] Estação Servidora (Best Server) ao Longo do Trajeto');
        grid on;

        % Subplot 2 – RSSI ao longo do percurso
        subplot(3,1,3);
        plot(1:numPoints, rssiMatrix, 'LineWidth', 1.2);
        xlabel('Ponto no trajeto');
        ylabel('RSSI (dBm)');
        title('[Gráfico] RSSI Recebido por Estação Base');
        legend({txSites.Name}, 'Location','best');
        grid on;

        % HANDOVER - identificação e impressão
        handoverPoints = find(diff(bestServerIndex) ~= 0);
        fprintf('\n--- Handovers Detected ---\n');
        for i = 1:length(handoverPoints)
            idx = handoverPoints(i);
            fprintf('Handover at point %d: %s → %s\n', ...
                idx+1, txSites(bestServerIndex(idx)).Name, txSites(bestServerIndex(idx+1)).Name);
        end

        %% Estudo de capacidade

        fprintf('\n--- Capacity Analysis with 11940 Users ---\n');

        % Número total de usuários
        nUsers = 128642;

        % Definir os pesos de distribuição de usuários (maior peso nas áreas urbanas - antenas 1, 2 e 9)
        % Vamos calcular as posições de todas as antenas primeiro
        txPositions = zeros(numel(txSites), 2);
        for i = 1:numel(txSites)
            txPositions(i, 1) = txSites(i).Latitude;
            txPositions(i, 2) = txSites(i).Longitude;
        end

        % Gerar usuários com distribuição não uniforme, concentrando nas áreas urbanas
        userLat = zeros(nUsers, 1);
        userLon = zeros(nUsers, 1);

        % Definir proporções aproximadas: 75% dos usuários nas áreas urbanas (antenas 1, 2 e 9)
        % e 25% distribuídos entre as demais antenas
        urbanUsers = round(0.75 * nUsers);
        ruralUsers = nUsers - urbanUsers;

        % Definir raios de influência (menores para áreas urbanas, maiores para rurais)
        urbanRadius = 0.015; % Raio menor para áreas urbanas (graus de lat/lon)
        ruralRadius = 0.06; % Raio maior para áreas rurais

        % Identificar quais são as antenas urbanas
        urbanAntennas = [1, 2];  % Adjust this to include only valid indices from your txPositions array
        ruralAntennas = setdiff(1:numel(txSites), urbanAntennas);

        % Distribuir usuários urbanos entre as antenas 1, 2 e 9
        % Pesos para distribuição desigual entre as antenas urbanas
        urbanWeights = [0.4, 0.6]; % Updated weights for just 2 antennas
        urbanDistribution = round(urbanUsers * urbanWeights);
        % Ajustar para garantir que soma seja exatamente urbanUsers
        urbanDistribution(end) = urbanUsers - sum(urbanDistribution(1:end-1));

        currentUserIndex = 1;

        for i = 1:length(urbanAntennas)
            antennaIndex = urbanAntennas(i);
            centerLat = txPositions(antennaIndex, 1);
            centerLon = txPositions(antennaIndex, 2);

            % Gerar usuários em torno desta antena urbana com distribuição gaussiana
            usersForThisAntenna = urbanDistribution(i);

            % Gerar coordenadas com distribuição gaussiana ao redor da antena
            r = urbanRadius * sqrt(rand(usersForThisAntenna, 1)); % Distribuição uniforme em área
            theta = 2 * pi * rand(usersForThisAntenna, 1);

            userLat(currentUserIndex:currentUserIndex+usersForThisAntenna-1) = centerLat + r .* cos(theta);
            userLon(currentUserIndex:currentUserIndex+usersForThisAntenna-1) = centerLon + r .* sin(theta);

            currentUserIndex = currentUserIndex + usersForThisAntenna;
        end

        % Distribuir usuários rurais entre as demais antenas
        ruralUsersPerAntenna = round(ruralUsers / length(ruralAntennas));

        for i = 1:length(ruralAntennas)
            antennaIndex = ruralAntennas(i);
            centerLat = txPositions(antennaIndex, 1);
            centerLon = txPositions(antennaIndex, 2);

            % Gerar usuários em torno desta antena rural
            usersForThisAntenna = ruralUsersPerAntenna;
            if i == length(ruralAntennas)
                % Ajustar o último grupo para garantir que o total seja exato
                usersForThisAntenna = ruralUsers - (i-1) * ruralUsersPerAntenna;
            end

            % Gerar coordenadas com distribuição gaussiana ao redor da antena
            r = ruralRadius * sqrt(rand(usersForThisAntenna, 1));
            theta = 2 * pi * rand(usersForThisAntenna, 1);

            userLat(currentUserIndex:currentUserIndex+usersForThisAntenna-1) = centerLat + r .* cos(theta);
            userLon(currentUserIndex:currentUserIndex+usersForThisAntenna-1) = centerLon + r .* sin(theta);

            currentUserIndex = currentUserIndex + usersForThisAntenna;
        end

        % Garantir que as coordenadas estão dentro dos limites do mapa
        userLat = max(min(userLat, latLim(2)), latLim(1));
        userLon = max(min(userLon, lonLim(2)), lonLim(1));

        % Criar objeto receptor para todos os usuários uma única vez
        rx = rxsite('Latitude', userLat, 'Longitude', userLon);

        % Calcular RSSI para todos usuários e todas estações base
        userRSSI = zeros(nUsers, numel(txSites));
        for i = 1:numel(txSites)
            userRSSI(:, i) = sigstrength(rx, txSites(i), pm);
        end

        % Definir larguras de banda (Hz) antes de calcular capacidade
        bandwidths = 100e6 * ones(1, numel(txSites));
        bandwidths([1, 2]) = 100e6;  % Only reference existing stations

        % Encontrar a melhor estação para cada usuário (maior RSSI)
        [bestRSSI, userAssignedBS] = max(userRSSI, [], 2);

        % Definir ruído (em dBm e converter para mW)
        noise_dBm = -100;
        noise_mW = 10^(noise_dBm / 10);

        % Converter RSSI para potência linear (mW)
        userPower_mW = 10.^(userRSSI / 10);

        % Calcular SINR para todos usuários e todas BS
        userSINR_all = zeros(nUsers, numel(txSites));
        for u = 1:nUsers
            for bs = 1:numel(txSites)
                signal = userPower_mW(u, bs);
                interference = sum(userPower_mW(u, :)) - signal;
                userSINR_all(u, bs) = signal / (interference + noise_mW);
            end
        end

        % SINR do BS associado
        userSINR = zeros(nUsers, 1);
        for u = 1:nUsers
            userSINR(u) = userSINR_all(u, userAssignedBS(u));
        end

        % Calcular capacidade individual por usuário com fórmula de Shannon (bps)
        capacityPerUser_shannon = bandwidths(userAssignedBS)' .* log2(1 + userSINR);

        % Calcular capacidade total por estação base
        capacityPerStation = zeros(1, numel(txSites));
        for i = 1:numel(txSites)
            usersInBS = (userAssignedBS == i);
            capacityPerStation(i) = sum(capacityPerUser_shannon(usersInBS));
        end

        % Calcular capacidade média por usuário por estação base
        capacityPerUser = zeros(1, numel(txSites));
        for i = 1:numel(txSites)
            usersInBS = (userAssignedBS == i);
            if sum(usersInBS) > 0
                capacityPerUser(i) = capacityPerStation(i) / sum(usersInBS);
            end
        end


        % Contagem de usuários por estação base
        capacityCount = histcounts(userAssignedBS, 1:(numel(txSites)+1));

        % Visualização dos resultados
        figure('Name', 'Capacity with 11940 Users');

        % Subplot 1: Distribuição de usuários no mapa, colorido pelo SINR em dB
        %subplot(1, 2, 1);
        %scatter(userLon, userLat, 10, 10*log10(userSINR), 'filled');
        %colormap(jet);
        %colorbar;
        %caxis([0 30]);
        %hold on;

        %% Plot antenas urbanas e rurais com marcadores e cores diferentes
        %urbanAntennas = [1, 2, 9];
        %for i = 1:numel(txSites)
        %    if ismember(i, urbanAntennas)
        %        plot(txSites(i).Longitude, txSites(i).Latitude, 'rs', 'MarkerSize', 12, 'MarkerFaceColor', 'red');
        %    else
        %        plot(txSites(i).Longitude, txSites(i).Latitude, 'k^', 'MarkerSize', 10, 'MarkerFaceColor', 'green');
        %    end
        %    text(txSites(i).Longitude, txSites(i).Latitude, [' BS' num2str(i)], 'FontSize', 9, 'FontWeight', 'bold');
        %end
        %xlim(lonLim);
        %ylim(latLim);
        %title('Distribuição de Usuários (cor: SINR dB)');
        %xlabel('Longitude');
        %ylabel('Latitude');
        %legend('Usuários', 'Estações Urbanas', 'Estações Rurais', 'Location', 'best');

        % Subplot 2: Capacidade média por usuário por estação (Mbps)
        subplot(1, 2, 2);
        bar(capacityPerUser / 1e6);
        xticks(1:numel(txSites));
        xticklabels({txSites.Name});
        title('Capacidade Média por Usuário');
        xlabel('Estação Base');
        ylabel('Capacidade (Mbps)');

        % Figura extra: Distribuição de usuários por estação base
        figure('Name', 'Distribuição de Usuários por Estação');
        bar(capacityCount);
        xticks(1:numel(txSites));
        xticklabels({txSites.Name});
        title('Distribuição de Usuários por Estação Base');
        xlabel('Estação Base');
        ylabel('Número de Usuários');
        hold on;
        urbanAntennas = [1, 2];  % Define only valid indices
        bar(urbanAntennas, capacityCount(urbanAntennas), 'r');
        legend('Estações Rurais', 'Estações Urbanas');

        % Impressão dos resultados
        fprintf('Capacity Statistics with %d users:\n', nUsers);
        totalUrbanUsers = 0;
        totalRuralUsers = 0;

        for i = 1:numel(txSites)
            if ismember(i, urbanAntennas)
                stationType = 'Urban';
                totalUrbanUsers = totalUrbanUsers + capacityCount(i);
            else
                stationType = 'Rural';
                totalRuralUsers = totalRuralUsers + capacityCount(i);
            end

            fprintf('Station %d (%s): %d users (%.1f%%) - Type: %s\n', ...
                i, txSites(i).Name, capacityCount(i), 100*capacityCount(i)/nUsers, stationType);

            if capacityCount(i) > 0
                fprintf('  Total capacity: %.2f Mbps\n', capacityPerStation(i)/1e6);
                fprintf('  Average capacity per user: %.2f Mbps\n', capacityPerUser(i)/1e6);
            else
                fprintf('  [No assigned users]\n');
            end
        end

        fprintf('\nFinal Summary:\n');
        fprintf('Total users: %d\n', nUsers);
        fprintf('Users in urban areas: %d (%.1f%%)\n', totalUrbanUsers, 100*totalUrbanUsers/nUsers);
        % Ajusta o separador se necessário (ex: ',' ou ';')lRuralUsers, 100*totalRuralUsers/nUsers);
        %data = readmatrix("Bernardo_Joao_CM.txt");  % Make sure the file is in the same directory

    case 2
        printf("=== STARTING SIMULATION ===\n");

        %% === INITIAL CONFIGURATION ===
        fprintf("Setting up base station and transmitter parameters...\n");

        % Coordinates of the base station
        latitude = 39.735250;
        longitude = -8.820638;

        % Transmitter parameters
        txHeight = 3;                  % Antenna height in meters
        Ptx_dBm = 0;                   % Transmission power in dBm
        Ptx_W = 10^((Ptx_dBm - 30)/10); % Convert dBm to Watts
        frequency = 3e9;               % Frequency: 3 GHz

        % Create transmitter site
        tx = txsite("Name", "Base Station", ...
            "Latitude", latitude, ...
            "Longitude", longitude, ...
            "AntennaHeight", txHeight, ...
            "TransmitterPower", Ptx_W, ...
            "TransmitterFrequency", frequency);

        %% === SETUP RAY TRACING MODEL ===
        fprintf("Creating ray tracing propagation model...\n");

        pm = propagationModel("raytracing", ...
            "Method", "sbr", ...
            "MaxNumReflections", 2, ...
            "SurfaceMaterial", "concrete");

        viewer1 = siteviewer("Name", "Ray Tracing Viewer", ...
            Buildings="estg.osm", ...
            Terrain="gmted2010");

        fprintf("Plotting simulated coverage using ray tracing...\n");

        coverage(tx, ...
            "PropagationModel", pm, ...
            "MaxRange", 2000, ...
            "Resolution", 10, ...
            "SignalStrengths", [-120 -100 -90 -80 -70 -60 -50]);

        %% === IMPORT MEASURED DATA ===
        fprintf("\nImporting measured data from file...\n");

        fid = fopen('gupo1.txt');
        C = textscan(fid, '%f%f%f','Delimiter',',');
        fclose(fid);

        Latitude = C{1};
        Longitude = C{2};
        signalMeasured = C{3};

        fprintf("Import completed. %d points loaded.\n", length(Latitude));

        %% === VISUALIZE MEASURED DATA ===
        fprintf("Plotting measured data on map...\n");

        viewer2 = siteviewer("Name", "Measured Data Viewer");

        tbl_measured = table(Latitude, Longitude, signalMeasured, ...
            'VariableNames', {'Latitude', 'Longitude', 'SignalStrength'});
        pd_measured = propagationData(tbl_measured);

        plot(pd_measured, "Colormap", "turbo", "LegendTitle", "Measured RSSI (dBm)");
        contour(pd_measured, "Colormap", "turbo", "LegendTitle", "Measured RSSI (dBm)");

        %% === COMPUTE SIMULATED RSSI ===
        fprintf("\nComputing simulated RSSI for each measurement point...\n");

        n = length(Latitude);
        simulatedRSSI = zeros(size(Latitude));

        for i = 1:n
            rx = rxsite("Latitude", Latitude(i), ...
                "Longitude", Longitude(i), ...
                "AntennaHeight", 1.5);

            simulatedRSSI(i) = sigstrength(rx, tx, pm);

            if mod(i, ceil(n/10)) == 0 || i == n
                fprintf("Progress: %d/%d points (%.1f%%)\n", i, n, 100*i/n);
            end
        end

        fprintf("Simulated RSSI computation completed.\n");

        %% === VISUALIZE SIMULATED RSSI ===
        fprintf("Plotting simulated data on map...\n");

        viewer3 = siteviewer("Name", "Simulated Data Viewer");

        tbl_simulated = table(Latitude, Longitude, simulatedRSSI, ...
            'VariableNames', {'Latitude', 'Longitude', 'SignalStrength'});
        pd_simulated = propagationData(tbl_simulated);

        plot(pd_simulated, "Colormap", "parula", "LegendTitle", "Simulated RSSI (dBm)");
        contour(pd_simulated, "Colormap", "parula", "LegendTitle", "Simulated RSSI (dBm)");

        %% === ERROR ANALYSIS ===
        fprintf("Computing error between measured and simulated RSSI...\n");

        errorRSSI = signalMeasured - simulatedRSSI;
        absError = abs(errorRSSI);

        meanError = mean(errorRSSI);
        rmse = sqrt(mean(errorRSSI.^2));
        mae = mean(absError);

        fprintf("\n=== ERROR STATISTICS ===\n");
        fprintf("Mean Error (Measured - Simulated): %.2f dB\n", meanError);
        fprintf("RMSE: %.2f dB\n", rmse);
        fprintf("MAE: %.2f dB\n", mae);

        % %% === PLOT RSSI ERROR HEATMAP ===
        % fprintf("Plotting heatmap of RSSI error...\n");
        %
        % viewer4 = siteviewer("Name", "RSSI Error Viewer");
        %
        % tbl_error = table(Latitude, Longitude, errorRSSI, ...
        %     'VariableNames', {'Latitude', 'Longitude', 'SignalStrength'});
        % pd_error = propagationData(tbl_error);
        %
        % plot(pd_error, "Colormap", "jet", "LegendTitle", "RSSI Error (dB)");
        % contour(pd_error, "Colormap", "jet", "LegendTitle", "RSSI Error (dB)");

        %% === PERCENTAGE ERROR MAP ===
        fprintf("Plotting percentage error heatmap...\n");

        % Avoid division by zero or values too close to 0
        minValidRSSI = -100; % filter unrealistic RSSI values
        validIdx = abs(signalMeasured) > 1 & signalMeasured < 0;

        % Preallocate and compute percentage error
        percentageError = zeros(size(signalMeasured));
        percentageError(validIdx) = abs(errorRSSI(validIdx) ./ signalMeasured(validIdx)) * 100;

        % Clamp large errors to improve visualization
        percentageError = min(percentageError, 100); % cap at 100%

        viewer5 = siteviewer("Name", "Percentage Error Viewer");

        tbl_percent = table(Latitude, Longitude, percentageError, ...
            'VariableNames', {'Latitude', 'Longitude', 'SignalStrength'});
        pd_percent = propagationData(tbl_percent);

        % Use 'hot' colormap but flip it so red is for high values (100% error)
        plot(pd_percent, "Colormap", flipud(hot), "LegendTitle", "RSSI % Error");
        contour(pd_percent, "Colormap", flipud(hot), "LegendTitle", "RSSI % Error");


        %% === ERROR HISTOGRAM ===
        fprintf("Plotting histogram of RSSI error...\n");

        figure;
        histogram(errorRSSI, 20);
        title("RSSI Error Distribution (Measured - Simulated)");
        xlabel("Error (dB)");
        ylabel("Frequency");
        grid on;

        fprintf("\n=== SIMULATION COMPLETE ===\n");

    otherwise
        error('Invalid option! Please run again and choose 1 or 2.');
end