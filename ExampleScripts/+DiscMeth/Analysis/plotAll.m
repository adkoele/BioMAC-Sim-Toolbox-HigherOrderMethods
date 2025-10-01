%======================================================================
%> @file    +DiscMeth/plotAll.m
%> @brief   Function to generate summary plots for all simulations.
%>
%> @details
%> This function creates plots from dynamic residuals and 
%> other simulation results. To provide visual 
%> comparisons across participants, movements, trials, and 
%> discretization methods. 
%>
%> @param   path2repo        String: Path to the BioMAC-Sim-Toolbox repository
%> @param   dynResTable      Table: Table containing dynamic residuals and 
%>                           other simulation results
%> @param   column_width_in  Double: Column width of the figure in inches
%> @param   height_in        Double: Height of the figure in inches
%> @param   nPar             Integer: Number of participants
%> @param   mov              String: Movement type(s)
%> @param   tra              String: Trial name(s)
%> @param   FontName         String: Font name for all text in the plot
%> @param   FontSize         Double: Font size
%> @param   LineWidth        Double: Line width for plotted lines
%> @param   dim              String: Dimension of the simulation ('2D' or '3D')
%>
%> @retval  None             The function generates plots but does not return 
%>                           output variables.
%>
%> @author  Alexander Weiss 
%> @date    August, 2025
%======================================================================



function plotAll(path2repo, dynResTable, column_width_in, height_in, nPar, mov, tra, FontName, FontSize, LineWidth, dim)
    % Convert first column (method names) to categorical
    methods = unique(dynResTable{:,'Method'}, 'stable');  % Extract unique method names in the order they appear
    nMethods = length(methods);
    
    % Extract unique nNodes values
    nNodes = unique(dynResTable{:,'nNodes'}, 'sorted');
    
    % Initialize matrices for errors
    errorT = nan(nPar-1,nMethods, length(nNodes));
    errorR = nan(nPar-1,nMethods, length(nNodes));
    errorS = nan(nPar-1,nMethods, length(nNodes));
    errorA = nan(nPar-1,nMethods, length(nNodes));
    serrorR = nan(nPar-1,nMethods, length(nNodes));
    wTime = nan(nPar-1,nMethods, length(nNodes));
    cTime = nan(nPar-1,nMethods, length(nNodes));
    
    %[~,uniqueIdx] = unique(dynResTable(:,{'Method', 'nNodes', 'Participant'}), 'rows');
    %dynResTable = dynResTable(uniqueIdx, :);
    
    for i = 1:nMethods
        for j = 1:length(nNodes)
            % Find the row corresponding to method and nNodes
            idx = strcmp(dynResTable{:,'Method'}, methods{i}) & dynResTable{:,4} == nNodes(j) ...
                & strcmp(dynResTable{:,'Movement'}, mov) &  strcmp(dynResTable{:,'Tracking'}, tra) ...
                & (dynResTable.Converged == 1 | strcmp(dynResTable.Method, 'ME'));        % Number of rows found
            nRows = sum(idx);
            if nRows > 0
                maxFill = min(nRows, nPar - 1); % Don't overflow preallocated arrays            % Assign values to first maxFill rows
                errorT(1:maxFill,i,j)   = dynResTable{idx,'IntErrorT'}(1:maxFill);
                errorR(1:maxFill,i,j)   = dynResTable{idx,'IntErrorR'}(1:maxFill);
                errorA(1:maxFill,i,j)   = dynResTable{idx,'IntErrorA'}(1:maxFill);
                errorS(1:maxFill,i,j)   = dynResTable{idx,'IntErrorS'}(1:maxFill);
                serrorR(1:maxFill,i,j)  = dynResTable{idx,'IntSecErrorR'}(1:maxFill);
                serrorT(1:maxFill,i,j)  = dynResTable{idx,'IntSecErrorT'}(1:maxFill);
                wTime(1:maxFill,i,j)    = dynResTable{idx,'WallTime'}(1:maxFill);
                cTime(1:maxFill,i,j)    = dynResTable{idx,'CpuTime'}(1:maxFill);
                traData(1:maxFill,i,j)  = dynResTable{idx,'tracking'}(1:maxFill);
                deltaT(1:maxFill,i,j)  = dynResTable{idx,'deltaT'}(1:maxFill);
                dynThreshold(1:maxFill,i,j)  = dynResTable{idx,'IntThreshold'}(1:maxFill);
            end
        end
    end
    
    
    % Define colors for each method
    colors = lines(nMethods);


     %% Length Errors
    figure();
    
   if strcmp(dim, '2D')
        % Set figure properties
        set(gcf, 'Units', 'Inches', 'Position', [1, 1, column_width_in/2, height_in/3]);
        set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [column_width_in/2, height_in/3]);
    else
        set(gcf, 'Units', 'Inches', 'Position', [1, 1, column_width_in/2, height_in/2]);
        set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [column_width_in/2, height_in/2]);
    end
    
    % First subplot: ErrorT
    %subplot(2,1,1);
    hold on;
    for i = 1:nMethods
        % Extract values for the current method: size [n, m]
        values = squeeze(errorS(:, i, :)); % size: [nNodes x mSamples]
    
        % Compute mean and standard deviation across samples (dim = 2)
        meanVals = nanmean(values, 1); % size: [nNodes x 1]
        stdVals  = nanstd(values, 0, 1); % standard deviation
        % Plot mean with error bars
        errorbar(nNodes, meanVals, stdVals, '-o', ...
            'Color', colors(i,:), 'LineWidth', LineWidth, ...
            'DisplayName', methods{i}, ...
            'CapSize', 4, 'MarkerSize', 6);%
    end
    
    hold off;
    xlabel('Number of nodes');
    ylabel('Dynamic error E^{1}');
    %title('Dynamic Error vs. nNodes');
    set(gca, 'YScale', 'log');
    set(gca,'FontName', FontName);
    set(gca,'FontSize', FontSize);
    set(gca,'LineWidth', round(LineWidth/2));
    % ylim([10e-8, 1]);
    
    %set(gca,'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', FontSize*2/3, 'Orientation','vertical');
    grid on;
    
    if size(nNodes,1) > 1
        print(gcf, fullfile([path2repo, filesep,'SErrvsNodes_', mov,'_' tra]), '-depsc', '-r300');
        print(gcf, fullfile([path2repo, filesep,'SErrvsNodes_', mov,'_' tra]), '-dpng', '-r300');
    end

     %% Acivation Errors
    figure();
    
    if strcmp(dim, '2D')
        % Set figure properties
        set(gcf, 'Units', 'Inches', 'Position', [1, 1, column_width_in/2, height_in/3]);
        set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [column_width_in/2, height_in/3]);
    else
        set(gcf, 'Units', 'Inches', 'Position', [1, 1, column_width_in/2, height_in/2]);
        set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [column_width_in/2, height_in/2]);
    end
    
    % First subplot: ErrorT
    %subplot(2,1,1);
    hold on;
    for i = 1:nMethods
        % Extract values for the current method: size [n, m]
        values = squeeze(errorA(:, i, :)); % size: [nNodes x mSamples]
    
        % Compute mean and standard deviation across samples (dim = 2)
        meanVals = nanmean(values, 1); % size: [nNodes x 1]
        stdVals  = nanstd(values, 0, 1); % standard deviation
        % Plot mean with error bars
        errorbar(nNodes, meanVals, stdVals, '-o', ...
            'Color', colors(i,:), 'LineWidth', LineWidth, ...
            'DisplayName', methods{i}, ...
            'CapSize', 4, 'MarkerSize', 6);%
    end
    
    hold off;
    xlabel('Number of nodes');
    ylabel('Dynamic error E^{1}');
    %title('Dynamic Error vs. nNodes');
    set(gca, 'YScale', 'log');
    set(gca,'FontName', FontName);
    set(gca,'FontSize', FontSize);
    set(gca,'LineWidth', round(LineWidth/2));
    % ylim([10e-8, 1]);
    
    %set(gca,'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', FontSize*2/3, 'Orientation','vertical');
    grid on;
    
    if size(nNodes,1) > 1
        print(gcf, fullfile([path2repo, filesep,'AErrvsNodes_', mov,'_' tra]), '-depsc', '-r300');
        print(gcf, fullfile([path2repo, filesep,'AErrvsNodes_', mov,'_' tra]), '-dpng', '-r300');
    end
    
    %% Translation Errors
    figure();
    
    
    if strcmp(dim, '2D')
        % Set figure properties
        set(gcf, 'Units', 'Inches', 'Position', [1, 1, column_width_in/2, height_in/3]);
        set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [column_width_in/2, height_in/3]);
    else
        set(gcf, 'Units', 'Inches', 'Position', [1, 1, column_width_in/2, height_in/2]);
        set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [column_width_in/2, height_in/2]);
    end
    
    % First subplot: ErrorT
    %subplot(2,1,1);
    hold on;
    for i = 1:nMethods
        % Extract values for the current method: size [n, m]
        values = squeeze(errorT(:, i, :)); % size: [nNodes x mSamples]
    
        % Compute mean and standard deviation across samples (dim = 2)
        meanVals = nanmean(values, 1); % size: [nNodes x 1]
        stdVals  = nanstd(values, 0, 1); % standard deviation
        % Plot mean with error bars
        errorbar(nNodes, meanVals, stdVals, '-o', ...
            'Color', colors(i,:), 'LineWidth', LineWidth, ...
            'DisplayName', methods{i}, ...
            'CapSize', 4, 'MarkerSize', 6);%
    end
    
    hold off;
    xlabel('Number of nodes');
    ylabel('Dynamic error E^{1} (m)');
    %title('Dynamic Error vs. nNodes');
    set(gca, 'YScale', 'log');
    set(gca,'FontName', FontName);
    set(gca,'FontSize', FontSize);
    set(gca,'LineWidth', round(LineWidth/2));
    % ylim([10e-8, 1]);
    
    %set(gca,'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', FontSize*2/3, 'Orientation','vertical');
    grid on;
    
    if size(nNodes,1) > 1
        print(gcf, fullfile([path2repo, filesep,'TErrvsNodes_', mov,'_' tra]), '-depsc', '-r300');
        print(gcf, fullfile([path2repo, filesep,'TErrvsNodes_', mov,'_' tra]), '-dpng', '-r300');
    end
    
    
    figure();
    
    if strcmp(dim, '2D')
        % Set figure properties
        set(gcf, 'Units', 'Inches', 'Position', [1, 1, column_width_in/2, height_in/3]);
        set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [column_width_in/2, height_in/3]);
    else
        set(gcf, 'Units', 'Inches', 'Position', [1, 1, column_width_in/2, height_in/2]);
        set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [column_width_in/2, height_in/2]);
    end
    
    
    % First subplot: ErrorT
    %subplot(2,1,1);
    hold on;
    for i = 1:nMethods
        % Extract values for the current method: size [n, m]
        values = squeeze(serrorT(:, i, :)); % size: [nNodes x mSamples]
    
        % Compute mean and standard deviation across samples (dim = 2)
        meanVals = nanmean(values, 1); % size: [nNodes x 1]
        stdVals  = nanstd(values, 0, 1); % standard deviation
        % Alternatively, for standard error:
        % stdVals = std(values, 0, 2) ./ sqrt(size(values,2));
    
        % Plot mean with error bars
        errorbar(nNodes, meanVals, stdVals, '-o', ...
            'Color', colors(i,:), 'LineWidth', LineWidth, ...
            'DisplayName', methods{i}, ...
            'CapSize', 4, 'MarkerSize', 6);
    end
    
    hold off;
    xlabel('Number of nodes');
    ylabel('Dynamic error E^{2} (m/s)');
    %title('Dynamic Error vs. nNodes');
    set(gca, 'YScale', 'log');
    set(gca,'FontName', FontName);
    set(gca,'FontSize', FontSize);
    set(gca,'LineWidth', round(LineWidth/2));
    %set(gca,'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', FontSize*2/3, 'Orientation','vertical');    grid on;
    
    if size(nNodes,1) > 1
        print(gcf, fullfile([path2repo, filesep,'STErrvsNodes_', mov,'_' tra]), '-depsc', '-r300');
        print(gcf, fullfile([path2repo, filesep,'STErrvsNodes_', mov,'_' tra]), '-dpng', '-r300');
    end
    
    
    
    figure();
    
     if strcmp(dim, '2D')
        % Set figure properties
        set(gcf, 'Units', 'Inches', 'Position', [1, 1, column_width_in/2, height_in/3]);
        set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [column_width_in/2, height_in/3]);
    else
        set(gcf, 'Units', 'Inches', 'Position', [1, 1, column_width_in/2, height_in/2]);
        set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [column_width_in/2, height_in/2]);
    end
    
    % First subplot: ErrorT
    %subplot(2,1,1);
    hold on;
    for i = 1:nMethods
        % Extract values for the current method: size [n, m]
        values = squeeze(errorR(:, i, :)); % size: [nNodes x mSamples]
        % Compute mean and standard deviation across samples (dim = 2)
        meanVals = nanmean(values, 1); % size: [nNodes x 1]
        stdVals  = nanstd(values, 0, 1); % standard deviation
        % Plot mean with error bars
        errorbar(nNodes, meanVals, stdVals, '-o', ...
            'Color', colors(i,:), 'LineWidth', LineWidth, ...
            'DisplayName',  methods{i}, ...
            'CapSize', 4, 'MarkerSize', 6);%
    end
    
    hold off;
    xlabel('Number of nodes');
    ylabel('Dynamic error E^{1} (rad)');
    %title('Dynamic Error vs. nNodes');
    set(gca, 'YScale', 'log');
    set(gca,'FontName', FontName);
    set(gca,'FontSize', FontSize);
    set(gca,'LineWidth', round(LineWidth/2));
    
    %set(gca,'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', FontSize*2/3, 'Orientation','vertical');    grid on;
    
    if size(nNodes,1) > 1
        print(gcf, fullfile([path2repo, filesep,'RErrvsNodes_', mov,'_' tra]), '-depsc', '-r300');
        print(gcf, fullfile([path2repo, filesep,'RErrvsNodes_', mov,'_' tra]), '-dpng', '-r300');
    end
    
    
    figure();
    
    if strcmp(dim, '2D')
        % Set figure properties
        set(gcf, 'Units', 'Inches', 'Position', [1, 1, column_width_in/2, height_in/3]);
        set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [column_width_in/2, height_in/3]);
    else
        set(gcf, 'Units', 'Inches', 'Position', [1, 1, column_width_in/2, height_in/2]);
        set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [column_width_in/2, height_in/2]);
    end
    
    % First subplot: ErrorT
    hold on;
    
    for i = 1:nMethods
        % Extract values for the current method: size [n, m]
        values = squeeze(serrorR(:, i, :)); % size: [nNodes x mSamples] 
        % Compute mean and standard deviation across samples (dim = 2)
        meanVals = nanmean(values, 1); % size: [nNodes x 1]
        stdVals  = nanstd(values, 0, 1); % standard deviation
        % Plot mean with error bars
        errorbar(nNodes, meanVals, stdVals, '-o', ...
            'Color', colors(i,:), 'LineWidth', LineWidth, ...
            'DisplayName', methods{i}, ...
            'CapSize', 4, 'MarkerSize', 6);
    end
    
    hold off;
    xlabel('Number of nodes');
    ylabel('Dynamic error E^{2} (rad/s)');
    %title('Dynamic Error vs. nNodes');
    set(gca, 'YScale', 'log');
    set(gca,'FontName', FontName);
    set(gca,'FontSize', FontSize);
    set(gca,'LineWidth', round(LineWidth/2));
    %set(gca,'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', FontSize*2/3, 'Orientation','vertical');    grid on;
    
    if size(nNodes,1) > 1
        print(gcf, fullfile([path2repo, filesep,'SRErrvsNodes_', mov,'_' tra]), '-depsc', '-r300');
        print(gcf, fullfile([path2repo, filesep,'SRErrvsNodes_', mov,'_' tra]), '-dpng', '-r300');
    end
    
    figure();
    
    if strcmp(dim, '2D')
        % Set figure properties
        set(gcf, 'Units', 'Inches', 'Position', [1, 1, column_width_in/2, height_in/3]);
        set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [column_width_in/2, height_in/3]);
    else
        set(gcf, 'Units', 'Inches', 'Position', [1, 1, column_width_in/2, height_in/2]);
        set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [column_width_in/2, height_in/2]);
    end
    
    
    % First subplot: ErrorT
    %subplot(2,1,1);
    hold on;
    
    for i = 1:nMethods
        % Extract values for the current method: size [n, m]
        values = squeeze(cTime(:, i, :)); % size: [nNodes x mSamples]
    
        % Compute mean and standard deviation across samples (dim = 2)
        meanVals = nanmean(values, 1); % size: [nNodes x 1]
        stdVals  = nanstd(values, 0, 1); % standard deviation
    
        % Plot mean with error bars
        errorbar(nNodes, meanVals, stdVals, '-o', ...
            'Color', colors(i,:), 'LineWidth', LineWidth, ...
            'DisplayName', methods{i}, ...
            'CapSize', 4, 'MarkerSize', 6);
    end
    
    hold off;
    xlabel('Number of nodes');
    ylabel('CPU Time (s)');
    %title('Dynamic Error vs. nNodes');
    set(gca, 'YScale', 'log');
    set(gca,'FontName', FontName);
    set(gca,'FontSize', FontSize);
    set(gca,'LineWidth', round(LineWidth/2));
    %set(gca,'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', FontSize*2/3, 'Orientation','vertical');    grid on;
    
    if size(nNodes,1) > 1
        print(gcf, fullfile([path2repo, filesep,'CpuTime_', mov,'_' tra]), '-depsc', '-r300');
        print(gcf, fullfile([path2repo, filesep,'CpuTime_', mov,'_' tra]), '-dpng', '-r300');
    end
    
    
    
    %CPU vs. Dyn Error
    figure();
    
    % Set figure properties
    set(gcf, 'Units', 'Inches', 'Position', [1, 1, column_width_in, height_in/3*2]);
    set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [column_width_in, height_in/3*2]);
    set(gcf,'WindowStyle','normal')
    
    hold on;
    
    % Preallocate for all points
    allCPU = [];
    allErr = [];
    allMethodID = [];
    allNodeIdx = [];
    
    markers = {'hexagram', 'o', 's', 'd', '^', 'v', 'pentagram'};  % One per node
    
    if strcmp(dim,'2D')
        node_counts = [10, 25, 50, 75, 100, 150, 200];
    else
        node_counts = [10, 25, 50, 75, 100];
    end
    
    textHandles = []; % store handles for later adjustment
    
    for i = 1:nMethods
    
        % Extract values for current method
        values = squeeze(cTime(:, i, :));      % [participants x nodes]
        values_e = squeeze(errorR(:, i, :));   % [participants x nodes]
    
        meanVals = nanmean(values, 1);         % [1 x nNodes]
        meanVals_e = nanmean(values_e, 1);     % [1 x nNodes]
    
        % Store for Pareto computation
        allCPU = [allCPU, meanVals];
        allErr = [allErr, meanVals_e];
        allMethodID = [allMethodID, repmat(i, 1, length(meanVals))];
        allNodeIdx = [allNodeIdx, 1:length(meanVals)];
    
        % Get color from plot order
        colorOrder = get(gca, 'ColorOrder');
        color = colorOrder(mod(i-1, size(colorOrder, 1)) + 1, :);
    
    
        for j = 1:length(meanVals)
            % Save point info for later
            thisCPU = meanVals(j);
            thisErr = meanVals_e(j);
    
            if j == 1
                % Markers with label
                scatter(thisCPU, thisErr, 200, ...%, 100*ones(size(thisCPU)),...
                    'Marker', markers{j}, ...
                    'MarkerFaceColor', color, ...
                    'MarkerEdgeColor', color, ...
                    'DisplayName', methods{i});
            else
                scatter(thisCPU, thisErr, 200,...%, 100*ones(size(thisCPU)),...
                    'Marker', markers{j}, ...
                    'MarkerEdgeColor', color,...
                    'MarkerFaceColor', color);
            end

            % % dummy data
            % scatter(nan, nan, 200,...%, 100*ones(size(thisCPU)),...
            %     'Marker', markers{j}, ...
            %     'MarkerEdgeColor', color,...
            %     'MarkerFaceColor', color);
            % 
            % %Text labelV

            % %Put 2,4, and 6 right
            %SR
           % if i == 2 || i == 4  || (i == 6 && (j < 4)) || (i == 5 && (j > 3 || j == 1)) ||(i == 3 && (j > 3))            % vcut   
               %CR

                %VC
            if i == 2 || (i == 4 && j > 3) || i == 1 || (i == 6) || (i == 3 && (j < 3))
                text(thisCPU*1.06, thisErr*1.05, ...
                    num2str(node_counts(j)), ...
                    'FontSize', FontSize, ...
                    'Color', color, ...
                    'HandleVisibility', 'off');
            else
                if j > 4
                    text(thisCPU*0.72, thisErr*1.05, ...
                        num2str(node_counts(j)), ...
                        'FontSize', FontSize, ...
                        'Color', color, ...
                        'HandleVisibility', 'off');
                else
                    text(thisCPU*0.79, thisErr*1.05, ...
                        num2str(node_counts(j)), ...
                        'FontSize',FontSize, ...
                        'Color', color, ...
                        'HandleVisibility', 'off');
                end
            end
        end
    end
    
    % === Axis and aesthetics ===
    ylabel('Dynamic error E^{1} (deg)');
    xlabel('CPU Time (s)');
    xlim([min(allCPU)*0.75, max(allCPU)*1.45]);
    ylim([min(allErr)*0.7, max(allErr)*1.3]);
    set(gca, 'YScale', 'log');
    set(gca, 'XScale', 'log');
    set(gca,'FontName', FontName);
    set(gca,'FontSize', FontSize);
    set(gca,'LineWidth', round(LineWidth/2));
    if strcmp(dim,'2D')
        legend(['','BE','','','','','','', 'ME', '','','','','','', 'LIIIc-2', '','','','','','','RIIa-2',  '','','','','','','LIIIc-3', '','','','','','','RIIa-3'],...
            'Location', 'southwest', 'FontSize', FontSize);%, 'Orientation','vertical');
    else
        methodHandles = findobj(gca, 'Type','Scatter');
        methodHandles = flipud(methodHandles);
        methodHandles = methodHandles(3:nMethods-1:size(methodHandles,1));
        % === Dummy handles for node counts (symbols) ===
        hSymbols = gobjects(length(node_counts),1);
        for j = 1:length(node_counts)
            hSymbols(j) = scatter(nan, nan, 200, ...
                'Marker', markers{j}, ...
                'MarkerEdgeColor', 'k', ...
                'MarkerFaceColor', 'k', ...
                'DisplayName', num2str(node_counts(j)));
        end

        for h = size(hSymbols,1)
            hSymbols(h).LineWidth = 2;
        end

        % === Combine into one legend ===
        allHandles = [methodHandles(:); hSymbols(:)];
        allLabels  = [methods(:); cellstr(num2str(node_counts(:)))];

        lgd = legend(allHandles, allLabels, ...
                      'NumColumns', 2, ...
                      'Location', 'southwest');           
    end
    grid on;
    
    set(gca, 'LooseInset', max(get(gca, 'TightInset'), 0.02));

   % Draw vertical line between columns
    pos = lgd.Position;   % [x y w h] in normalized figure units
    xSep = pos(1) + pos(3)/5*3;    % halfway across legend
    annotation('line', [xSep xSep], [pos(2) pos(2)+pos(4)], ...
    'Color','k','LineWidth',1.5);
     pos = lgd.Position;   % [x y w h] in normalized figure units
    xSep = pos(1) + pos(3)/5*3;    % halfway across legend
    annotation('line', [xSep xSep], [pos(2) pos(2)+pos(4)], ...
    'Color','k','LineWidth',1.5);
    
    if size(nNodes,1) > 1
        print(gcf, fullfile([path2repo, filesep,'CpuTimeDynErr_', mov,'_' tra]), '-depsc', '-r300');
        print(gcf, fullfile([path2repo, filesep,'CpuTimeDynErr_', mov,'_' tra]), '-dpng', '-r300');
    end
end