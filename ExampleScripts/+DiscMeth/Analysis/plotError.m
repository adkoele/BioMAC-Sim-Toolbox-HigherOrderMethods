%======================================================================
%> @file    +DiscMeth/plotError.m
%> @brief   Plot discretization method errors and dynamic residuals.
%>
%> @details
%> This function visualizes the error distribution and dynamic residuals 
%> of different discretization methods across participants, nodes, and 
%> joints.
%>
%>
%> @param   path2repo       String: Path to the repository 
%> @param   resultRunning   Struct: Struct with simulation results for a 
%>                                   specific discretization method and trial.
%> @param   dynResAll       Table: Table containing error metrics and 
%>                                   residual values across methods.
%> @param   column_width_in Double: Plot width in inches (for figure scaling).
%> @param   height_in       Double: Plot height in inches (for figure scaling).
%> @param   numInterp       Int: Number of interpolation points between 
%>                                   collocation nodes.
%> @param   indJoint        Int: Index of the joint being analyzed.
%> @param   startNode       Int: First node index for analysis window.
%> @param   endNode         Int: Last node index for analysis window.
%> @param   FontName        String: Font family used for figure text.
%> @param   FontSize        Double: Font size for axis labels and legends.
%> @param   LineWidth       Double: Line width for plotted curves.
%> @param   plotAllTime     Int: Selects error plotting mode:
%>                            - 0: All
%>                            - 1: One stage methods
%>                            - 2: Two stage methods
%>                            - 3: Three stage methods
%> @param   stdSumVals      Double array Standard deviation values of angles 
%>                                   across participants.
%> @param   stdSumValsdt    Double array: Standard deviation values of derivatives 
%>                                   across participants.
%> @param   indJointSim     Int: Index of the simulated joint in the data table.
%> @param   movement        String Name of the movement being analyzed.
%>
%>
%>
%> @author  Alexander Weiss
%> @date    August, 2025
%======================================================================


function plotError(path2repo, resultRunning, dynResAll, column_width_in, height_in, numInterp, indJoint, startNode, endNode, FontName, FontSize, LineWidth, plotAllTime, stdSumVals, stdSumValsdt, indJointSim, movement)
    
    colours = lines(6);
    if plotAllTime==0
        allFlag = '_Short';
    elseif plotAllTime==1
        allFlag = '_OneStage';
        colours = colours(1:2,:);
        color = lines();
    elseif plotAllTime==2
        allFlag = '_TwoStage';
        colours = colours(3:4,:);
    elseif plotAllTime==3
        allFlag = '_ThreeStage';
        colours = colours(5:6,:);
    end
    %resultPlotFolder = [path2repo, filesep ];
    X = resultRunning.X;
    nNodes = resultRunning.problem.nNodesDur;
    h = X(resultRunning.problem.idx.dur)/(nNodes-1);


    
    newSize = numInterp*(nNodes-1);
    newSizeLong = numInterp*(nNodes);
    timeVec = 0:h:X(resultRunning.problem.idx.dur);
    timeVec = timeVec(2:end);
    
    Xq = linspace(0,nNodes-1,newSize);
    XqLong = linspace(0,nNodes-1,newSizeLong);
    XqMin = linspace(0,nNodes-1,newSize-(1*numInterp));
    XNode = 1:nNodes-1;

    thresOneStd = interp1((1:size(stdSumVals,2)),stdSumVals(indJointSim,:),Xq);
    thresTwoStd = 2*interp1((1:size(stdSumVals,2)),stdSumVals(indJointSim,:),Xq);

    thresOneStd(1:100) = thresOneStd(101);
    threTwoStd(1:100) = thresTwoStd(101);

    thresOneStd(size(Xq,2)+1+1:size(Xq,2)+100) = thresOneStd(size(Xq,2));
    thresTwoStd(size(Xq,2)+1+1:size(Xq,2)+100) = thresTwoStd(size(Xq,2));

    thresOneStddt = interp1((1:size(stdSumValsdt,2)),stdSumValsdt(indJointSim,:),Xq);
    thresTwoStddt = 2*interp1((1:size(stdSumValsdt,2)),stdSumValsdt(indJointSim,:),Xq);
    % Plot single joint errors
    figure();
    
    % Set figure properties
    set(gcf, 'Units', 'Inches', 'Position', [1, 1, column_width_in/2, height_in/2]);
    set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [column_width_in/2, height_in/2]);
    
    for iPlot = 1:size(dynResAll,1)
        hold on;
        toPlot = dynResAll.DynRes{iPlot,1};
        toPlotSum = ((toPlot(:,indJoint)));
        plot(toPlotSum, 'LineWidth',LineWidth/2, 'Color', colours(iPlot,:));
        if contains(resultRunning.problem.model.states.name{indJoint}, 'pelvis_t') & ~strcmp(resultRunning.problem.model.states.name{indJoint}, 'pelvis_tilt')
            disp([dynResAll.Method{iPlot,1}, '_trans: ' num2str(dynResAll.ErrorT(iPlot))]);
            ylabel('Dynamic error e_{q}^{1}) (m/s)');
        else
            disp([dynResAll.Method{iPlot,1}, '_rad: ' num2str(dynResAll.ErrorR(iPlot))]);
            ylabel('Dynamic error e_{q}^{1} (rad/s)');
        end      
    end
    plot(XNode*numInterp, zeros(size(XNode,2)) , 'ro');
   

    plot(thresOneStd, 'LineWidth', 2, 'LineStyle', '--', 'Color', 'r');
    plot(thresTwoStd, 'LineWidth', 2, 'Color', 'r');
    plot(-thresOneStd, 'LineWidth', 2, 'LineStyle', '--', 'Color', 'r');
    plot(-thresTwoStd, 'LineWidth', 2, 'Color', 'r');
    
    if plotAllTime==0
    else
        ylim([-8, 8]);
        %xlim([size(Xq,2)-(nNodes-startNode-1)*numInterp, size(Xq,2)-(nNodes-endNode-2)*numInterp]);%size(Xq,2)]);
    end

    set(gca,'XTick',numInterp:numInterp:size(Xq,2));
    set(gca,'XTickLabel',1:1:dynResAll.nNodes(iPlot)-1);
    set(gca,'FontName', FontName);
    set(gca,'FontSize', FontSize);
    set(gca,'LineWidth', LineWidth);
   % set(gca,'FontWeight', 'bold');
    legMeth = reshape([dynResAll.Method], size(dynResAll.Method));
    legend('1SD', '2SD', '', '', legMeth, 'Location', 'southeast', 'FontSize', FontSize);%'Orientation', 'horizontal'
    xlabel('Nodes');
    hold on;
    
    if contains(resultRunning.problem.model.states.name{indJoint}, 'pelvis_t') & ~strcmp(resultRunning.problem.model.states.name{indJoint}, 'pelvis_tilt')
        print(gcf, fullfile([path2repo,filesep,movement,'_JointTransErr', allFlag]), '-depsc', '-r300');
        print(gcf, fullfile([path2repo,filesep,movement,'_JointTransErr', allFlag]), '-dpng', '-r300');
    else
        print(gcf, fullfile([path2repo,filesep,movement,'_JointDynErr', allFlag]), '-depsc', '-r300');
        print(gcf, fullfile([path2repo,filesep,movement,'_JointDynErr', allFlag]), '-dpng', '-r300');
    end


    %% Plot single joint errors normalized
    figure();
    
    % Set figure properties
    set(gcf, 'Units', 'Inches', 'Position', [1, 1, column_width_in/2, height_in/2]);
    set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [column_width_in/2, height_in/2]);
    
    for iPlot = 1:size(dynResAll,1)
        hold on;
        toPlot = dynResAll.DynRes{iPlot,1};
        toPlotSum = ((toPlot(:,indJoint)));
        plot(toPlotSum./thresOneStd', 'LineWidth',LineWidth/2 , 'Color', colours(iPlot,:));
        if contains(resultRunning.problem.model.states.name{indJoint}, 'pelvis_t') & ~strcmp(resultRunning.problem.model.states.name{indJoint}, 'pelvis_tilt')
            disp([dynResAll.Method{iPlot,1}, '_trans: ' num2str(dynResAll.ErrorT(iPlot))]);
            ylabel('Dynamic error e_{q}^{1} (SD)');
        else
            disp([dynResAll.Method{iPlot,1}, '_rad: ' num2str(dynResAll.ErrorR(iPlot))]);
            ylabel('Dynamic error e_{q}^{1} (SD)');
        end      
    end
    yline(1, 'LineWidth', 2, 'LineStyle', '--', 'Color', 'r');
    yline(2, 'LineWidth', 2, 'Color', 'r');
    yline(-1, 'LineWidth', 2, 'LineStyle', '--', 'Color', 'r');
    yline(-2, 'LineWidth', 2, 'Color', 'r');
    if plotAllTime == 0
    else
        ylim([-8 8]);
    end
    plot(XNode(1:end-1)*numInterp, zeros(size(XNode,2)-1) , 'ro');

    %ylim([-8, 8]);
    set(gca,'XTick',numInterp:numInterp:size(Xq,2));
    set(gca,'XTickLabel',1:1:dynResAll.nNodes(iPlot)-1);
    set(gca,'FontName', FontName);
    set(gca,'FontSize', FontSize);
    set(gca,'LineWidth', LineWidth);
   % set(gca,'FontWeight', 'bold');
    legMeth = reshape([dynResAll.Method], size(dynResAll.Method));
    legend('1SD', '2SD', '', '', legMeth, 'Location', 'southeast', 'FontSize', FontSize);%'Orientation', 'horizontal'
    xlabel('Nodes');


    if contains(resultRunning.problem.model.states.name{indJoint}, 'pelvis_t') & ~strcmp(resultRunning.problem.model.states.name{indJoint}, 'pelvis_tilt')
        print(gcf, fullfile([path2repo,filesep,movement,'_JointTransErrNorm', allFlag]), '-depsc', '-r300');
        print(gcf, fullfile([path2repo,filesep,movement,'_JointTransErrNorm', allFlag]), '-dpng', '-r300');
    else
        print(gcf, fullfile([path2repo,filesep,movement,'_JointDynErrNorm', allFlag]), '-depsc', '-r300');
        print(gcf, fullfile([path2repo,filesep,movement,'_JointDynErrNorm', allFlag]), '-dpng', '-r300');
    end

    
    %% Second order single joint error
    figure();
    % Set figure properties
    set(gcf, 'Units', 'Inches', 'Position', [1, 1, column_width_in/2, height_in/2]);
    set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [column_width_in/2, height_in/2]);
    %subplot(2,1,2);
    %subtitle('First order error angles' , 'FontSize', FontSize, 'FontWeight', 'bold', 'FontName', FontName);
    for iPlot = 1:size(dynResAll,1)
        hold on;
        toPlot = dynResAll.DynResDot{iPlot,1};
        toPlotSum = ((toPlot(:,indJoint)));
        plot(toPlotSum, 'LineWidth',LineWidth/2, 'Color', colours(iPlot,:));
        if contains(resultRunning.problem.model.states.name{indJoint}, 'pelvis_t') & ~strcmp(resultRunning.problem.model.states.name{indJoint}, 'pelvis_tilt')
            disp([dynResAll.Method{iPlot,1}, '_transSec: ' num2str(dynResAll.SErrorT(iPlot))]);
            ylabel('Dynamic error e^{2} (m/s²)');
        else
            disp([dynResAll.Method{iPlot,1}, '_radSec: ' num2str(dynResAll.SErrorR(iPlot))]);
            ylabel('Dynamic error e^{2} (rad/s²)');
        end
        %disp([dynResAll.Method{iPlot,1}, '_radSEC: ' num2str(dynResAll.SErrorR(iPlot))]);
        %ylabel('Second order Dynamic error (rad/s^2)');
    end

    % plot(thresOneStddt, 'LineWidth', 2, 'LineStyle', '--', 'Color', 'r');
    % plot(thresTwoStddt, 'LineWidth', 2, 'Color', 'r');
    % plot(-thresOneStddt, 'LineWidth', 2, 'LineStyle', '--', 'Color', 'r');
    % plot(-thresTwoStddt, 'LineWidth', 2, 'Color', 'r');

    plot(XNode*numInterp, zeros(size(XNode,2)) , 'ro');
    if plotAllTime==0
       % xlim([size(Xq,2)-(nNodes-startNode-1)*numInterp, size(Xq,2)-(nNodes-endNode-2)*numInterp]);%size(Xq,2)]);
    end
    ylim([-8, 8]);
    set(gca,'XTick',numInterp:numInterp:size(Xq,2));
    set(gca,'XTickLabel',1:1:dynResAll.nNodes(iPlot)-1);
    set(gca,'FontName', FontName);
    set(gca,'FontSize', FontSize);
    set(gca,'LineWidth', LineWidth);
    %set(gca,'FontWeight', 'bold');
    legMeth = reshape([dynResAll.Method], size(dynResAll.Method));
    legend(legMeth, 'Location', 'southeast', 'FontSize', FontSize);%'Orientation', 'horizontal'
    xlabel('Nodes');
    hold off;
    
    if contains(resultRunning.problem.model.states.name{indJoint}, 'pelvis_t') & ~strcmp(resultRunning.problem.model.states.name{indJoint}, 'pelvis_tilt')
        print(gcf, fullfile([path2repo,filesep,movement,'_JointTransSecErr', allFlag]), '-depsc', '-r300');
        print(gcf, fullfile([path2repo,filesep,movement,'_JointTransSecErr', allFlag]), '-dpng', '-r300');
    else
        print(gcf, fullfile([path2repo,filesep,movement,'_JointDynSecErr', allFlag]), '-depsc', '-r300');
        print(gcf, fullfile([path2repo,filesep,movement,'_JointDynSecErr', allFlag]), '-dpng', '-r300');
    end

    % Plot summed error angles
    figure();
    
    % Set figure properties
    set(gcf, 'Units', 'Inches', 'Position', [1, 1, column_width_in/2, height_in/2]);
    set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [column_width_in/2, height_in/2]);
    
    for iPlot = 1:size(dynResAll,1)
        hold on;
        toPlot = dynResAll.SumErrorR{iPlot,1};
        plot(toPlot, 'LineWidth',LineWidth/2, 'Color', colours(iPlot,:));
        disp([dynResAll.Method{iPlot,1}, '_rad: ' num2str(dynResAll.IntErrorR(iPlot))]);
        ylabel('Summed dynamic error e^{1} (rad/s)');
    end
    plot(XNode*numInterp, zeros(size(XNode,2)) , 'ro');
    if plotAllTime==0
       % xlim([size(Xq,2)-(nNodes-startNode-1)*numInterp, size(Xq,2)-(nNodes-endNode-2)*numInterp]);%size(Xq,2)]);
    end
      
    %yline(2 * mean(stdSumVals));
    set(gca,'XTick',numInterp:numInterp:size(Xq,2));
    set(gca,'XTickLabel',1:1:dynResAll.nNodes(iPlot)-1);
    set(gca,'FontName', FontName);
    set(gca,'FontSize', FontSize);
    set(gca,'LineWidth', round(LineWidth));
    %set(gca,'FontWeight', 'bold');
    legMeth = reshape([dynResAll.Method], size(dynResAll.Method));
    legend(legMeth, 'Location', 'northeast', 'FontSize', FontSize);%'Orientation', 'horizontal'
    xlabel('Nodes');
    hold off;
    
    print(gcf, fullfile([path2repo,filesep,movement,'_SumDynErr', allFlag]), '-depsc', '-r300');
    print(gcf, fullfile([path2repo,filesep,movement,'_SumDynErr', allFlag]), '-dpng', '-r300');
    
    figure();
    % Set figure properties
    set(gcf, 'Units', 'Inches', 'Position', [1, 1, column_width_in/2, height_in/2]);
    set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [column_width_in/2, height_in/2]);
    %subplot(2,1,2);
    %subtitle('First order error angles' , 'FontSize', FontSize, 'FontWeight', 'bold', 'FontName', FontName);
    for iPlot = 1:size(dynResAll,1)
        hold on;
        toPlot = dynResAll.SumSecErrorR{iPlot,1};
        plot(toPlot, 'LineWidth',LineWidth/2, 'Color', colours(iPlot,:));
        disp([dynResAll.Method{iPlot,1}, '_radSEC: ' num2str(dynResAll.IntSecErrorR(iPlot))]);
        ylabel('Summed dynamic error e^{2} (rad/s²)');
    end
    
     plot(XNode*numInterp, zeros(size(XNode,2)) , 'ro');
    
    if plotAllTime==0
      %  xlim([size(Xq,2)-(nNodes-startNode-1)*numInterp, size(Xq,2)-(nNodes-endNode-2)*numInterp]);%size(Xq,2)]);
    end
    
    set(gca,'XTick',numInterp:numInterp:size(Xq,2));
    set(gca,'XTickLabel',1:1:dynResAll.nNodes(iPlot)-1);
    set(gca,'FontName', FontName);
    set(gca,'FontSize', FontSize);
    set(gca,'LineWidth', LineWidth);
    %set(gca,'FontWeight', 'bold');
    legMeth = reshape([dynResAll.Method], size(dynResAll.Method));
    legend(legMeth, 'Location', 'northeast', 'FontSize', FontSize);%'Orientation', 'horizontal'
    xlabel('Nodes');
    hold off;
    
    print(gcf, fullfile([path2repo,filesep,movement,'_SumDynSecErr', allFlag]), '-depsc', '-r300');
    print(gcf, fullfile([path2repo,filesep,movement,'_SumDynSecErr', allFlag]), '-dpng', '-r300');


    % Plot summed error translations
    figure();
    
    % Set figure properties
    set(gcf, 'Units', 'Inches', 'Position', [1, 1, column_width_in/2, height_in/2]);
    set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [column_width_in/2, height_in/2]);
    
    for iPlot = 1:size(dynResAll,1)
        hold on;
        toPlot = dynResAll.SumErrorT{iPlot,1};
        plot(toPlot, 'LineWidth',LineWidth/2, 'Color', colours(iPlot,:));
        disp([dynResAll.Method{iPlot,1}, '_trans: ' num2str(dynResAll.IntErrorT(iPlot))]);
        ylabel('Summed dynamic error e^{1} (m/s)');
    end
    plot(XNode*numInterp, zeros(size(XNode,2)) , 'ro');
    if plotAllTime==0
       % xlim([size(Xq,2)-(nNodes-startNode-1)*numInterp, size(Xq,2)-(nNodes-endNode-2)*numInterp]);%size(Xq,2)]);
    end
      
    set(gca,'XTick',numInterp:numInterp:size(Xq,2));
    set(gca,'XTickLabel',1:1:dynResAll.nNodes(iPlot)-1);
    set(gca,'FontName', FontName);
    set(gca,'FontSize', FontSize);
    set(gca,'LineWidth', LineWidth);
    %set(gca,'FontWeight', 'bold');
    legMeth = reshape([dynResAll.Method], size(dynResAll.Method));
    legend(legMeth, 'Location', 'northeast', 'FontSize', FontSize);%'Orientation', 'horizontal'
    xlabel('Nodes');
    hold off;
    
    print(gcf, fullfile([path2repo,filesep,movement,'_SumTransErr', allFlag]), '-depsc', '-r300');
    print(gcf, fullfile([path2repo,filesep,movement,'_SumTransErr', allFlag]), '-dpng', '-r300');
    
    figure();
    % Set figure properties
    set(gcf, 'Units', 'Inches', 'Position', [1, 1, column_width_in/2, height_in/2]);
    set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [column_width_in/2, height_in/2]);
    %subplot(2,1,2);
    %subtitle('First order error angles' , 'FontSize', FontSize, 'FontWeight', 'bold', 'FontName', FontName);
    for iPlot = 1:size(dynResAll,1)
        hold on;
        toPlot = dynResAll.SumSecErrorT{iPlot,1};
        plot(toPlot, 'LineWidth',LineWidth/2, 'Color', colours(iPlot,:));
        disp([dynResAll.Method{iPlot,1}, '_transSEC: ' num2str(dynResAll.IntSecErrorT(iPlot))]);
        ylabel('Summed dynamic error e^{2} (m/s²)');
    end
    plot(XNode*numInterp, zeros(size(XNode,2)) , 'ro');
    if plotAllTime==0
       % xlim([size(Xq,2)-(nNodes-startNode-1)*numInterp, size(Xq,2)-(nNodes-endNode-2)*numInterp]);%size(Xq,2)]);
    end
       
    set(gca,'XTick',numInterp:numInterp:size(Xq,2));
    set(gca,'XTickLabel',1:1:dynResAll.nNodes(iPlot)-1);
    set(gca,'FontName', FontName);
    set(gca,'FontSize', FontSize);
    set(gca,'LineWidth', LineWidth);
    %set(gca,'FontWeight', 'bold');
    legMeth = reshape([dynResAll.Method], size(dynResAll.Method));
    legend(legMeth, 'Location', 'northeast',  'FontSize', FontSize);%'Orientation', 'horizontal'
    xlabel('Nodes');
    hold off;
    
   print(gcf, fullfile([path2repo,filesep,movement,'_SumTransSecErr', allFlag]), '-depsc', '-r300');
    print(gcf, fullfile([path2repo,filesep,movement,'_SumTransSecErr', allFlag]), '-dpng', '-r300');


%     %% Muscles
% 
%      figure();
%     % Set figure properties
%     set(gcf, 'Units', 'Inches', 'Position', [1, 1, column_width_in/2, height_in/2]);
%     set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [column_width_in/2, height_in/2]);
%     %subplot(2,1,2);
%     %subtitle('First order error angles' , 'FontSize', FontSize, 'FontWeight', 'bold', 'FontName', FontName);
%     for iPlot = 1:size(dynResAll,1)
%         hold on;
%         toPlot = dynResAll.SumErrorA{iPlot,1};
%         plot(toPlot, 'LineWidth',LineWidth);
%         disp([dynResAll.Method{iPlot,1}, ': ' num2str(dynResAll.IntErrorA(iPlot))]);
%         ylabel('Summed activation error');
%     end
% 
%      plot(XNode*numInterp, zeros(size(XNode,2)) , 'ro');
% 
%     if plotAllTime==0
%         xlim([size(Xq,2)-(nNodes-startNode-1)*numInterp, size(Xq,2)-(nNodes-endNode-2)*numInterp]);%size(Xq,2)]);
%     end
% 
%     set(gca,'XTick',numInterp:numInterp:size(Xq,2));
%     set(gca,'XTickLabel',1:1:dynResAll.nNodes(iPlot)-1);
%     set(gca,'FontName', FontName);
%     set(gca,'FontSize', FontSize);
%     set(gca,'LineWidth', LineWidth);
%     %set(gca,'FontWeight', 'bold');
%     legMeth = reshape([dynResAll.Method], size(dynResAll.Method));
%     legend(legMeth, 'Location', 'northeast', 'FontSize', FontSize);%'Orientation', 'horizontal'
%     xlabel('Nodes');
%     hold off;
% 
%          figure();
%     % Set figure properties
%     set(gcf, 'Units', 'Inches', 'Position', [1, 1, column_width_in/2, height_in/2]);
%     set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [column_width_in/2, height_in/2]);
%     %subplot(2,1,2);
%     %subtitle('First order error angles' , 'FontSize', FontSize, 'FontWeight', 'bold', 'FontName', FontName);
%     for iPlot = 1:size(dynResAll,1)
%         hold on;
%         toPlot = dynResAll.SumErrorS{iPlot,1};
%         plot(toPlot, 'LineWidth',LineWidth);
%         disp([dynResAll.Method{iPlot,1}, ': ' num2str(dynResAll.IntErrorS(iPlot))]);
%         ylabel('Summed length error');
%     end
% 
%      plot(XNode*numInterp, zeros(size(XNode,2)) , 'ro');
% 
%     if plotAllTime==0
%         xlim([size(Xq,2)-(nNodes-startNode-1)*numInterp, size(Xq,2)-(nNodes-endNode-1)*numInterp]);%size(Xq,2)]);
%     end
% 
%     set(gca,'XTick',numInterp:numInterp:size(Xq,2));
%     set(gca,'XTickLabel',1:1:dynResAll.nNodes(iPlot)-1);
%     set(gca,'FontName', FontName);
%     set(gca,'FontSize', FontSize);
%     set(gca,'LineWidth', LineWidth);
%     %set(gca,'FontWeight', 'bold');
%     legMeth = reshape([dynResAll.Method], size(dynResAll.Method));
%     legend(legMeth, 'Location', 'northeast', 'FontSize', FontSize);%'Orientation', 'horizontal'
%     xlabel('Nodes');
%     hold off;
% 
end