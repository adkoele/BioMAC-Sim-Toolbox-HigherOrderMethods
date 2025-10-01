%======================================================================
%> @file    +DiscMeth/analyzeANOVA.m
%> @brief   Perform ANOVA analysis on discretization method results.
%>
%> @details
%> This function applies statistical analysis (ANOVA) to evaluate 
%> differences in dynamic residuals and error metrics between 
%> discretization methods. 
%>
%> @param   path2repo   String: Path to the repository.
%> @param   dynResTable Table: Table with dynamic residuals, errors, 
%>                              and timing information from simulations.
%> @param   movement    String: Name of the movement (e.g. 'straightslowrunning').
%> @param   dim         String: Model dimensionality: '2D' or '3D'.
%> @param   FontName    String: Font family used in plots.
%> @param   FontSize    Double: Font size for plot labels and legends.
%>
%>
%> @author  Alexander Weiss
%> @date    August, 2025
%======================================================================


function analyzeANOVA(path2repo, dynResTable, movement, dim, FontName, FontSize)
    % Convert categorical variables if needed
    dynResTable.Method = categorical(dynResTable.Method);
    dynResTable.Tracking = categorical(dynResTable.Tracking);
    dynResTable.Movement = categorical(dynResTable.Movement);
    dynResTable.nNodes = categorical(dynResTable.nNodes);
    dynResTable.Participant = categorical(dynResTable.Participant);% Create wide table (reshape if needed) so that each method is a separate column
   
    %% Trans Error
    disp('Trans FO Error');
    wideTable = unstack(dynResTable(:, {'Participant','nNodes', 'Method', 'IntErrorT'}), ...
        'IntErrorT', 'Method');
    % Define within-subject factor (Method)
    withinDesign = table(categorical({'BE','ME','RIIa_2','LIIIc_2','RIIa_3','LIIIc_3'})', ...
                         'VariableNames', {'Method'});% Fit the repeated measures model
    rm = fitrm(wideTable, ...
        'BE,ME,RIIa_2,RIIa_3,LIIIc_2,LIIIc_3 ~ nNodes', ...
        'WithinDesign', withinDesign);

    ranova_tbl_errorT = ranova(rm, 'WithinModel', 'Method');% Show results
   % disp(ranova_tbl_errorT);
    res_error_errorT = multcompare(rm, 'Method', 'ComparisonType', 'bonferroni');

    % Extract raw p-values
    p_raw = ranova_tbl_errorT.pValue;
    terms = ranova_tbl_errorT.Properties.RowNames;% Bonferroni correction
    p_bonf = min(p_raw * numel(p_raw), 1);  % Cap at 1% Holm-Bonferroni (step-down)
    [~, sortIdx] = sort(p_raw);
    rank = 1:length(p_raw);
    p_holm = zeros(size(p_raw));
    for i = 1:length(p_raw)
        p_holm(sortIdx(i)) = min((numel(p_raw)-rank(i)+1) * p_raw(sortIdx(i)), 1);
    end% Add to table
    ranova_tbl_errorT.pBonferroni = p_bonf;
    ranova_tbl_errorT.pHolm = p_holm;% Display
    disp(ranova_tbl_errorT(:, {'SumSq', 'DF', 'F', 'pValue', 'pBonferroni', 'pHolm'}));

    lme_errorT = fitlme(dynResTable, ...
        'IntErrorT ~ Method * nNodes + (1|Participant)');
    %disp(anova(lme_errorT));



      %% Rad Error
      disp('Joint Angle FO Error');
    wideTable = unstack(dynResTable(:, {'Participant','nNodes', 'Method', 'IntErrorR'}), ...
        'IntErrorR', 'Method');
    % Define within-subject factor (Method)
    withinDesign = table(categorical({'BE','ME','RIIa_2','LIIIc_2','RIIa_3','LIIIc_3'})', ...
                         'VariableNames', {'Method'});% Fit the repeated measures model
    rm = fitrm(wideTable, ...
        'BE,ME,RIIa_2,RIIa_3,LIIIc_2,LIIIc_3 ~ nNodes', ...
        'WithinDesign', withinDesign);

    ranova_tbl_errorR = ranova(rm, 'WithinModel', 'Method');% Show results
    %disp(ranova_tbl_errorR);
    res_errorR = multcompare(rm, 'Method', 'ComparisonType', 'bonferroni');

    lme_errorR = fitlme(dynResTable, ...
        'IntErrorR ~ Method * nNodes + (1|Participant)');
   % disp(anova(lme_errorR));

       % Extract raw p-values
    p_raw = ranova_tbl_errorR.pValue;
    terms = ranova_tbl_errorR.Properties.RowNames;% Bonferroni correction
    p_bonf = min(p_raw * numel(p_raw), 1);  % Cap at 1% Holm-Bonferroni (step-down)
    [~, sortIdx] = sort(p_raw);
    rank = 1:length(p_raw);
    p_holm = zeros(size(p_raw));
    for i = 1:length(p_raw)
        p_holm(sortIdx(i)) = min((numel(p_raw)-rank(i)+1) * p_raw(sortIdx(i)), 1);
    end% Add to table
    ranova_tbl_errorR.pBonferroni = p_bonf;
    ranova_tbl_errorR.pHolm = p_holm;% Display
    disp(ranova_tbl_errorR(:, {'SumSq', 'DF', 'F', 'pValue', 'pBonferroni', 'pHolm'}));


    %% Second Order
    %% Trans Error
     disp('Translation SO Error');
    wideTable = unstack(dynResTable(:, {'Participant','nNodes', 'Method', 'IntSecErrorT'}), ...
        'IntSecErrorT', 'Method');
    % Define within-subject factor (Method)
    withinDesign = table(categorical({'BE','ME','RIIa_2','LIIIc_2','RIIa_3','LIIIc_3'})', ...
                         'VariableNames', {'Method'});% Fit the repeated measures model
    rm = fitrm(wideTable, ...
        'BE,ME,RIIa_2,RIIa_3,LIIIc_2,LIIIc_3 ~ nNodes', ...
        'WithinDesign', withinDesign);

    ranova_tbl_SerrorT = ranova(rm, 'WithinModel', 'Method');% Show results
   % disp(ranova_tbl_errorT);
    res_error_SerrorT = multcompare(rm, 'Method', 'ComparisonType', 'bonferroni');

    % Extract raw p-values
    p_raw = ranova_tbl_SerrorT.pValue;
    terms = ranova_tbl_SerrorT.Properties.RowNames;% Bonferroni correction
    p_bonf = min(p_raw * numel(p_raw), 1);  % Cap at 1% Holm-Bonferroni (step-down)
    [~, sortIdx] = sort(p_raw);
    rank = 1:length(p_raw);
    p_holm = zeros(size(p_raw));
    for i = 1:length(p_raw)
        p_holm(sortIdx(i)) = min((numel(p_raw)-rank(i)+1) * p_raw(sortIdx(i)), 1);
    end% Add to table
    ranova_tbl_SerrorT.pBonferroni = p_bonf;
    ranova_tbl_SerrorT.pHolm = p_holm;% Display
    disp(ranova_tbl_SerrorT(:, {'SumSq', 'DF', 'F', 'pValue', 'pBonferroni', 'pHolm'}));

    lme_SerrorT = fitlme(dynResTable, ...
        'IntSecErrorT ~ Method * nNodes + (1|Participant)');
    %disp(anova(lme_errorT));

      %% Rad Error
   disp('Joint angle SO Error');
    wideTable = unstack(dynResTable(:, {'Participant','nNodes', 'Method', 'IntSecErrorR'}), ...
        'IntSecErrorR', 'Method');
    % Define within-subject factor (Method)
    withinDesign = table(categorical({'BE','ME','RIIa_2','LIIIc_2','RIIa_3','LIIIc_3'})', ...
                         'VariableNames', {'Method'});% Fit the repeated measures model
    rm = fitrm(wideTable, ...
        'BE,ME,RIIa_2,RIIa_3,LIIIc_2,LIIIc_3 ~ nNodes', ...
        'WithinDesign', withinDesign);

    ranova_tbl_SerrorR = ranova(rm, 'WithinModel', 'Method');% Show results
    %disp(ranova_tbl_errorR);
    res_SerrorR = multcompare(rm, 'Method', 'ComparisonType', 'bonferroni');

    lme_SerrorR = fitlme(dynResTable, ...
        'IntSecErrorR ~ Method * nNodes + (1|Participant)');
   % disp(anova(lme_errorR));

       % Extract raw p-values
    p_raw = ranova_tbl_SerrorR.pValue;
    terms = ranova_tbl_SerrorR.Properties.RowNames;% Bonferroni correction
    p_bonf = min(p_raw * numel(p_raw), 1);  % Cap at 1% Holm-Bonferroni (step-down)
    [~, sortIdx] = sort(p_raw);
    rank = 1:length(p_raw);
    p_holm = zeros(size(p_raw));
    for i = 1:length(p_raw)
        p_holm(sortIdx(i)) = min((numel(p_raw)-rank(i)+1) * p_raw(sortIdx(i)), 1);
    end% Add to table
    ranova_tbl_SerrorR.pBonferroni = p_bonf;
    ranova_tbl_SerrorR.pHolm = p_holm;% Display
    disp(ranova_tbl_SerrorR(:, {'SumSq', 'DF', 'F', 'pValue', 'pBonferroni', 'pHolm'}));


    %% Activation Error
     disp('Activation Error');
    wideTable = unstack(dynResTable(:, {'Participant','nNodes', 'Method', 'IntErrorA'}), ...
        'IntErrorA', 'Method');
    % Define within-subject factor (Method)
    withinDesign = table(categorical({'BE','ME','RIIa_2','LIIIc_2','RIIa_3','LIIIc_3'})', ...
                         'VariableNames', {'Method'});% Fit the repeated measures model
    rm = fitrm(wideTable, ...
        'BE,ME,RIIa_2,RIIa_3,LIIIc_2,LIIIc_3 ~ nNodes', ...
        'WithinDesign', withinDesign);

    ranova_tbl_SerrorR = ranova(rm, 'WithinModel', 'Method');% Show results
    %disp(ranova_tbl_errorR);
    res_SerrorR = multcompare(rm, 'Method', 'ComparisonType', 'bonferroni');

    lme_SerrorR = fitlme(dynResTable, ...
        'IntErrorA ~ Method * nNodes + (1|Participant)');
   % disp(anova(lme_errorR));

       % Extract raw p-values
    p_raw = ranova_tbl_SerrorR.pValue;
    terms = ranova_tbl_SerrorR.Properties.RowNames;% Bonferroni correction
    p_bonf = min(p_raw * numel(p_raw), 1);  % Cap at 1% Holm-Bonferroni (step-down)
    [~, sortIdx] = sort(p_raw);
    rank = 1:length(p_raw);
    p_holm = zeros(size(p_raw));
    for i = 1:length(p_raw)
        p_holm(sortIdx(i)) = min((numel(p_raw)-rank(i)+1) * p_raw(sortIdx(i)), 1);
    end% Add to table
    ranova_tbl_SerrorR.pBonferroni = p_bonf;
    ranova_tbl_SerrorR.pHolm = p_holm;% Display
    disp(ranova_tbl_SerrorR(:, {'SumSq', 'DF', 'F', 'pValue', 'pBonferroni', 'pHolm'}));


     %% CE Length Error
    disp('CE Length Error');
    wideTable = unstack(dynResTable(:, {'Participant','nNodes', 'Method', 'IntErrorS'}), ...
        'IntErrorS', 'Method');
    % Define within-subject factor (Method)
    withinDesign = table(categorical({'BE','ME','RIIa_2','LIIIc_2','RIIa_3','LIIIc_3'})', ...
                         'VariableNames', {'Method'});% Fit the repeated measures model
    rm = fitrm(wideTable, ...
        'BE,ME,RIIa_2,RIIa_3,LIIIc_2,LIIIc_3 ~ nNodes', ...
        'WithinDesign', withinDesign);

    ranova_tbl_SerrorR = ranova(rm, 'WithinModel', 'Method');% Show results
    %disp(ranova_tbl_errorR);
    res_SerrorR = multcompare(rm, 'Method', 'ComparisonType', 'bonferroni');

    lme_SerrorR = fitlme(dynResTable, ...
        'IntErrorS ~ Method * nNodes + (1|Participant)');
   % disp(anova(lme_errorR));

       % Extract raw p-values
    p_raw = ranova_tbl_SerrorR.pValue;
    terms = ranova_tbl_SerrorR.Properties.RowNames;% Bonferroni correction
    p_bonf = min(p_raw * numel(p_raw), 1);  % Cap at 1% Holm-Bonferroni (step-down)
    [~, sortIdx] = sort(p_raw);
    rank = 1:length(p_raw);
    p_holm = zeros(size(p_raw));
    for i = 1:length(p_raw)
        p_holm(sortIdx(i)) = min((numel(p_raw)-rank(i)+1) * p_raw(sortIdx(i)), 1);
    end% Add to table
    ranova_tbl_SerrorR.pBonferroni = p_bonf;
    ranova_tbl_SerrorR.pHolm = p_holm;% Display
    disp(ranova_tbl_SerrorR(:, {'SumSq', 'DF', 'F', 'pValue', 'pBonferroni', 'pHolm'}));

    %% CPU Time
  disp('CPU');
    wideTable = unstack(dynResTable(:, {'Participant','nNodes', 'Method', 'CpuTime'}), ...
        'CpuTime', 'Method');
    % Define within-subject factor (Method)
    withinDesign = table(categorical({'BE','ME','RIIa_2','LIIIc_2','RIIa_3','LIIIc_3'})', ...
                         'VariableNames', {'Method'});% Fit the repeated measures model
    rm = fitrm(wideTable, ...
        'BE,ME,RIIa_2,RIIa_3,LIIIc_2,LIIIc_3 ~ nNodes', ...
        'WithinDesign', withinDesign);

    ranova_tbl_cpu = ranova(rm, 'WithinModel', 'Method');% Show results
   % disp(ranova_tbl_cpu);

    res_error_CPU = multcompare(rm, 'Method', 'ComparisonType', 'bonferroni');

    lme_cpu = fitlme(dynResTable, ...
        'CpuTime ~ Method * nNodes + (1|Participant)');
end