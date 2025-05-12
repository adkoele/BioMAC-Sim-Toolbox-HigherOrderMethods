function plotStabRegion(discMeth)

    FS = 15;
    [A,b,c, name] = DiscMeth.getButcherTableau(discMeth);
    b = b';
    
    s = length(b);
    I = eye(s);
    one = ones(s, 1);

    figure();
   
    if strcmp(discMeth, 'RIIa-3')
        [x, y] = meshgrid(linspace(-10,15,600), linspace(-10,10,600));
    else
        [x, y] = meshgrid(linspace(-6,6,600), linspace(-6,6,600));
    end
    z = x + 1i*y;
    R = zeros(size(z));

    % Evaluate R(z) over grid
    for i = 1:numel(z)
        R(i) = 1 + z(i) * (b * ((I - z(i)*A) \ one));
    end
    %1 ./ (1 - z);  % Example: Backward Euler
    mask = abs(real(abs(R)))< 1;
    h = imagesc(x(1,:), y(:,1), abs(R));
    axis xy; colorbar; colormap jet;
    caxis([0 1]); % <<< LIMIT THE COLOR RANGE
    hold on;
    contour(x, y, abs(R), [1 1], 'k', 'LineWidth', 2); % stability boundary
    % Apply the mask
    set(h, 'AlphaData', mask);  % Transparent where mask is false
    xlabel('Re'); ylabel('Im');
    title(['Stability Region of ', name],'FontSize', FS);
    set(gca,'FontSize',FS);
    saveas(gcf,['/Users/alexanderweiss/Documents/Papers/Discretization/Images/StabRegion',discMeth,'.jpg']);
end