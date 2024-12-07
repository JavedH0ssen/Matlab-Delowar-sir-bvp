%% Final Shortened Complex one
clc;
tic;
aa =1;
bb = 2;
cc = 3;
print_state=0;
sol = jobs(aa);
sol1 = jobs(bb);
sol2 = jobs(cc);

name = '\Gamma';
ame=eraseBetween(name,1,1);
% main(aa, bb, cc, name,sol,sol1,sol2)

% function main(aa, bb, cc, name,sol,sol1,sol2)
   close all;
     fig=figure(1);
%     fig = figure('Menu','none','ToolBar','none'); 
%     ah = axes('Units','Normalize','Position',[0 0 1 1]);
    hold on
    
    plotting(2, "f'", sol, sol1, sol2, aa, bb, cc, name)
    

    
    xlim([0 5])

    addArrowAnnotations(fig);
    ax(1) = gca;
    hold off
    fig=figure(2);

    hold on
    plotting(4, "g", sol, sol1, sol2, aa, bb, cc, name)
    xlim([0 7.2])
    addArrowAnnotations(fig);
    legend('Location', 'southeast');
    ax(2) = gca;
    hold off
%     subplot(3, 2, 3)
    fig=figure(3);

    hold on
    plotting(6, "\Xi", sol, sol1, sol2, aa, bb, cc, name)
    xh = get(gca,'ylabel') % handle to the label object
    p = get(xh,'position') % get the current position property
    p(2) = 2*p(2) ;        % double the distance, 
                       % negative values put the label below the axis
    set(xh,'position',p)   % set the new position
    xlim([0 5.5])
    addArrowAnnotations(fig);
    legend('Location', 'southeast');
    ax(3) = gca;
    hold off
%     subplot(3, 2, 5)
    fig=figure(4);

    hold on
    plotting(8, "\Theta", sol, sol1, sol2, aa, bb, cc, name)
    xlim([0 3.2])
    addArrowAnnotations(fig);
    ax(4) = gca;
    hold off
%     subplot(3, 2, 6)
    fig=figure(5);

    hold on
    plotting(10, "\Phi", sol, sol1, sol2, aa, bb, cc, name)
    xlim([0 3.9])
    addArrowAnnotations(fig);
    ax(5) = gca;
    hold off
  if print_state~=0  
    for x =1:5
       exportgraphics(ax(x),append(ame,string(x),".jpg")) 
    end
  end
   
% end
toc
function plotting(line, titles, sol, sol1, sol2, aa, bb, cc, name)
    hold on          
    a = plot(sol.x, sol.y(line, :), '-r');
    b = plot(sol1.x, sol1.y(line, :), '--g');
    c = plot(sol2.x, sol2.y(line, :), '-.b');
    lgd=legend([a; b; c], {[name '= ' num2str(aa)],[name '= ' num2str(bb)],[name '= ' num2str(cc)]}, 'Location', 'northeast');
    lgd.FontSize = 11.5;
    data = [sol.y(line, :), sol1.y(line, :), sol2.y(line, :)];
    minValue = min(data);
    maxValue = max(data);
    ylim([minValue, maxValue+(maxValue-minValue)/100])
    % title(titles)
    % xlabel(['eta(\eta) -------{\bf\rightarrow}'], 'Interpreter', 'tex')
    % ylabel([titles, '------{\bf\rightarrow}'], 'Interpreter', 'tex')
    ax = gca;  % Get the current axes
    ax.FontSize = 10;  % Set the font size to 10
    xlabel('\eta','FontWeight','bold','FontSize', 13, 'Interpreter', 'tex')
    ylabel([titles],'FontWeight','bold','FontSize', 13, 'Interpreter', 'tex')
%     RC = .1; FW = .5;
%     M = 0.5; BI = 1.5; BE = 1.5;R = 0.6;EC = 0.3; GR = 8; GAMMA = 1;  AE = 1 + BE * BI; 
%     GM = 6;  S = 0.5; S0 = 1; K = .5; SC = .6;PR=.71;
%     str = sprintf(['Pr = %.1f, \\betai = %.1f\n\\betae = %.1f,...' ...
%         ' M = %.1f\nR = %.1f, Ec = %.1f,\n\\Gamma = %.1f,Gr=%.1f']...
%         , PR, BI, BE, M, R, EC, GAMMA,GR);
%     text(0, 0, str, 'FontSize', 7, 'Interpreter', 'tex');  % Adjust the coordinates as needed
    hYLabel = get(gca,'YLabel');
    set(hYLabel,'Units','normalized');
    set(hYLabel,'rotation',0)


%     num=sol.y(line, :);
%     s=  char( extractAfter( string(num),'.') );
%     length(s)
    label_n=[-0.02 0 0];
    if (abs(minValue-maxValue)<.7)
        label_n=[0 0 0];
    end
    set(hYLabel,'Position',get(hYLabel,'Position') + label_n);
    
    box on
    hold off

end

function out = jobs(VARABLE)
    RC = .1; FW = .5;
    M = .1; BI = 1.5; BE = 1.5;R =.6;EC =.3; GR = 8; GAMMA = VARABLE;  AE = 1 + BE * BI; 
    GM = 6; PR = .71  ; S0 = 1; K = .5; SC = .6;
    D = 1.5;EPS =.8; ST = .5; STT = .5;
    NEBLA = 2; LAMB = .8; short = AE ^ 2 + BE ^ 2;
    %options = bvpset('RelTol',0.0000001,'Stats','on')
    %     x=sol.x;
    %     y=sol.y;
    %     sol2=bvpinit(sol,[0 6]);
    %     sol=bvp5c(@bvp2d,@bc2d,sol2);
    sol1 = bvpinit(linspace(0,9, 500), [0 0 0 0 0 0 0 0 0 0 0]);
    sol = bvp4c(@bvp2d, @bc2d, sol1);
    sol = bvp4c(@bvp2d, @bc2d, sol);
    out = sol;

    function yvector = bvp2d(~, y)

        yy3 =- (D * y(7) + GR * y(8) + GM * y(10) + y(1) * y(3) / (EPS ^ 2) - ...
            (K + M * AE / short) * y(2) + (R - (M * BE / short)) * y(4) - GAMMA * y(2) ^ 2) / ((1 + D) / EPS);
        yy5 =- (y(1) * y(5) / EPS ^ 2 - (K + M * AE / short) * y(4) - ...
            (R - M * BE / short) * y(2) - GAMMA * y(4) ^ 2) / ((1 + D) / EPS);
        yy7 =- (y(6) * y(2) + y(1) * y(7) - 2 * LAMB * y(6) - LAMB * y(3)) / NEBLA;
        Holder =- ((1 + D) * PR * EC * (y(3) ^ 2 + y(5) ^ 2) + EC * M * PR * (y(2) ^ 2 + ...
            y(4) ^ 2) / short + PR * y(1) * y(9) - ST * PR * y(2));
        yy9 = Holder;
        yy11 =- (y(1) * SC * y(11) - STT * SC * y(2) + S0 * SC * Holder - RC * SC * y(10));
        yvector = [y(2); y(3); yy3; y(5); yy5; y(7); yy7; y(9); yy9; y(11); yy11];
    end

    function residual = bc2d(y0, yinf)

        residual = [y0(2) - 1; y0(1) - FW; y0(4); y0(6) + y0(3) / 2; y0(8) - 1 + ST / 2; ...
                      y0(10) - 1 + STT / 2; ...
                      yinf(2); yinf(4); yinf(6); yinf(8); yinf(10)];
    end

end
function addArrowAnnotations(figure1)
annotation(figure1,'arrow',[0.217327459618208 0.499265785609398],...
    [0.03 0.03]);

% Create arrow
annotation(figure1,'arrow',[0.06 0.06],...
    [0.158663865546219 0.514285714285714]);
end