%% Final Shortened Complex one
clc; clear;
% sol=jobs(2,200);
% sol=jobs(2,50);
% sol=jobs(2,100);
aa = .1;
bb = .5;
cc = 1;
name = 'Fw';
main(aa, bb, cc, name)

function main(aa, bb, cc, name)
    sol = jobs(aa);
    sol1 = jobs(bb);
    sol2 = jobs(cc);
    figure1 = figure('Name', name);
    hold on
    subplot(3, 2, 1)
    plotting(2, "f'", sol, sol1, sol2, aa, bb, cc, name)
    xlim([0 4.5])
    subplot(3, 2, 2)
    plotting(4, "g", sol, sol1, sol2, aa, bb, cc, name)
    subplot(3, 2, 3)
    plotting(6, "\Xi", sol, sol1, sol2, aa, bb, cc, name)
    xlim([0 5.5])
    subplot(3, 2, 5)
    plotting(8, "\Theta", sol, sol1, sol2, aa, bb, cc, name)
    xlim([0 3.2])
    subplot(3, 2, 6)
    plotting(10, "\Phi", sol, sol1, sol2, aa, bb, cc, name)
    xlim([0 3.2])
    hold off
    % Call the function to add arrow annotations
    addArrowAnnotations(figure1);
end

function plotting(line, titles, sol, sol1, sol2, aa, bb, cc, name)
    hold on          
    a = plot(sol.x, sol.y(line, :), '-r');
    b = plot(sol1.x, sol1.y(line, :), '--g');
    c = plot(sol2.x, sol2.y(line, :), '-.b');
    % title(titles)
    % xlabel(['eta(\eta) -------{\bf\rightarrow}'], 'Interpreter', 'tex')
    % ylabel([titles, '------{\bf\rightarrow}'], 'Interpreter', 'tex')
    ax = gca;  % Get the current axes
    ax.FontSize = 5;  % Set the font size to 10
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
    set(hYLabel,'Position',get(hYLabel,'Position') + [-.0401 0.27 0]);
 
    box on
    hold off

end

function out = jobs(VARABLE)
    RC = .1; FW = .5;
    M = VARABLE; BI = 1.5; BE = 1.5;R = 0.6;EC = 0.3; GR = 8; GAMMA = 1;  AE = 1 + BE * BI; 
    GM = 6; PR = .71  ;S = 0.5; S0 = 1; K = .5; SC = .6;
    D = 2;EPS =.8; ST = .5; STT = .5;
    NEBLA = 2; LAMB = .8; short = AE ^ 2 + BE ^ 2; XX = 1 + D;
    %options = bvpset('RelTol',0.0000001,'Stats','on')
    %     x=sol.x;
    %     y=sol.y;
    %     sol2=bvpinit(sol,[0 6]);
    %     sol=bvp5c(@bvp2d,@bc2d,sol2);
    sol1 = bvpinit(linspace(0, 6, 70), [0 0 0 0 0 0 0 0 0 0 0]);
    sol = bvp4c(@bvp2d, @bc2d, sol1);
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

% Create arrow
annotation(figure1,'arrow',[0.152610441767068 0.268983863844123],...
    [0.367067226890757 0.368067226890757]);

% Create arrow
annotation(figure1,'arrow',[0.152610441767068 0.268983863844123],...
    [0.678991596638655 0.678991596638656]);

% Create arrow
annotation(figure1,'arrow',[0.152610441767068 0.268983863844123],...
    [0.0679075630252105 0.0689075630252107]);

% Create arrow
annotation(figure1,'arrow',[0.593038821954485 0.711678769639513],...
    [0.667226890756303 0.667226890756304]);

% Create arrow
annotation(figure1,'arrow',[0.593038821954485 0.711678769639513],...
    [0.0679075630252104 0.0689075630252106]);

% Create arrow
annotation(figure1,'arrow',[0.0933503836317142 0.0933503836317142],...
    [0.741857142857146 0.862184873949583]);

% Create arrow
annotation(figure1,'arrow',[0.0893343193746859 0.0893343193746859],...
    [0.449420168067228 0.569747899159665]);

% Create arrow
annotation(figure1,'arrow',[0.523391222641952 0.523391222641952],...
    [0.7317731092437 0.852100840336137]);

% Create arrow
annotation(figure1,'arrow',[0.0891130954961223 0.0891130954961223],...
    [0.143537815126052 0.263865546218489]);

% Create arrow
annotation(figure1,'arrow',[0.536722786648495 0.536722786648495],...
    [0.141857142857144 0.262184873949581]);
end
