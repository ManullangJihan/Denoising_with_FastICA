
function Plot_Signal(x, y, oname2, vlgn, lgn, xlbl, ylbl, tlt, ord, closefig, txt, x1, y1, str)

figure;
plot(x, y, 'LineWidth', 0.5);
pbaspect([3.8 2 1]);
set(gca, 'XColor', 'black', 'YColor', 'black', 'LineWidth', 1.2, 'GridAlpha', 0.1);
set(gca, 'FontSize', 11, 'FontName', 'Times New Roman', 'FontWeight', 'bold');

if txt
    text(x1, y1, str, 'Color','black','FontSize', 12, 'LineStyle' , '-', 'BackgroundColor', 'white', 'EdgeColor', 'black');
end

if ord
    set(ax2, 'YLim', [-1.5, 1.5]);
end

grid on

% xy label
xlabel(xlbl, 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
ylabel(ylbl, 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');

% title label
title(tlt, 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');

% legend
if vlgn
    alg = legend(lgn, 'FontName', 'Rockwell', 'FontSize', 8, 'FontWeight', 'normal', 'Location', 'northeast');
    set(alg, 'LineWidth', 0.1, 'Color', 'white');
end

export_fig (oname2, '-jpg', '-r200', '-a4', '-painters', '-transparent');

if closefig
    close;
end

warning('off')
