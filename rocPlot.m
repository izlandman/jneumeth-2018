% pull out ploting function from compute_eer

function eer = rocPlot(FNR, FPR, labels, r_title)

stock_colors = [0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250];
line_opt = {'-','--',':','-','--',':'};
n_line_opts = numel(line_opt);
l_size = 4;
m_size = 4;
f_size = 30;

% enable variable sizes of FNR and FPR to be given
[r_FNR, c_FNR] = size(FNR);
[r_FPR, c_FPR] = size(FPR);

% save eers
eer = zeros(r_FNR,1);

if( r_FNR ~= r_FPR )
    fprintf('Number of rows does not match!\n');
    return;
end

fnr = zeros(numel(FNR{1}),r_FNR);
fpr = fnr;
if( r_FNR > 6 )
    types = 3;
else
    types = 2;
end
for i=1:r_FNR
    % find eer
    difs = FNR{i} - FPR{i};
    idx1 = find(difs < 0, 1, 'last');
    idx2 = find(difs>= 0, 1 );
    x = [FNR{i}(idx1); FPR{i}(idx1)];
    y = [FNR{i}(idx2); FPR{i}(idx2)];
    a = ( x(1) - x(2) ) / ( y(2) - x(2) - y(1) + x(1) );
    eer(i) = 100 * ( x(1) + a * ( y(1) - x(1) ) );
    
    % produce plot
    fnr(:,i) = icdf(FNR{i});
    fpr(:,i) = icdf(FPR{i});
    % handle plotting
    % line_flag = mod(i-1,n_line_opts)+1;
    line_flag = floor( (i-1) / types ) + 1;
    h = plot(fpr(:,i), fnr(:,i),line_opt{line_flag},'LineWidth',l_size);
    h.Color = stock_colors(1+mod(i-1,3),:);
    hold on;
end

min_lim = min( min(fnr(:)), min(fpr(:)) );
max_lim = max( max(fnr(:)), max(fpr(:)) );

% plot EER line?
line([min_lim, max_lim],[min_lim, max_lim],'color','black');


xtick = [0.001, 0.01, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.70];
xticklabel = num2str(xtick * 100, '%g\n');
xticklabel = textscan(xticklabel, '%s'); xticklabel = xticklabel{1};
set (gca, 'xtick', icdf(xtick));
set (gca, 'xticklabel', xticklabel);
x_limits = icdf([0.0001 0.90]);
xlim( x_limits );
% xlim( [min_lim max_lim] );
xlabel ('False Positive Rate (%)');
xtickangle(-45);

ytick = xtick;
yticklabel = num2str(ytick * 100, '%g\n');
yticklabel = textscan(yticklabel, '%s'); yticklabel = yticklabel{1};
set (gca, 'ytick', icdf(ytick));
set (gca, 'yticklabel', yticklabel);
y_limits = icdf([0.0001 0.90]);
ylim( y_limits );
% ylim( [min_lim max_lim] );
ylabel ('False Negative Rate (%)')

legend(labels);
title(r_title);
grid on;
box on;
axis square;
axis manual;
set(gca, 'FontSize', f_size);

end

function y = icdf(x)
% computes the inverse of cumulative distribution function in x
y = -sqrt(2).*erfcinv(2 * ( x + eps));
y(isinf(y)) = nan;
end