function h = ebar_plot(data, index, peak, m_size)
% use the native matlab colors
% stock_colors = [0    0.4470    0.7410;
%     0.8500    0.3250    0.0980;
%     0.9290    0.6940    0.1250;
%     0.4940    0.1840    0.5560;
%     0.4660    0.6740    0.1880;
%     0.3010    0.7450    0.9330;
%     0.6350    0.0780    0.1840];
% stock_colors = [0    0.4470    0.7410;
%     0.8500    0.3250    0.0980;
%     0.9290    0.6940    0.1250;
%     0    0.4470    0.7410;
%     0.8500    0.3250    0.0980;
%     0.9290    0.6940    0.1250;
%     0    0.4470    0.7410;
%     0.8500    0.3250    0.0980;
%     0.9290    0.6940    0.1250];

stock_colors = {'r','b'};
line_opt = {'s','s','o','o','p','p'};
color_index = mod(index-1,2)+1;
mu = data(1);
sigma = data(2);
offset = 0.1;
x = [index-offset index+offset index+offset index-offset];
y = [mu-sigma mu-sigma mu+sigma mu+sigma];
C = stock_colors{color_index};
h = patch(x,y,C); hold on;
% plot([index index index], [mu-sigma mu mu+sigma], '+-', ...
%    'LineWidth', l_size, 'MarkerSize', m_size);hold on;
plot(index,mu,'w.','MarkerSize',m_size+10);
plot(index,peak(2), [line_opt{index} stock_colors{color_index}], ...
    'MarkerSize', m_size, 'MarkerFaceColor', stock_colors{color_index});
plot(index, peak(1), [line_opt{index} stock_colors{color_index}], ...
    'MarkerSize', m_size, 'MarkerFaceColor', stock_colors{color_index});
yticks([0:0.1:1]);
end