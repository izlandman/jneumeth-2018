% plot data from 10s, 5s, 2s, and 1s epoch based results
function epochVccr(epoch_10, epoch_5, epoch_2, epoch_1)

% collect data from each experiment
% data, [x y], epoch
ccr_data = cell(3,2,4);
ccr_sub = cell(3,2,4);
eer_data = ccr_data;
eer_sub = ccr_data;
labels = cell(3,4);
% cells! (1) sub_ses_data (2) sub_data (3) ses_data
[epoch_10_m, epoch_10_g, epoch_10_i ]= expTriFolder(epoch_10,0);
[epoch_5_m, epoch_5_g, epoch_5_i ]= expTriFolder(epoch_5,0);
[epoch_2_m, epoch_2_g, epoch_2_i ]= expTriFolder(epoch_2,0);
[epoch_1_m, epoch_1_g, epoch_1_i ]= expTriFolder(epoch_1,0);

% error bar plot data!
eM = uniteEpochs(2,epoch_10_m{1},epoch_5_m{1},epoch_2_m{1},epoch_1_m{1})/100;
eG = uniteEpochs(2,epoch_10_g{1},epoch_5_g{1},epoch_2_g{1},epoch_1_g{1})/100;
eI = uniteEpochs(2,epoch_10_i{1},epoch_5_i{1},epoch_2_i{1},epoch_1_i{1})/100;
cM = uniteEpochs(1,epoch_10_m{1},epoch_5_m{1},epoch_2_m{1},epoch_1_m{1});
cG = uniteEpochs(1,epoch_10_g{1},epoch_5_g{1},epoch_2_g{1},epoch_1_g{1});
cI = uniteEpochs(1,epoch_10_i{1},epoch_5_i{1},epoch_2_i{1},epoch_1_i{1});
%  average subject based classification

% figure('Name','Stretched Plot','NumberTitle','off');
% grid on; hold on;
% err_m = ebarPlot(eM,2);
% err_g = ebarPlot(eG,2);
% err_i = ebarPlot(eI,2);
% 
% figure('Name','Stretched Sub-Plots','NumberTitle','off');
% grid on; hold on;
% subplot(311);
% ebarPlot(eM,2);
% subplot(312);
% ebarPlot(eG,2);
% subplot(313);
% ebarPlot(eI,2);

% this is neat
l_size = 2;
f_size = 30;
e_max = max(max(cat(1,eM(1,:),eG(1,:),eI(1,:))))+...
    max(max(cat(1,eM(2,:),eG(2,:),eI(2,:))));
e_max = 0.05*ceil(e_max/0.05);
figure('Name',[epoch_10 ' | Stretch Sub-Plot of the Future'],...
    'NumberTitle','off','units','normalized','outerposition',[0  0 1 1]);
subplot(311);grid on; title('Mahalanobis Performance');
errorYY(eM,cM,l_size,f_size,e_max);
subplot(312);grid on; title('GMM-UBM Performance');
errorYY(eG,cG,l_size,f_size,e_max);
subplot(313);grid on; title('I-Vector Performance');
errorYY(eI,cI,l_size,f_size,e_max);
% add all sorts of labels here!
h0 = findobj(gcf,'type','axes');
p = get(gca,'position');
h_top = axes('position',[p(1) p(2)+p(4) p(3) 0]);
h_top.XAxisLocation = 'top';
h_top.XLim = h0(1).XLim;
h_top.XTick = [1:1:9];
h_top.XTickLabel = ({'2','','','','UBM','','','','512'});
h_top.FontSize = f_size-4;
% h_bot.XLabel('Epoch Size (seconds)');
h_bot.FontSize = f_size;
h_bot = axes('position', [p(1:3) 0]);
h_bot.XLim = h0(1).XLim;
h_bot.XTick = [5:9:36];
h_bot.XTickLabel = [10, 5, 2, 1];
h_bot.TickLength(1) = 0;
xlabel('Epoch Size (seconds)');
h_bot.FontSize = f_size;

% l_size = 3;
% figure('name','EER','NumberTitle','off');
% hold on; grid on;
% ebarPlotPrep(epoch_10_g{1},l_size);
% ebarPlotPrep(epoch_10_i{1},l_size);
% ebarPlotPrep(epoch_10_m{1},l_size);
% ebarPlotPrep(epoch_5_g{1},l_size);
% ebarPlotPrep(epoch_5_i{1},l_size);
% ebarPlotPrep(epoch_5_m{1},l_size);
% ebarPlotPrep(epoch_2_g{1},l_size);
% ebarPlotPrep(epoch_2_i{1},l_size);
% ebarPlotPrep(epoch_2_m{1},l_size);
% ebarPlotPrep(epoch_1_g{1},l_size);
% ebarPlotPrep(epoch_1_i{1},l_size);
% ebarPlotPrep(epoch_1_m{1},l_size);

% Subject-Session Performance
%% 10s
name = '10s';
[ccr_data(:,:,1), labels(:,1)] = extractSubSes(epoch_10_m{1},...
    epoch_10_g{1},epoch_10_i{1},name,1);
[eer_data(:,:,1), labels(:,1)] = extractSubSes(epoch_10_m{1},...
    epoch_10_g{1},epoch_10_i{1},name,2);
[ccr_sub(:,:,1),eer_sub(:,:,1)] = extractSub(epoch_10_m{2},...
    epoch_10_g{2},epoch_10_i{2},name);
%% 5s
name = '5s';
[ccr_data(:,:,2), labels(:,2)] = extractSubSes(epoch_5_m{1},...
    epoch_5_g{1},epoch_5_i{1},name,1);
[eer_data(:,:,2), labels(:,2)] = extractSubSes(epoch_5_m{1},...
    epoch_5_g{1},epoch_5_i{1},name,2);
[ccr_sub(:,:,2),eer_sub(:,:,2)] = extractSub(epoch_5_m{2},...
    epoch_5_g{2},epoch_5_i{2},name);
%%  2s
name = '2s';
[ccr_data(:,:,3), labels(:,3)] = extractSubSes(epoch_2_m{1},...
    epoch_2_g{1},epoch_2_i{1},name,1);
[eer_data(:,:,3), labels(:,3)] = extractSubSes(epoch_2_m{1},...
    epoch_2_g{1},epoch_2_i{1},name,2);
[ccr_sub(:,:,3),eer_sub(:,:,3)] = extractSub(epoch_2_m{2},...
    epoch_2_g{2},epoch_2_i{2},name);
%% 1s
name = '1s';
[ccr_data(:,:,4), labels(:,4)] = extractSubSes(epoch_1_m{1},...
    epoch_1_g{1},epoch_1_i{1},name,1);
[eer_data(:,:,4), labels(:,4)] = extractSubSes(epoch_1_m{1},...
    epoch_1_g{1},epoch_1_i{1},name,2);
[ccr_sub(:,:,4),eer_sub(:,:,4)] = extractSub(epoch_1_m{2},...
    epoch_1_g{2},epoch_1_i{2},name);

% build newScore data!
mixtures = numel(ccr_data{end,end,end});
time_steps = size(ccr_data,3);
mahal_newScores = zeros(1,2,time_steps);
gmm_newScores = zeros(mixtures,2,time_steps);
ivector_newScores = zeros(mixtures,2,time_steps);
m_plot_scores = zeros(1,2,time_steps);
g_plot_scores = zeros(mixtures,2,time_steps);
i_plot_scores = zeros(mixtures,2,time_steps);

for i=1:time_steps
    mahal_newScores(:,1,i) = newScore(ccr_data{1,1,i},ccr_sub{1,1,i});
    mahal_newScores(:,2,i) = newScore(ccr_sub{1,1,i}',eer_sub{1,1,i}');
    m_plot_scores(:,1,i) = newScore(ccr_data{1,1,i},ccr_data{1,1,i});
    gmm_newScores(:,1,i) = newScore(ccr_data{2,1,i},ccr_sub{2,1,i});
    gmm_newScores(:,2,i) = newScore(ccr_sub{2,1,i},eer_sub{2,1,i});
    g_plot_scores(:,1,i) = newScore(ccr_data{2,1,i},ccr_data{2,1,i});
    ivector_newScores(:,1,i) = newScore(ccr_data{3,1,i},ccr_sub{3,1,i});
    ivector_newScores(:,2,i) = newScore(ccr_sub{3,1,i},eer_sub{3,1,i});
    i_plot_scores(:,1,i) = newScore(ccr_data{3,1,i},ccr_data{3,1,i});
end

% add additional plot focusing on the higher mixtures and combine onto a
% single plot
m_size_diss = 25;
l_size_diss = 1.5;
f_size_diss = 30;
line_sty = {'s:' , 'o:', 'p:'};
% marker_sty = {'o', 's', 'p'};
color_sty = {'r','b'};
m_range = [5 9];

figure('numbertitle','off','name','Grouped Results');
% first plot
x1_min = 0.09;
y1_min = 0.55;
x_axis = 0.95;
y_axis = 0.935;
axes('Position',[x1_min y1_min x_axis-x1_min y_axis-y1_min]);
title('CRR and EER verus Epoch and Mixture');
hold on;
ylabel('CRR');
plot(epoch_2_m{1}(:,1,1),[ color_sty{2} line_sty{1}],'LineWidth',...
    l_size_diss,'MarkerSize',m_size_diss);
plot(epoch_2_g{1}(:,1,1),[ color_sty{2} line_sty{2}],'LineWidth',...
    l_size_diss,'MarkerSize',m_size_diss);
plot(epoch_2_i{1}(:,1,1),[ color_sty{2} line_sty{3}],'LineWidth',...
    l_size_diss,'MarkerSize',m_size_diss);
plot(epoch_1_m{1}(:,1,1),[ color_sty{1} line_sty{1}],'LineWidth',...
    l_size_diss,'MarkerSize',m_size_diss);
plot(epoch_1_g{1}(:,1,1),[ color_sty{1} line_sty{2}],'LineWidth',...
    l_size_diss,'MarkerSize',m_size_diss);
plot(epoch_1_i{1}(:,1,1),[ color_sty{1} line_sty{3}],'LineWidth',...
    l_size_diss,'MarkerSize',m_size_diss);
xticklabels({' ',' ',' ',' ',' ',' ',' ',' ',' '});
xlim(m_range);
y_bot = min([min(epoch_1_m{1}(m_range(1):m_range(2),1,1)),...
    min(epoch_1_g{1}(m_range(1):m_range(2),1,1)),...
    min(epoch_1_i{1}(m_range(1):m_range(2),1,1)),...
    min(epoch_2_m{1}(m_range(1):m_range(2),1,1)),...
    min(epoch_2_g{1}(m_range(1):m_range(2),1,1)),...
    min(epoch_2_i{1}(m_range(1):m_range(2),1,1))]);
y_top = max([max(epoch_1_m{1}(m_range(1):m_range(2),1,1)),...
    max(epoch_1_g{1}(m_range(1):m_range(2),1,1)),...
    max(epoch_1_i{1}(m_range(1):m_range(2),1,1)),...
    max(epoch_2_m{1}(m_range(1):m_range(2),1,1)),...
    max(epoch_2_g{1}(m_range(1):m_range(2),1,1)),...
    max(epoch_2_i{1}(m_range(1):m_range(2),1,1))]);
ylim([y_bot*.98 y_top*1.01]);
grid on;
set(gca,'fontsize',f_size_diss);

% second plot
y_axis = y1_min-0.02;
y1_min = 0.03;
axes('Position',[x1_min y1_min x_axis-x1_min y_axis-y1_min]);
hold on;
grid on;
ylabel('EER');
plot(epoch_2_m{1}(:,2,1)/100,[ color_sty{2} line_sty{1}],'LineWidth',...
    l_size_diss,'MarkerSize',m_size_diss);
plot(epoch_2_g{1}(:,2,1)/100,[ color_sty{2} line_sty{2}],'LineWidth',...
    l_size_diss,'MarkerSize',m_size_diss);
plot(epoch_2_i{1}(:,2,1)/100,[ color_sty{2} line_sty{3}],'LineWidth',...
    l_size_diss,'MarkerSize',m_size_diss);
plot(epoch_1_m{1}(:,2,1)/100,[ color_sty{1} line_sty{1}],'LineWidth',...
    l_size_diss,'MarkerSize',m_size_diss);
plot(epoch_1_g{1}(:,2,1)/100,[ color_sty{1} line_sty{2}],'LineWidth',...
    l_size_diss,'MarkerSize',m_size_diss);
plot(epoch_1_i{1}(:,2,1)/100,[ color_sty{1} line_sty{3}],'LineWidth',...
    l_size_diss,'MarkerSize',m_size_diss);
% labeling
xticklabels({'2','4','8','16','32','64','128','256','512'});
xlim(m_range);
y_bot = min([min(epoch_1_m{1}(m_range(1):m_range(2),2,1)),...
    min(epoch_1_g{1}(m_range(1):m_range(2),2,1)),...
    min(epoch_1_i{1}(m_range(1):m_range(2),2,1)),...
    min(epoch_2_m{1}(m_range(1):m_range(2),2,1)),...
    min(epoch_2_g{1}(m_range(1):m_range(2),2,1)),...
    min(epoch_2_i{1}(m_range(1):m_range(2),2,1))]);
y_top = max([max(epoch_1_m{1}(m_range(1):m_range(2),1,1)),...
    max(epoch_1_g{1}(m_range(1):m_range(2),2,1)),...
    max(epoch_1_i{1}(m_range(1):m_range(2),2,1)),...
    max(epoch_2_m{1}(m_range(1):m_range(2),2,1)),...
    max(epoch_2_g{1}(m_range(1):m_range(2),2,1)),...
    max(epoch_2_i{1}(m_range(1):m_range(2),2,1))]);
ylim([y_bot*.98/100 y_top*1.01/100]);
legend({'Mahal 2s','GMM-UBM 2s','IVEC 2s','Mahal 1s', 'GMM-UBM 1s','IVEC 1s'},...
    'orientation','horizontal','location','southoutside');
grid on;
xlabel('Mixture size');
set(gca,'fontsize',f_size_diss);
%% if plot_data is built, run a function to make the plots!
% this allows CCR and EER to be buil in the same format since the data is
% mirrored

% % Pure Subject Matching
% doublePlotCcr(ccr_sub,labels,'CCR');
% doublePlotEer(eer_sub,labels,'EER');
% 
% % CCR
% doublePlotCcr(ccr_data,labels,'CCR');
% 
% % EER
% doublePlotEer(eer_data,labels,'EER');
% 
% % circles?
% circlePlot(m_plot_scores,g_plot_scores,i_plot_scores);
end

% plot error bar data of EER along with CCR on secondary axis!
function errorYY(e_data,c_data,l_size,f_size,e_max)

% left axis
yyaxis left
e = ebarPlot(e_data,l_size);
ylim([0 e_max]);ylabel('EER');
yticks([0:(e_max)/5:e_max*1.1]);

% right axis
yyaxis right
e = ebarPlot(c_data,l_size);
ylim([0 1]);ylabel('CCR');
yticks([0,0.2,0.4,0.6,0.8,1]);

xlim([0 36+1]);
xticks([1,10,19,28]);
xticklabels({'10','5','2','1'});
set(gca,'FontSize',f_size);

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

set(gca,'xticklabel',{[]});
end

function result = uniteEpochs(dim,e10,e5,e2,e1)
[mix,~,~] = size(e10);
result = zeros(2,mix*4);

for i=1:2
    result(i,:) = cat(1,e10(:,dim,i),e5(:,dim,i),e2(:,dim,i),e1(:,dim,i));
end

end

function ebarPlotPrep(data,l_size)

result = [data(:,2,1),data(:,2,2)]';
ebarPlot(result,l_size);
end

function e = ebarPlot(data,l_size)
e = errorbar(data(1,:),data(2,:));
e.LineStyle = 'none';
e.Marker = 's';
e.LineWidth = l_size;
end

function doublePlotEer(plot_data,labels,plot_name)

% plot some data!
labels_list = labels(:);
iterations = numel(labels_list);
line_type = {'-','--','-.',':'};
marker_type = {'*','sq','o'};
l_width = 2;
m_size = 12;
f_size = 18;
x1_min = 0.05;
y1_min = 0.16;
x_axis = 0.95;
y_axis = 0.95;

% plot it
figure('name',plot_name,'numbertitle','off');
axes('Position',[x1_min y1_min x_axis-x1_min y_axis-y1_min]);

for i=1:iterations
    % organize averaged subject CCR epoch_xx_x{1}
    index = floor( (i-1)/3 ) + 1;
    step = mod( (i-1),3 ) + 1;
    line_style = [line_type{index} marker_type{step}];
    x_data = plot_data{step,2,index};
    y_data = plot_data{step,1,index};
    plot(x_data,y_data,line_style,'LineWidth',l_width,...
        'MarkerSize',m_size);
    hold on;
end
xticklabels(int2str(2.^(plot_data{end,end,end})));
xlabel('UBM Mixture Size');
ylabel(plot_name);
box off;
grid on;
title([plot_name ' as a function of epoch and UBM mixture size']);
set(gca,'FontSize',f_size);
Leg = legend(labels_list,'Orientation', 'horizontal');
Leg.Box = 'off';
set(Leg, 'position', [ 0.01 0.01 0.99 0.07]);
end

function doublePlotCcr(plot_data,labels,plot_name)

% plot some data!
labels_list = labels(:);
iterations = numel(labels_list);
line_type = {'-','--','-.',':'};
marker_type = {'*','sq','o'};
l_width = 2;
m_size = 12;
f_size = 18;
x1_min = 0.05;
y1_min = 0.16;
x_axis = 0.70;
y_axis = 0.95;

% plot it
figure('name',plot_name,'numbertitle','off');
axes('Position',[x1_min y1_min x_axis-x1_min y_axis-y1_min]);

for i=1:iterations
    % organize averaged subject CCR epoch_xx_x{1}
    index = floor( (i-1)/3 ) + 1;
    step = mod( (i-1),3 ) + 1;
    line_style = [line_type{index} marker_type{step}];
    x_data = plot_data{step,2,index};
    y_data = plot_data{step,1,index};
    h1(i) = plot(x_data,y_data,line_style,'LineWidth',l_width,...
        'MarkerSize',m_size);
    hold on;
end
xticklabels(int2str(2.^(plot_data{end,end,end})));
xlabel('UBM Mixture Size');
ylabel(plot_name);
box off;
grid on;
title([plot_name ' as a function of epoch and UBM mixture size']);
set(gca,'FontSize',f_size);

% inset plot!
axes('Position', ...
    [x_axis+x1_min y1_min 1-x_axis-x1_min*1.5 y_axis-y1_min]);
for i=1:iterations
    index = floor( (i-1)/3 ) + 1;
    step = mod( (i-1),3 ) + 1;
    line_style = [line_type{index} marker_type{step}];
    x_data = plot_data{step,2,index};
    y_data = plot_data{step,1,index};
    % rescale data?
    new_y_data = plotRescale(y_data,0.50);
    plot(x_data,new_y_data,line_style,'LineWidth',l_width,...
        'MarkerSize',m_size);
    hold on;
end
xlim([1 numel(plot_data{end,end,end})]);
xticks(1:numel(plot_data{end,end,end}));
xticklabels(int2str(2.^(plot_data{end,end,end})));
yticks([ 0.1 0.2 0.3 0.4 0.5 0.6 0.8 1.1 1.5 2.0]);
yticklabels({'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'});
ylim([0.50 1.6]);
grid on;

Leg = legend(h1, labels_list, 'Orientation', 'horizontal');
set(Leg, 'position', [ 0.01 0.01 0.99 0.07]);
Leg.Box = 'off';
set(gca,'FontSize',f_size);

end

function [result,scale] = plotRescale(data, base)
iter = numel(data);
result = data;
summit = 1;
steps = uint8((summit - base)*10);
scale = zeros(steps,1);
for i=1:iter
    step = base;
    scale = 0.10;
    offset = base;
    while step <= summit
        if( data(i) >= step && step+0.1 > data(i) )
            result(i) = offset + (data(i)-step)/0.10*scale;
        end
        offset = offset + scale;
        step = step + 0.10;
        scale = scale + 0.10;
    end
end

end

function [ccr, eer, names] = extractSub(mahal_data,gmm_data,ivector_data,...
    name)
names = cell(3,1);
types = {'Mahal ', 'GMM ', 'I-Vector '};
names(:,1) = strcat(types,name)';
ccr = cell(3,2);
eer = cell(3,2);

ccr{1,1} = repmat(mean(mahal_data(:,1,1)),2,1);
ccr{2,1} = mean(gmm_data(:,:,1))';
ccr{3,1} = mean(ivector_data(:,:,1))';

eer{1,1} = repmat(mean(mahal_data(:,1,2)),2,1);
eer{2,1} = mean(gmm_data(:,:,2))';
eer{3,1} = mean(ivector_data(:,:,2))';

mixtures = numel(ivector_data(1,:,1));
ccr{1,2} = [1 mixtures]';
ccr{2,2} = [1:mixtures]';
ccr{3,2} = [1:mixtures]';
eer(:,2) = ccr(:,2);
end

function [results, names] = extractSubSes(mahal_data,gmm_data,...
    ivector_data,name,type)
names = cell(3,1);
types = {'Mahal ', 'GMM ', 'I-Vector '};
names(:,1) = strcat(types,name)';

results = cell(3,2);
% EER is 1 to 100, CCR is 0.0 to 1.0!
if( type == 2 )
    mahal_data = mahal_data ./ 100;
    gmm_data = gmm_data ./ 100;
    ivector_data = ivector_data ./ 100;
end
results{1,1} = [repmat(mahal_data(type,1,1),1,2)];
results{2,1} = gmm_data(:,type,1);
results{3,1} = ivector_data(:,type,1);

mixtures = numel(ivector_data(:,1,1));
results{1,2} = [1 mixtures];
results{2,2} = [1:mixtures]';
results{3,2} = [1:mixtures]';

end

function [sub_ses, sub, ses] = dataGet(folder)

sub_ses = iVectorBinary([folder filesep 'sub_ses.bin']);
sub = iVectorBinary([folder filesep 'sub.bin']);
ses = iVectorBinary([folder filesep 'ses.bin']);

end