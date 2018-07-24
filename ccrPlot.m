% plot multiple ccr_bench results. provide a vector of folder (experiment)
% numbers and a cell of associated labels for the plot

function ccrPlot(condition)

% statically define the labels. do this by methodology, LaRocca, GMM,
% I-Vector and feature type Cepstrum, PSD, COH

switch condition
    case {'All'}
        labels = {'PSD-EO-FULL' 'PSD-EO-F' 'PSD-EO-C' 'PSD-EO-P' ...
            'PSD-EC-FULL' 'PSD-EC-F' 'PSD-EC-C' 'PSD-EC-P' ...
            'COH-EO-FULL' 'COH-EO-F' 'COH-EO-C' 'COH-EO-P' ...
            'COH-EC-FULL' 'COH-EC-F' 'COH-EC-C' 'COH-EC-P'};
        folders = [0904 0905 0906 0907 0908 0909 0910 0911 0912 0913 ...
            0914 0915 0916 0917 0918 0919];
        x_title = {'All'};
    case {'FigureF'}
        labels = {'PSD-EO-F' 'PSD-EC-F' 'COH-EO-F' 'COH-EC-F'};
        folders = [0905 0909 0913 0917];
        x_title = {'Frontal'};
    case {'FigureC'}
        labels = {'PSD-EO-C' 'PSD-EC-C' 'COH-EO-C' 'COH-EC-C'};
        folders = [0906 0910 0914 0918];
        x_title = {'Central'};
    case {'FigureP'}
        labels = {'PSD-EO-P' 'PSD-EC-P' 'COH-EO-P' 'COH-EC-P'};
        folders = [0907 0911 0915 0919];
        x_title = {'Parietal'};
    case {'Figure1'}
        labels = {'PSD-EO-MHAL' 'PSD-EC-MHAL' 'COH-EO-MHAL' ...
            'COH-EC-MHAL' 'CEP-EO-GMM' 'CEP-EC-GMM' 'CEP-EO-IVEC' ...
            'CEP-EC-IVEC'};
        folders = [0904 0908 0912 0916 930 931 900 901];
        x_title = {'Baseline Performance'};
    case {'Figure2'}
        labels = {'PSD-EO' 'COH-EO' 'CEP-EO'};
        folders = [920 926 930];
        x_title = {'GMM-UBM Classification'};
    case {'Figure2a'}
        labels = {'PSD-EC' 'COH-EC' 'CEP-EC'};
        folders = [921 927 931];
        x_title = {'GMM-UBM Classification'};
    case {'Figure2b'}
        labels = {'PSD-EO' 'PSD-EC' 'COH-EO' 'COH-EC' 'CEP-EO' 'CEP-EC'};
        folders = [920 921 926 927 930 931];
        x_title = {'GMM-UBM Classification'};
    case {'Figure3'}
        labels = {'PSD-EO' 'COH-EO' 'CEP-EO'};
        folders = [922 928 900];
        x_title = {'I-Vector Classification'};
    case {'Figure3a'}
        labels = {'PSD-EC' 'COH-EC' 'CEP-EC'};
        folders = [923 929 901];
        x_title = {'I-Vector Classification'};
    case {'Figure3b'}
        labels = {'PSD-EO' 'PSD-EC' 'COH-EO' 'COH-EC' 'CEP-EO' 'CEP-EC'};
        folders = [922 923 928 929 900 901];
        x_title = {'I-Vector Classification'};
    case {'Figure4'}
        labels = {'PSD-EO' 'COH-EO' 'CEP-EO'};
        folders = [0904 0912 924];
        x_title = {'Mahalanobis Classification'};
    case {'Figure4a'}
        labels = {'PSD-EC' 'COH-EC' 'CEP-EC'};
        folders = [0905 0913 925];
        x_title = {'Mahalanobis Classification'};
    case {'Figure4b'}
        labels = {'PSD-EO' 'PSD-EC' 'COH-EO' 'COH-EC' 'CEP-EO' 'CEP-EC'};
        folders = [904 905 912 913 924 925];
        x_title = {'Mahalanobis Classification'};
    case {'Figure5'}
        labels = {'PSD-MHAL','COH-MHAL','CEP-MHAL','PSD-GMM','COH-GMM',...
            'CEP-GMM','PSD-IVEC','COH-IVEC','CEP-IVEC'};
        folders = [0904 0912 924 920 926 930 922 928 900];
        x_title = {'EO Feature Sets'};
    case {'alt-cov'}
        labels = {'cep-mhal-eo' 'cep-mhal-ec' 'psd-mhal-eo' ...
            'psd-mhal-ec' 'coh-mhal-eo' 'coh-mhal-ec'};
        folders = [934 935 936 937 938 939];
        x_title = {'Alt-COV'};
    otherwise
        fprintf('What is this madness!\n');
        % assume folder of full experiments was given
end



% generate all full titles
p_title = ['mCRR from ' x_title{1}];
e_title = ['EER from ' x_title{1}];
r_title = ['mEER from ' x_title{1}];
n_folders = numel(folders);
n_labels = numel(labels);

if( n_folders ~= n_labels )
    fprintf('n_folders and n_labels are not equal!.\n');
end

% determine mixture size!
mixture_size = 5;

% collect folder data
plot_data = cell(n_folders,5);
ebar_data = zeros(n_folders,4);
FNR = cell(n_folders,1);
FPR = cell(n_folders,1);
for i=1:n_folders
    f_name = ['exp_' num2str(folders(i),'%04d') filesep 'CCR.bin'];
    ccr = iVectorBinary(f_name);
    f_name = ['exp_' num2str(folders(i),'%04d') filesep 'CCR_bench.bin'];
    ccr_bench = iVectorBinary(f_name);
    f_name = ['exp_' num2str(folders(i),'%04d') filesep 'mCCR.bin'];
    mccr = iVectorBinary(f_name);
    f_name = ['exp_' num2str(folders(i),'%04d') filesep 'EER.bin'];
    eer = iVectorBinary(f_name);
    eer = eer/100;
    f_name = ['exp_' num2str(folders(i),'%04d') filesep 'mEER.bin'];
    meer = iVectorBinary(f_name);
    meer = meer/100;
    % collect COR data
    f_name = ['exp_' num2str(folders(i),'%04d') filesep 'mFNR.bin'];
    fnr = iVectorBinary(f_name);
    f_name = ['exp_' num2str(folders(i),'%04d') filesep 'mFPR.bin'];
    fpr = iVectorBinary(f_name);
    
    [channels,elements] = size(mccr);
    if( elements > 1 )
        plot_data{i,2} = mccr(:,mixture_size);
        plot_data{i,4} = meer(:,mixture_size);
    else
        plot_data{i,2} = mccr;
        plot_data{i,4} = meer;
        
    end
    
    [channels,elements] = size(ccr_bench);
    if( elements > 4 )
        ccr_index = (mixture_size-1)*4+1;
        plot_data{i,1} = ccr_bench(ccr_bench(:,ccr_index)~=0,ccr_index);
    else
        plot_data{i,1} = ccr_bench(ccr_bench(:,1)~=0,1);
    end
    
    [~,elements] = size(fnr);
    if( elements > 1 )
        FNR{i} = fnr(:,mixture_size);
        FPR{i} = fpr(:,mixture_size);
    else
        FNR{i} = fnr;
        FPR{i} = fpr;
    end
    
    ebar_data(i,1) = mean(plot_data{i,2});
    ebar_data(i,2) = std(plot_data{i,2});
    ebar_data(i,3) = mean(plot_data{i,4});
    ebar_data(i,4) = std(plot_data{i,4});
    plot_data{i,3} = [max(ccr(:)), min(ccr(:))];
    plot_data{i,5} = [max(eer(:)), min(eer(:))];
end

line_opt = {'rs--','bs--','ro--','bo--','rp--','bp--','d--','+--'};
n_line_opts = numel(line_opt);

figure('numbertitle','off','name',['mCRR Evaluation: ' condition]);
x1_min = 0.07;
y1_min = 0.175;
x_axis = 0.70;
y_axis = 0.935;
f_size = 30;
l_size = 4;
m_size = 15;

axes('Position',[x1_min y1_min x_axis-x1_min y_axis-y1_min]);
for i=1:n_folders
    line_flag = mod(i-1,n_line_opts)+1;
    h1(i) = plot( plot_data{i,1}, line_opt{line_flag}, 'LineWidth', ...
        l_size, 'MarkerSize', m_size);
    hold on; grid on;
end
box off;
ylim([0 1]);
yticks([0:0.1:1]);
xlim([0.5 max(cellfun(@numel,plot_data(:,1)))+0.5]);
title(p_title);
ylabel('CRR (%)');
xlabel('Elements');
set(gca, 'FontSize', f_size);

% inset plot?
axes('Position', ...
    [x_axis+x1_min/10 y1_min 1-x_axis-x1_min/4 y_axis-y1_min]);
for i=1:n_folders
    ebar_plot(ebar_data(i,1:2),i,plot_data{i,3},m_size);
    hold on;
end
set(gca,'xticklabels',[]);
set(gca,'yticklabels',[]);
box off;
grid on;
ylim([0 1]);
xlim([0.75 n_folders+0.25]);
xticks([1:n_folders]);
% xticklabels(labels);
% xtickangle(45);
% ylabel('mCCR (%)');
title('CRR Distribution');
set(gca, 'FontSize', f_size);
% legend?
Leg = legend(h1, labels, 'Orientation', 'horizontal');
set(Leg, 'position', [ 0.01 0.005 0.99 0.07]);
Leg.Box = 'off';
% Leg2 = legend(h2, labels, 'Orientation', 'horizontal');
% set(Leg2, 'position', [ 0.01 0.01 0.99 0.07]);

% eer plot!
figure('numbertitle','off','name',['eer evaluation: ' condition]);
ylabel('mEER (%)');
for i=1:n_folders
    h2(i) = ebar_plot(ebar_data(i,3:4),i,plot_data{i,5},m_size);
    hold on; grid on;
end
ylim([0 1]);
xlim([0.75 n_folders+0.25]);
xticks([1:n_folders]);
set(gca,'xticklabels',[]);
title(e_title);
set(gca, 'FontSize', f_size);
% legend?
Leg2 = legend(h2, labels, 'Orientation', 'horizontal');
set(Leg2, 'position', [ 0.01 0.005 0.99 0.07]);
Leg2.Box = 'off';

% ROC plot
figure('numbertitle','off','name',['ROC: ' condition]);
eers = rocPlot(FNR,FPR,labels,r_title);

end