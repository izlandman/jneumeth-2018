% process data from an experimentTrio results folder
function [mahal_data, gmm_data, ivector_data] = ...
    expTriFolder(folder_name,plot_flag)

mahal_folder = [folder_name filesep 'MAHAL'];
gmm_folder = [folder_name filesep 'GMM'];
ivector_folder = [folder_name filesep 'IVECTOR'];

% acquire raw data
[sub_ses_m, sub_m, ses_m] = dataGet(mahal_folder);
n_session = size(ses_m,2) / 2;
n_subject = size(sub_m,1) / (size(sub_ses_m,1)/2);
[sub_ses_g, sub_g, ses_g] = dataGet(gmm_folder);
[sub_ses_i, sub_i, ses_i] = dataGet(ivector_folder);

% acquire processed data
[SUB_SES_m, SUB_m, SES_m] = singleFolderData(sub_ses_m, sub_m, ses_m, ...
    n_subject, n_session);
[SUB_SES_g, SUB_g, SES_g] = singleFolderData(sub_ses_g, sub_g, ses_g, ...
    n_subject, n_session);
[SUB_SES_i, SUB_i, SES_i] = singleFolderData(sub_ses_i, sub_i, ses_i, ...
    n_subject, n_session);

% scale out Mahal results based upon mixture count
mixtures = size(SUB_SES_g,1);
SUB_SES_m = repmat(SUB_SES_m,mixtures,1);
SES_m = repmat(SES_m,1,1,mixtures);

if( plot_flag > 0 )
    elements = 2.^(1:mixtures);
    x_labels = int2str(elements(:));
    names = {'Mahal','GMM','I-Vector'};
    p_title = ['CCR of ' int2str(plot_flag) 's epochs'];
    l_size = 2;
    m_size = 14;
    % only one of each value for MAHAL, use as baseline for plots
    
    % Overall CCR Performance
    figure('name',[folder_name ' | Subject Average CCR'],'NumberTitle','off');
    grid on;
    hold on;
    this_plot([1 mixtures],[SUB_SES_m(1,1,1) SUB_SES_m(1,1,1)],...
        '-',l_size,m_size);
    this_plot([1:mixtures],SUB_SES_g(:,1,1),'sq-',l_size,m_size);
    this_plot([1:mixtures],SUB_SES_i(:,1,1),'o-',l_size,m_size);
    xticks([1:mixtures]);
    xticklabels(x_labels);
    legend(names);
    title(p_title);
    set(gca,'FontSize',22);
    ylim([0 1]);
    
    % Session Based CCR Performance
    figure('name',[folder_name ' | Sesssion Average CCR'],'NumberTitle','off');
    grid on;
    hold on;
    plot([1 mixtures],[SES_m(1,1) SES_m(1,1)],'--','LineWidth',2,'MarkerSize',14);
    plot([1 mixtures],[SES_m(3,1) SES_m(3,1)],'--','LineWidth',2,'MarkerSize',14);
    plot(squeeze(SES_g(1,1,:)),'sq-','LineWidth',2,'MarkerSize',14);
    plot(squeeze(SES_g(3,1,:)),'sq-','LineWidth',2,'MarkerSize',14);
    plot(squeeze(SES_i(1,1,:)),'o-','LineWidth',2,'MarkerSize',14);
    plot(squeeze(SES_i(3,1,:)),'o-','LineWidth',2,'MarkerSize',14);
    names_1 = strcat(names, ' Internal');
    names_2 = strcat(names, ' External');
    n_names = [names_1; names_2];
    legend(n_names(:));
    title(['Internal vs External Session ' p_title]);
    xticks([1:mixtures]);
    xticklabels(x_labels);
    set(gca,'FontSize',22);
    ylim([0 1]);
    
    % Subject Based (mean) CCR Performance
    figure('name',[folder_name ' | Subject CCR'],'NumberTitle','off');
    grid on;
    hold on;
    plot(SUB_m(:,1,1),'--');
    plot(squeeze(SUB_g(:,9,1)),'sq-');
    plot(squeeze(SUB_i(:,9,1)),'o-');
    legend(names);
    title(p_title);
    set(gca,'FontSize',22);
    ylim([0 1]);
else
    mahal_data = {SUB_SES_m, SUB_m, SES_m};
    gmm_data = {SUB_SES_g, SUB_g, SES_g};
    ivector_data = {SUB_SES_i, SUB_i, SES_i};
end
end

function [sub_ses, sub, ses] = dataGet(folder)

sub_ses = iVectorBinary([folder filesep 'sub_ses.bin']);
sub = iVectorBinary([folder filesep 'sub.bin']);
ses = iVectorBinary([folder filesep 'ses.bin']);

end

function this_plot(x,y,line_style,l_size,m_size)
plot(x,y,line_style,'LineWidth',l_size,'MarkerSize',m_size);
end