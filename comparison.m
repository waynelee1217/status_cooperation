%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code description: This code is used to plot the difference of Mark's
% analytics and our analytics given a certain simulation setting.  

% Author: Hsuan-Wei Wayne Lee
% Contact information: hwwaynelee@sinica.edu.tw
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc; 
close all;

n = 10;
c_prob = 1/n;
trial = 100000;

fc_n = 9;
fc_results = zeros(11,2);


kk = 1;
for fc_n = 0:n
    
    fc_n

    trial_results = zeros(n, 2); % number of leaders at the top & failure times

    for jj = 1:trial

    flag = 0;
    count_F = 0;

    while flag == 0

        node_matrix1 = zeros(n, 3);

        for ii = 1:n

            node_matrix1(ii,1) = ii; % node id

            if ismember(ii, 1:fc_n)         
                node_matrix1(ii,2) = 1; % being a cooperator
            end

            if node_matrix1(ii,2) == 1 && rand < c_prob
                node_matrix1(ii,3) = 1; % going to the top level 
            end

        end

        if sum(node_matrix1(:,3))>1
            count_F = count_F + 1;
        end

        if sum(node_matrix1(:,3))==0
            trial_results(jj,1) = 0;
            flag = 1;
        end

        if sum(node_matrix1(:,3))==1
            trial_results(jj,1) = 1;
            flag = 1;
        end

        trial_results(jj,2) = count_F;

    end


    end

    %trial_results


    fc_results(kk,1) = sum(trial_results(:,1))/trial; % p1

    fc_results(kk,2) = sum(trial_results(:,2))/trial;
    
    kk = kk + 1;

end

fc_results


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mark's estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fc_n = 0:n;

p0_mark = ((n-1)/n).^fc_n;

p1_mark = 1 - p0_mark;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ours estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p0_raw = ((n-1)/n).^fc_n;

p1_raw = fc_n.*(1/n).*((n-1)/n).^(fc_n-1);

p1_ours = p1_raw./(p0_raw+p1_raw);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% plot(fc_n,fc_results(:,1),'-o','LineWidth',2,'Color','k')
% hold on
% plot(fc_n,p1_mark,'-o','LineWidth',2,'Color','r')
% hold on
% plot(fc_n,p1_ours,'-o','LineWidth',2,'Color','b')
% hold on


diff_mark = p1_mark-fc_results(:,1)';
diff_ours = p1_ours-fc_results(:,1)';

% plot(fc_n,fc_results(:,1),'-kp',... 
%     'LineWidth',0.001,...
%     'MarkerSize',40,...
%     'MarkerEdgeColor','none',...
%     'MarkerFaceColor',[0.8 0.8 0.8])
% hold on


% sz = 40;
% scatter(fc_n,fc_results(:,1),sz,...
%               'MarkerEdgeColor',[1 1 1],...
%               'MarkerFaceColor',[0.8 0.8 0.8],...
%               'LineWidth',0.01,'p')

sz = 1000;
scatter(fc_n,fc_results(:,1),sz,'p','MarkerEdgeColor',[1 1 1],...
              'MarkerFaceColor',[0.8 0.8 0.8],...
              'LineWidth',0.1)

hold on

plot(fc_n,p1_mark,'-ro',... 
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','r')
hold on

plot(fc_n,p1_ours,'-b^',...   
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','b')
hold on




set(gca, 'XTick', [0 2 4 6 8 10],'FontSize',16)
set(gca, 'YTick', [0 0.25 0.5 0.75 1],'FontSize',16)

axis([-0.2 10.2 -0.01 1.01])

xlabel('Number of status cooperators','FontSize',20,'FontWeight','bold')
ylabel('Probability of status hierarchy','FontSize',20,'FontWeight','bold')

box on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots: failure times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot(fc_n,fc_results(:,2),'-ro',... 
%     'LineWidth',2,...
%     'MarkerSize',10,...
%     'MarkerEdgeColor','k',...
%     'MarkerFaceColor','r')
% hold on
% 
% set(gca, 'XTick', [0 2 4 6 8 10],'FontSize',16)
% set(gca, 'YTick', [0 0.1 0.2 0.3 0.4],'FontSize',16)
% 
% axis([-0.2 10.2 -0.01 0.41])
% 
% xlabel('Number of status cooperators','FontSize',20,'FontWeight','bold')
% ylabel('Percentage of failed attempt','FontSize',20,'FontWeight','bold')




