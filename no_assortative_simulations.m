%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code description: This code is used to plot the status cooperation with
% hierarchy pursing with no assortative mixing. Here we plot the
% averaged simulation results over 100,000 trials. The sizes of groups are
% set to be 4, 6, 8, and 10. 

% Author: Hsuan-Wei Wayne Lee
% Contact information: hwwaynelee@sinica.edu.tw
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameter setting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 10;          % group size (temporary)
fc = (n-1)/n;    % fraction of C (temporary), put this here is dangerous since we will change n later
tau = 0;         % tau
c_prob = 1/n;    % prob. for a C going to top level, put this here is dangerous since we will change n later

b = 1;
c = 0.1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trial = 100000;
WC_WD_all = zeros(2091,trial);
colors = jet(4);

for n = 4:2:10
    
    n
    
    fc = (n-1)/n; 
    c_prob = 1/n;

for kk = 1:trial
  
    
    WC_WD_vec = [];

    for c = 0:0.01:0.50

        for fc = 0:0.025:1

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % WC
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            node_matrix1 = zeros(n, 4);

            for ii = 1:n

                node_matrix1(ii,1) = ii; % node id

                if ii == 1
                   node_matrix1(ii,2) = 1; % being a cooperator
                end

                if rand < fc && ii > 1
                   node_matrix1(ii,2) = 1; % being a cooperator
                end

                if node_matrix1(ii,2) == 1 && rand < c_prob
                   node_matrix1(ii,3) = 1; % going to the top level 
                end

            end

            x = sum(node_matrix1(:,3)); % number of nodes on the top level

            H_value1 = H_n(n,x);

            for ii = 1:n

                if node_matrix1(ii,2) == 1 && rand < H_value1 
                   node_matrix1(ii,4) = 1; % contributing 
                end

            end
            
            count_C_not_contribute = 0;
            
            for ii = 1:n

                if node_matrix1(ii,2) == 1 && node_matrix1(ii,4) == 0 
                   count_C_not_contribute = count_C_not_contribute + 1;
                end

            end
            

            %WC = c*(count_C_not_contribute/n)+sum(node_matrix1(:,4))*b/n; % payoff
            WC = c*(1-H_value1)+sum(node_matrix1(:,4))*b/n;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % WD
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            node_matrix2 = zeros(n, 4);

            for ii = 1:n

                node_matrix2(ii,1) = ii; % node id
           
                if  ii == 1 
                   node_matrix2(ii,2) = 0; % being a D
                end

                if rand < fc && ii > 1
                   node_matrix2(ii,2) = 1; % being a cooperator               
                end

                if node_matrix2(ii,2) == 1 && rand < c_prob
                   node_matrix2(ii,3) = 1; % going to the top level 
                end

            end

            x = sum(node_matrix2(:,3)); % number of nodes on the top level

            H_value2 = H_n(n,x);

            for ii = 1:n

                if node_matrix2(ii,2) == 1 && rand < H_value2  
                   node_matrix2(ii,4) = 1; % contributing 
                end

            end

            WD = c + sum(node_matrix2(:,4))*b/n; % payoff

            WC_WD = WC - WD;

            WC_WD_vec = [WC_WD_vec WC_WD];
            

        end

    end
    
    WC_WD_all(:,kk) = WC_WD_vec;
    kk;
    
    
end

WC_WD_vec_ave = mean(WC_WD_all,2);
WC_WD_matrix = reshape(WC_WD_vec_ave,[41,51])';




set(gca,'Xtick',1:8:41,'XTickLabel',{'0','0.2','0.4','0.6','0.8','1.0'},'FontSize',16)
set(gca,'Ytick',[1 11 21 31 41 51],'YTickLabel',{'0','0.1','0.2','0.3','0.4','0.5'},'FontSize',16) 


% xlabel('fc','FontSize',16,'FontWeight','bold')
% ylabel('c/b','FontSize',16,'FontWeight','bold')
% title('n = 3','FontSize',16,'FontWeight','bold')
% hcb=colorbar
% title(hcb,'WC - WD')


value = [];
for ii = 1:41
    [M,I] = min(abs(WC_WD_matrix(:,ii)));  
    value_new = I;   
    value = [value value_new];  
end


xx = 1:5:201;
yyy = value;
%plot(xx, yyy,'LineWidth',2,'Color',colors((n-2)/2, :))

sz = 60;
scatter(xx,yyy,sz,'MarkerEdgeColor',[0 0 0],...
              'MarkerFaceColor',colors((n-2)/2, :),...
              'LineWidth',0.3)

hold on


end

set(gca,'Xtick',1:40:201,'XTickLabel',{'0','0.2','0.4','0.6','0.8','1.0'},'FontSize',16)
set(gca,'Ytick',[1 11 21 31 41 51],'YTickLabel',{'0','0.1','0.2','0.3','0.4','0.5'},'FontSize',16) 

axis([1 201 0 52])


xlabel('fc','FontSize',16,'FontWeight','bold')
ylabel('c/b','FontSize',16,'FontWeight','bold')
title('MATLAB simulation','FontSize',16,'FontWeight','bold')
title('Comparison of simulation and analytic for different group size n','FontSize',16,'FontWeight','bold')

