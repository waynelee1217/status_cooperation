%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code description: This code is used to plot the status cooperation with
% hierarchy pursing with assortative mixing (\tau). Here we plot the
% averaged simulation results over 100,000 trials. The sizes of groups are
% set to be 10. 

% Author: Hsuan-Wei Wayne Lee
% Contact information: hwwaynelee@sinica.edu.tw
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear;
% clc; 
% close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameter setting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hold on

n = 100000;

trial = 10;
b = 1;

WC_WD_all = zeros(2501,trial);
colors = jet(5);

fc = (n-1)/n; 
c_prob = 1/n;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mm = 1;

for tau = 0:0.25:1
    
    tau

for kk = 1:trial
  
    
    WC_WD_vec = [];

    for c = 0:0.01:0.60

        for fc = 0:0.025:1

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % WC
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

node_matrix1 = zeros(n, 4);

cc1 = rand;
cc2 = rand;

front_prob = 1/n;
front1 = rand;
front2 = rand;

            for ii = 1:n

                node_matrix1(ii,1) = ii; % node id
                

%                 if rand < fc
%                    node_matrix1(ii,2) = 1; % being a cooperator
%                 end
%                 
%                 if ii == 1
%                    node_matrix1(ii,2) = 1; % being a cooperator
%                 end

                
                if cc1 < fc && ii == 1
                    node_matrix1(ii,2) = 1; % being a cooperator
                end
                
                if front1 < front_prob && ii == 1
                    node_matrix1(ii,2) = 1; % being a cooperator
                end


                aa1 = rand;
                bb1 = rand;

                if node_matrix1(1,2) == 1 && (ii ~= 1)
                    if aa1 < tau 
                       node_matrix1(ii,2) = 1; % being a cooperator 

                    end

                    if (aa1 > tau) && (bb1 < fc)
                       node_matrix1(ii,2) = 1; % being a cooperator

                    end
                end
                
                if node_matrix1(1,2) ~= 1 && (ii ~= 1)
                    if aa1 < tau 
                       node_matrix1(ii,2) = 0; % being a defector 
  
                    end

                    if (aa1 > tau) && (bb1 < fc)
                       node_matrix1(ii,2) = 1; % being a cooperator

                    end
                end
                
    
                
                if front1 > front_prob && ii == n % last person is the protagonist 
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
            
            WC = c*(1-H_value1)+sum(node_matrix1(:,4))*b/n;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % WD
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            node_matrix2 = zeros(n, 4);

            for ii = 1:n

                node_matrix2(ii,1) = ii; % node id
               
%                 if rand < fc
%                    node_matrix2(ii,2) = 1; % being a cooperator               
%                 end
%                 
%                 if  ii == 1 
%                    node_matrix2(ii,2) = 0; % being a D
%                 end



                
                if cc2 < fc && ii == 1
                    node_matrix2(ii,2) = 1; % being a cooperator
                end
                
                if front2 < front_prob && ii == 1
                    node_matrix2(ii,2) = 0; % being a D
                end


                aa2 = rand;
                bb2 = rand;

                if node_matrix2(1,2) == 1 && (ii ~= 1)
                    if aa2 < tau 
                       node_matrix2(ii,2) = 1; % being a cooperator 
                    end

                    if (aa2 > tau) && (bb2 < fc)
                       node_matrix2(ii,2) = 1; % being a cooperator
                    end
                end
                
                if node_matrix2(1,2) ~= 1 && (ii ~= 1)
                    if aa2 < tau 
                       node_matrix2(ii,2) = 0; % being a defector 
                    end

                    if (aa2 > tau) && (bb2 < fc)
                       node_matrix2(ii,2) = 1; % being a cooperator
                    end
                end
                

                if front2 > front_prob && ii == n % last person is the protagonist 
                   node_matrix2(ii,2) = 0; % being a D
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
WC_WD_matrix = reshape(WC_WD_vec_ave,[41,61])';


value = [];
for ii = 1:41
    [M,I] = min(abs(WC_WD_matrix(:,ii)));  
    value_new = I;   
    value = [value value_new];  
end

%axis([0 42 0 71])
xx = 1:5:201;
yyy = value;
%plot(xx, yyy,'LineWidth',2,'Color',colors(mm, :))

sz = 60;
scatter(xx,yyy,sz,'MarkerEdgeColor',[0 0 0],...
              'MarkerFaceColor',colors(mm, :),...
              'LineWidth',0.3)


hold on

mm = mm + 1;

end

set(gca,'Xtick',1:40:201,'XTickLabel',{'0','0.2','0.4','0.6','0.8','1.0'},'FontSize',16)
set(gca,'Ytick',[1 11 21 31 41 51 61],'YTickLabel',{'0','0.1','0.2','0.3','0.4','0.5','0.6'},'FontSize',16) 

axis([1 201 0 62])

box on
hold on

xlabel('fc','FontSize',16,'FontWeight','bold')
ylabel('c/b','FontSize',16,'FontWeight','bold')
%title('MATLAB simulation (n = 10)','FontSize',16,'FontWeight','bold')

