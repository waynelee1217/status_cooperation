%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Date 12.22.2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;
close all;

hold on

n = 5; % you probably want to fix n

fc = 0.8;
b = 1; % fix b = 1 here, since we only care about c/b
c = 1;



color = jet(4);

for n = 4:2:10
    
    n

% average payoff for D

diff_WC_WD = zeros(51,201);

count_col = 1;

for fc = 0:0.005:1
    
    count_row = 1;
    
    for c = 0:0.01:0.50
        
    triple_sum_D = 0;

    sum_i = 0;

        for ii = 0:(n-1)
            
            sum_j = 0;
            
            for jj = 0:ii
                
               sum_k = 0;
               
               for kk = 0:ii
                   
                   sum_k = sum_k + nchoosek(ii,kk)*H_n(n,jj)^kk*(1-H_n(n,jj))^(ii-kk)*(kk/n);
                   
               end
               
               sum_j = sum_j + nchoosek(ii,jj)*(1/n)^jj*((n-1)/n)^(ii-jj)*sum_k;
               
            end
            
            sum_i = sum_i + nchoosek(n-1,ii)*fc^ii*(1-fc)^(n-1-ii)*(sum_j);
        end
        
    triple_sum_D = sum_i;

    W_D = c + triple_sum_D;

        % average payoff for C
        
    first_term = 0;
        
    sum_i = 0;

    for ii = 0:(n-1)
            
            sum_j = 0;
            
            for jj = 0:(ii+1)
                
                sum_j = sum_j + nchoosek((ii+1),jj)*(1/n)^jj*((n-1)/n)^(ii+1-jj)*(1-H_n(n,jj));

            end
            
            sum_i = sum_i + nchoosek(n-1,ii)*fc^ii*(1-fc)^(n-1-ii)*sum_j;
    end

    first_term = (sum_i)*c;
    
    
    sum_ii = 0;

    for ii = 0:(n-1)
            
            sum_j = 0;
            
            for jj = 0:(ii+1)
                
               sum_k = 0;
                
               for kk = 0:(ii+1)
                   
                   sum_k = sum_k + nchoosek(ii+1,kk)*H_n(n,jj)^kk*(1-H_n(n,jj))^(ii+1-kk)*((kk)/n);

               end
               
               sum_j = sum_j + nchoosek((ii+1),jj)*(1/n)^jj*((n-1)/n)^(ii+1-jj)*sum_k;
               
            end
            
            sum_ii = sum_ii + nchoosek(n-1,ii)*fc^ii*(1-fc)^(n-1-ii)*(sum_j);
 
    end
        
    second_term = sum_ii;


    W_C = first_term + second_term;
        
        W_diff = W_C - W_D;
        
        diff_WC_WD(count_row,count_col) = W_diff;
        
        count_row = count_row + 1;
    
    end
    
    count_col = count_col + 1;
    
end


% imagesc(diff_WC_WD)
% set(gca,'YDir','normal')

diff_zero = zeros(1,0);

for ii = 1:201
  
  [M, I] = min(abs(diff_WC_WD(:,ii)));
  diff_zero = [diff_zero,I];
  
end


hold on
value_b = diff_zero;
xx = 1:201;
yy = value_b;
yyy_n = value_b;
plot(xx, yyy_n,'LineWidth',2,'Color',color((n-2)/2,:))
%plot(xx, yyy_n,'o','LineWidth',2,'Color','k')


end

set(gca,'Xtick',1:40:201,'XTickLabel',{'0','0.2','0.4','0.6','0.8','1.0'},'FontSize',16) 
set(gca,'Ytick',[1 11 21 31 41 51],'YTickLabel',{'0','0.1','0.2','0.3','0.4','0.5'},'FontSize',16) 

axis([1 201 0 52])

xlabel('fc','FontSize',16,'FontWeight','bold')
ylabel('c/b','FontSize',16,'FontWeight','bold')
title('MATLAB analysis','FontSize',16,'FontWeight','bold')
box on

% diff_WC_WD_matrix = reshape(diff_WC_WD,[21,26]);
% diff_WC_WD_matrix = diff_WC_WD_matrix';
% imagesc(diff_WC_WD_matrix)
% set(gca,'YDir','normal')


% set(gca,'Xtick',1:20:101,'XTickLabel',{'0','0.2','0.4','0.6','0.8','1.0'},'FontSize',16) 
% set(gca,'Ytick',[1 301 601 901 1201 1501],'YTickLabel',{'0','0.1','0.2','0.3','0.4','0.5'},'FontSize',16) 
% 
% 
% xlabel('','FontSize',16,'FontWeight','bold')
% ylabel('c/b','FontSize',16,'FontWeight','bold')
% title('n = 3','FontSize',16,'FontWeight','bold')
% hcb=colorbar
% title(hcb,'WC - WD')
% 
% hold on
% value = [];
% for mm = 1:101
%     [M,I] = min(abs(diff_WC_WD_matrix(:,mm)));  
%     value_new = 0.02 * I;   
%     value = [value value_new];  
% end
% 
% value_b = value./30;
% xx = 1:101;
% yy = value_b;
% yyy_n3 = value_b * 1500;
% plot(xx, yyy_n3,'LineWidth',2,'Color','r')