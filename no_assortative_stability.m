%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code description: This code is used to plot the status cooperation with
% hierarchy pursing with no assortative mixing. Here we plot the
% analytical approximation range that cooperation is a stable strategy 
% in social dilemmas with different \tau. The sizes of groups are
% set to be from 2 to 10. 

% Author: Hsuan-Wei Wayne Lee
% Contact information: hwwaynelee@sinica.edu.tw
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameter setting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 5; % you probably want to fix n
N = 10;

b = 2;
c = 1;

upper_bound = [];

for n = 2:N

    fc = 1;

    first_term = 0;
        
    sum_i = 0;

    for ii = 0:(n-1)
            
            sum_j = 0;
            
            for jj = 2:(ii+1)
                
                sum_j = sum_j + nchoosek((ii+1),jj)*(1/n)^jj*((n-1)/n)^(ii+1-jj)*(1-H_n(n,jj));

            end
            
            sum_i = sum_i + nchoosek(n-1,ii)*fc^ii*(1-fc)^(n-1-ii)*(nchoosek((ii+1),0)*(1/n)^0*((n-1)/n)^(ii+1)+sum_j);
    end

    first_term = (fc^(n-1)*(1/n)^n + sum_i)*c;
    
    
    sum_ii = 0;

    for ii = 0:(n-1)
            
            sum_j = 0;
            
            for jj = 2:(ii+1)
                
               sum_k = 0;
                
               for kk = 0:(ii+1)
                   
                   sum_k = sum_k + nchoosek(ii+1,kk)*H_n(n,jj)^kk*(1-H_n(n,jj))^(ii+1-kk)*((kk)/n);

               end
               
               sum_j = sum_j + nchoosek((ii+1),jj)*(1/n)^jj*((n-1)/n)^(ii+1-jj)*sum_k;
               
            end
            
            sum_ii = sum_ii + nchoosek(n-1,ii)*fc^ii*(1-fc)^(n-1-ii)*(nchoosek((ii+1),1)*(1/n)*((n-1)/n)^(ii)*((ii+1)/n)+sum_j);
 
    end
        
    second_term = sum_ii;


    W_C = first_term + second_term;
    

    Ac = fc^(n-1)*(1/n)^n + sum_i;
    Ab = second_term;


    fc = (n-1)/n;
    fc;

    % average payoff for D
    
    triple_sum_D = 0;

    sum_i = 0;

        for ii = 1:(n-1)
            
            sum_j = 0;
            
            for jj = 2:ii
                
               sum_k = 0;
               
               for kk = 1:ii
                   
                   sum_k = sum_k + nchoosek(ii,kk)*H_n(n,jj)^kk*(1-H_n(n,jj))^(ii-kk)*(kk/n);
                   
               end
               
               sum_j = sum_j + nchoosek(ii,jj)*(1/n)^jj*((n-1)/n)^(ii-jj)*sum_k;
               
            end
            
            sum_i = sum_i + nchoosek(n-1,ii)*fc^ii*(1-fc)^(n-1-ii)*(ii*(1/n)*((n-1)/n)^(ii-1)*(ii/n)+sum_j);
        end
        
    triple_sum_D = sum_i;

    W_D = c + triple_sum_D;

    Bc = 1;
    Bb = triple_sum_D;


    upper_bound = [upper_bound (Bb - Ab)/(Ac - 1)];
    
end

L = 2:N;


mark_upper = (L.^(L-1))./(L.^L-(L-1).^L);
lower_boud = 1./L;


plot(L,upper_bound,'b--o','LineWidth',3,'Color','r')
hold on


plot(L,mark_upper,'b--o','LineWidth',3,'Color','k')
hold on

plot(L,lower_boud,'b--o','LineWidth',3,'Color','b')
hold on

axis([1.95 10.05 -0.05 1.55])
box on

xlabel('n','FontSize',20,'FontWeight','bold')
ylabel('c/b','FontSize',20,'FontWeight','bold')
%set(gca, 'XTick', [0 1 2 3 4 5 6 7 8 9 10],'FontSize',16)
set(gca, 'XTick', [1 2 3 4 5 6 7 8 9 10],'FontSize',16)
%set(gca, 'XTick', [10 20 30 40 50],'FontSize',16)

title('Cooperation is a stable strategy in social dilemmas','FontSize',20,'FontWeight','bold')









