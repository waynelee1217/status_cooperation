%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code description: This code is used to plot the status cooperation with
% hierarchy pursing with assortative mixing (\tau). Here we plot the
% analytical approximation range that cooperation is a stable strategy 
% in social dilemmas with different \tau. The sizes of groups are
% set to be from 2 to 10. 

% Author: Hsuan-Wei Wayne Lee
% Contact information: hwwaynelee@sinica.edu.tw
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameter setting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 10;

L = 2:N;


mark_upper = (L.^(L-1))./(L.^L-(L-1).^L);
lower_boud = 1./L;

plot(L,mark_upper,'LineWidth',1.5,'Color','k')
hold on


plot(L,lower_boud,'LineWidth',1.5,'Color','b')
hold on



n = 5; % you probably want to fix n
N = 10;

tau = 0;
b = 1;
c = 1;

color = jet(12);

BBB = 1;

for tau = 0:0.25:1
    
    upper_bound = [];

for n = 2:N

    fc = 1;
    
    % average payoff for C
   

    first_term = 0;
        
    sum_i = 0;

    for ii = 0:(n-1)
            
            sum_j = 0;
            
            for jj = 0:(ii+1)
                
                sum_j = sum_j + nchoosek((ii+1),jj)*(1/n)^jj*((n-1)/n)^(ii+1-jj)*(1-H_n(n,jj));

            end
            
            sum_i = sum_i + (((ii+1)/n)*(1/(ii+1)*nchoosek(n-1,ii)*(tau+(1-tau)*(fc+eps))^(ii)*((1-tau)*(1-(fc+eps)))^(n-1-ii)+...
                           ii/(ii+1)*(fc+eps)*nchoosek(n-1,ii)*(tau+(1-tau)*(fc+eps))^(ii-1)*((1-tau)*(1-(fc+eps)))^(n-1-ii))+...
                           ((n-ii-1)/n)*((1-(fc+eps))*nchoosek(n-1,ii)*((1-tau)*(fc+eps))^(ii)*(tau+(1-tau)*(1-(fc+eps)))^(n-2-ii)))*sum_j;   %(n-2-ii)??
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
            
            sum_ii = sum_ii + (((ii+1)/n)*(1/(ii+1)*nchoosek(n-1,ii)*(tau+(1-tau)*(fc+eps))^(ii)*((1-tau)*(1-(fc+eps)))^(n-1-ii)+...
                            ii/(ii+1)*(fc+eps)*nchoosek(n-1,ii)*(tau+(1-tau)*(fc+eps))^(ii-1)*((1-tau)*(1-(fc+eps)))^(n-1-ii))+...
                           ((n-ii-1)/n)*((1-(fc+eps))*nchoosek(n-1,ii)*((1-tau)*(fc+eps))^(ii)*(tau+(1-tau)*(1-(fc+eps)))^(n-2-ii)))*(sum_j);  %(n-2-ii)??
 
    end
        
    second_term = sum_ii;


    W_C = first_term + second_term;
    

    Ac = sum_i;
    Ab = second_term;
    

    


    
    % average payoff for D
    
    
    
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
            
            sum_i = sum_i + ((ii/n)*((fc+eps)*nchoosek(n-1,ii)*(tau+(1-tau)*(fc+eps))^(ii-1)*((1-tau)*(1-(fc+eps)))^(n-ii-1))+...
                    (n-ii)/n*((1/(n-ii)*nchoosek(n-1,ii)*((1-tau)*(fc+eps))^(ii)*(tau+(1-tau)*(1-(fc+eps)))^(n-ii-1))+((n-ii-1)/(n-ii)*(1-(fc+eps))*nchoosek(n-1,ii)*((1-tau)*(fc+eps))^(ii)*(tau+(1-tau)*(1-(fc+eps)))^(n-ii-2))))*(sum_j);
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


plot(L,upper_bound,'b--o','LineWidth',3,'Color',color(12-BBB,:))
hold on



axis([1.95 10.05 -0.05 0.85])
box on

xlabel('n','FontSize',20,'FontWeight','bold')
ylabel('c/b','FontSize',20,'FontWeight','bold')
%set(gca, 'XTick', [0 1 2 3 4 5 6 7 8 9 10],'FontSize',16)
set(gca, 'XTick', [1 2 3 4 5 6 7 8 9 10],'FontSize',16)

BBB = BBB + 1;

end


% plot(L,mark_upper,'LineWidth',3,'Color','k')
% hold on
% 
% 
% plot(L,lower_boud,'LineWidth',3,'Color','b')
% hold on




title('Cooperation is a stable strategy in social dilemmas with different \tau','FontSize',20,'FontWeight','bold')









