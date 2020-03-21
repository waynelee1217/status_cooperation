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

n = 5; % you probably want to fix n
N = 10;

tau = 0;
b = 1;
c = 1;

color = jet(12);

BBB = 1;

for tau = 0:0.2:1
    
    upper_bound = [];

for n = 2:N

    fc = 1;
    
    % average payoff for C

    double_sum_C = 0;
        
    sum_i = 0;

    for ii = 0:(n-1)
            
            sum_j = 0;
            
            for jj = 2:(ii+1)
                
                sum_j = sum_j + nchoosek((ii+1),jj)*(1/n)^jj*((n-1)/n)^(ii+1-jj)*(1-H_n(n,jj));

            end
            
            sum_i = sum_i + (((ii+1)/n)*(1/(ii+1)*nchoosek(n-1,ii)*(tau+(1-tau)*fc)^(ii)*((1-tau)*(1-fc))^(n-1-ii)+...
                            ii/(ii+1)*fc*nchoosek(n-1,ii)*(tau+(1-tau)*fc)^(ii-1)*((1-tau)*(1-fc))^(n-1-ii))+...
                            ((n-ii-1)/n)*((1-fc)*nchoosek(n-1,ii)*((1-tau)*fc)^(ii)*(tau+(1-tau)*(1-fc))^(n-2-ii)))*(nchoosek((ii+1),0)*(1/n)^0*((n-1)/n)^(ii+1)+sum_j);
    end

    double_sum_C = sum_i;
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    triple_sum_C = 0;
    
    
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
            
            sum_ii = sum_ii + (((ii+1)/n)*(1/(ii+1)*nchoosek(n-1,ii)*(tau+(1-tau)*fc)^(ii)*((1-tau)*(1-fc))^(n-1-ii)+...
                                ii/(ii+1)*fc*nchoosek(n-1,ii)*(tau+(1-tau)*fc)^(ii-1)*((1-tau)*(1-fc))^(n-1-ii))+...
                                ((n-ii-1)/n)*((1-fc)*nchoosek(n-1,ii)*((1-tau)*fc)^(ii)*(tau+(1-tau)*(1-fc))^(n-2-ii)))*(nchoosek((ii+1),1)*(1/n)*((n-1)/n)^(ii)*((ii+1)/n)+sum_j);
 
    end
        
    triple_sum_C = sum_ii;
    

    W_C = c*(((1/n)*((tau+(1-tau)*fc)^(n-1))+((n-1)/n)*(fc*(tau+(1-tau)*fc)^(n-2)))*(1/n)^n*((n-1)/(n))^0 + double_sum_C) + triple_sum_C;

    Ac = ((1/n)*((tau+(1-tau)*fc)^(n-1))+((n-1)/n)*(fc*(tau+(1-tau)*fc)^(n-2)))*(1/n)^n*((n-1)/(n))^0 + double_sum_C;
    Ab = triple_sum_C;

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
            
            sum_i = sum_i + ((ii/n)*(fc*nchoosek(n-1,ii)*(tau+(1-tau)*fc)^(ii-1)*((1-tau)*(1-fc))^(n-ii-1))+...
                   (n-ii)/n*((1/(n-ii)*nchoosek(n-1,ii)*((1-tau)*fc)^(ii)*(tau+(1-tau)*(1-fc))^(n-ii-1))+((n-ii-1)/(n-ii)*(1-fc)*nchoosek(n-1,ii)*((1-tau)*fc)^(ii)*(tau+(1-tau)*(1-fc))^(n-ii-2))))*...
                   (nchoosek(ii,1)*(1/n)*((n-1)/n)^(ii-1)*(ii/n)+sum_j);
               
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



axis([1.95 10.05 -0.05 2.05])
box on

xlabel('n','FontSize',20,'FontWeight','bold')
ylabel('c/b','FontSize',20,'FontWeight','bold')
%set(gca, 'XTick', [0 1 2 3 4 5 6 7 8 9 10],'FontSize',16)
set(gca, 'XTick', [1 2 3 4 5 6 7 8 9 10],'FontSize',16)

BBB = BBB + 1;

end


plot(L,mark_upper,'b--o','LineWidth',3,'Color','k')
hold on


plot(L,lower_boud,'b--o','LineWidth',3,'Color','b')
hold on




title('Cooperation is a stable strategy in social dilemmas with different \tau','FontSize',20,'FontWeight','bold')









