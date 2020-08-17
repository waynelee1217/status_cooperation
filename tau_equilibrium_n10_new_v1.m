clear;
clc;
close all;

hold on


b = 1; 
c = 1;

tau = 0;

color = jet(5);

mm = 1;

for n = 10:10
    
    for tau = 0:0.25:1
        
        tau


        % average payoff for D

        diff_WC_WD = zeros(61,201);

        count_col = 1;

        for fc = 0:0.005:1 %0:0.005:1

            count_row = 1;

            for c = 0:0.01:0.6

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

                        %sum_i = sum_i + nchoosek(n-1,ii)*(fc+eps)^ii*(1-(fc+eps))^(n-1-ii)*(sum_j);
                        sum_i = sum_i + ((ii/n)*((fc+eps)*nchoosek(n-1,ii)*(tau+(1-tau)*(fc+eps))^(ii-1)*((1-tau)*(1-(fc+eps)))^(n-ii-1))+...
                    (n-ii)/n*((1/(n-ii)*nchoosek(n-1,ii)*((1-tau)*(fc+eps))^(ii)*(tau+(1-tau)*(1-(fc+eps)))^(n-ii-1))+((n-ii-1)/(n-ii)*(1-(fc+eps))*nchoosek(n-1,ii)*((1-tau)*(fc+eps))^(ii)*(tau+(1-tau)*(1-(fc+eps)))^(n-ii-2))))*(sum_j);
                        
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

                        %sum_i = sum_i + nchoosek(n-1,ii)*(fc+eps)^ii*(1-(fc+eps))^(n-1-ii)*sum_j;
                           sum_i = sum_i + (((ii+1)/n)*(1/(ii+1)*nchoosek(n-1,ii)*(tau+(1-tau)*(fc+eps))^(ii)*((1-tau)*(1-(fc+eps)))^(n-1-ii)+...
                           ii/(ii+1)*(fc+eps)*nchoosek(n-1,ii)*(tau+(1-tau)*(fc+eps))^(ii-1)*((1-tau)*(1-(fc+eps)))^(n-1-ii))+...
                           ((n-ii-1)/n)*((1-(fc+eps))*nchoosek(n-1,ii)*((1-tau)*(fc+eps))^(ii)*(tau+(1-tau)*(1-(fc+eps)))^(n-2-ii)))*sum_j; 
                        
                        
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

                        %sum_ii = sum_ii + nchoosek(n-1,ii)*(fc+eps)^ii*(1-(fc+eps))^(n-1-ii)*(sum_j);
                           sum_ii = sum_ii + (((ii+1)/n)*(1/(ii+1)*nchoosek(n-1,ii)*(tau+(1-tau)*(fc+eps))^(ii)*((1-tau)*(1-(fc+eps)))^(n-1-ii)+...
                            ii/(ii+1)*(fc+eps)*nchoosek(n-1,ii)*(tau+(1-tau)*(fc+eps))^(ii-1)*((1-tau)*(1-(fc+eps)))^(n-1-ii))+...
                           ((n-ii-1)/n)*((1-(fc+eps))*nchoosek(n-1,ii)*((1-tau)*(fc+eps))^(ii)*(tau+(1-tau)*(1-(fc+eps)))^(n-2-ii)))*(sum_j);

                end

                second_term = sum_ii;

                W_C = first_term + second_term ;

                W_diff = W_C - W_D;

                diff_WC_WD(count_row,count_col) = W_diff;

                count_row = count_row + 1;

            end

            count_col = count_col + 1;

        end



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
        plot(xx, yyy_n,'LineWidth',3,'Color',color(mm, :))
        
        mm = mm + 1;
        
    end


end

set(gca,'Xtick',1:40:201,'XTickLabel',{'0','0.2','0.4','0.6','0.8','1.0'},'FontSize',16) 
set(gca,'Ytick',[1 11 21 31 41 51 61],'YTickLabel',{'0','0.1','0.2','0.3','0.4','0.5','0.6'},'FontSize',16) 

axis([1 201 0 62])

xlabel('(fc+eps)','FontSize',16,'FontWeight','bold')
ylabel('c/b','FontSize',16,'FontWeight','bold')
title({'Comparison of simulation and analytics','for different \tau values (n = 10)'},'FontSize',16,'FontWeight','bold')

box on
hold on
