%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code description: Here we plot the analytical approximation range that 
% cooperation is a stable strategy in social dilemmas with different \tau.
% The sizes of groups are set to be from 3 to 10. 

% Author: Hsuan-Wei Wayne Lee
% Contact information: hwwaynelee@sinica.edu.tw
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






clear;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameter setting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

color = jet(5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n = 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('n = 3')
left= 0.2;
bottom=0.78;
width=0.18;
height=0.18;

axes('Position',[left bottom width height])

tau = 0;

color = jet(5);
b = 1; 
c = 1;

mm = 1;

for n = 3:3
    
    for tau = 0:0.25:1

        % average payoff for D

        diff_WC_WD = zeros(71,201);

        count_col = 1;

        for fc = 0:0.005:1

            count_row = 1;

            for c = 0:0.01:0.7

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
        plot(xx, yyy_n,'LineWidth',2,'Color',color(mm,:))
        
        mm = mm + 1;
        
    end


end

set(gca,'Xtick',1:40:201,'XTickLabel',{'0','0.2','0.4','0.6','0.8','1.0'},'FontSize',16) 
set(gca,'Ytick',[1 11 21 31 41 51 61 71],'YTickLabel',{'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7'},'FontSize',16) 

axis([1 201 0 72])

%xlabel('fc','FontSize',16,'FontWeight','bold')
ylabel('c/b','FontSize',16,'FontWeight','bold')
title('n = 3 ','FontSize',16,'FontWeight','bold')
box on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n = 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('n = 4')
left= 0.46;
bottom=0.78;
width=0.18;
height=0.18;

axes('Position',[left bottom width height])


tau = 0;

color = jet(5);

mm = 1;

for n = 4:4
    
    for tau = 0:0.25:1

        % average payoff for D

        diff_WC_WD = zeros(71,201);

        count_col = 1;

        for fc = 0:0.005:1

            count_row = 1;

            for c = 0:0.01:0.7

                % average payoff for D

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
        plot(xx, yyy_n,'LineWidth',2,'Color',color(mm,:))
        
        mm = mm + 1;
        
    end


end

set(gca,'Xtick',1:40:201,'XTickLabel',{'0','0.2','0.4','0.6','0.8','1.0'},'FontSize',16) 
set(gca,'Ytick',[1 11 21 31 41 51 61 71],'YTickLabel',{'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7'},'FontSize',16) 

axis([1 201 0 72])

%xlabel('fc','FontSize',16,'FontWeight','bold')
%ylabel('c/b','FontSize',16,'FontWeight','bold')
title('n = 4 ','FontSize',16,'FontWeight','bold')
box on



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n = 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('n = 5')
left= 0.2;
bottom=0.54;
width=0.18;
height=0.18;

axes('Position',[left bottom width height])


tau = 0;

color = jet(5);

mm = 1;

for n = 5:5
    
    for tau = 0:0.25:1

        % average payoff for D

        diff_WC_WD = zeros(71,201);

        count_col = 1;

        for fc = 0:0.005:1

            count_row = 1;

            for c = 0:0.01:0.7

                % average payoff for D

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
        plot(xx, yyy_n,'LineWidth',2,'Color',color(mm,:))
        
        mm = mm + 1;
        
    end


end

set(gca,'Xtick',1:40:201,'XTickLabel',{'0','0.2','0.4','0.6','0.8','1.0'},'FontSize',16) 
set(gca,'Ytick',[1 11 21 31 41 51 61 71],'YTickLabel',{'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7'},'FontSize',16) 

axis([1 201 0 72])

%xlabel('fc','FontSize',16,'FontWeight','bold')
ylabel('c/b','FontSize',16,'FontWeight','bold')
title('n = 5 ','FontSize',16,'FontWeight','bold')
box on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n = 6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('n = 6')
left= 0.46;
bottom=0.54;
width=0.18;
height=0.18;

axes('Position',[left bottom width height])


tau = 0;

color = jet(5);

mm = 1;

for n = 6:6
    
    for tau = 0:0.25:1

        % average payoff for D

        diff_WC_WD = zeros(71,201);

        count_col = 1;

        for fc = 0:0.005:1

            count_row = 1;

            for c = 0:0.01:0.7

                % average payoff for D

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
        plot(xx, yyy_n,'LineWidth',2,'Color',color(mm,:))
        
        mm = mm + 1;
        
    end


end

set(gca,'Xtick',1:40:201,'XTickLabel',{'0','0.2','0.4','0.6','0.8','1.0'},'FontSize',16) 
set(gca,'Ytick',[1 11 21 31 41 51 61 71],'YTickLabel',{'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7'},'FontSize',16) 

axis([1 201 0 72])

%xlabel('fc','FontSize',16,'FontWeight','bold')
%ylabel('c/b','FontSize',16,'FontWeight','bold')
title('n = 6 ','FontSize',16,'FontWeight','bold')
box on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n = 7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('n = 7')
left= 0.2;
bottom=0.3;
width=0.18;
height=0.18;


axes('Position',[left bottom width height])


tau = 0;

color = jet(5);

mm = 1;

for n = 7:7
    
    for tau = 0:0.25:1

        % average payoff for D

        diff_WC_WD = zeros(71,201);

        count_col = 1;

        for fc = 0:0.005:1

            count_row = 1;

            for c = 0:0.01:0.7

                % average payoff for D

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
        plot(xx, yyy_n,'LineWidth',2,'Color',color(mm,:))
        
        mm = mm + 1;
        
    end


end

set(gca,'Xtick',1:40:201,'XTickLabel',{'0','0.2','0.4','0.6','0.8','1.0'},'FontSize',16) 
set(gca,'Ytick',[1 11 21 31 41 51 61 71],'YTickLabel',{'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7'},'FontSize',16) 

axis([1 201 0 72])

%xlabel('fc','FontSize',16,'FontWeight','bold')
ylabel('c/b','FontSize',16,'FontWeight','bold')
title('n = 7 ','FontSize',16,'FontWeight','bold')
box on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n = 8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('n = 8')
left= 0.46;
bottom=0.3;
width=0.18;
height=0.18;

axes('Position',[left bottom width height])


tau = 0;

color = jet(5);

mm = 1;

for n = 8:8
    
    for tau = 0:0.25:1

        % average payoff for D

        diff_WC_WD = zeros(71,201);

        count_col = 1;

        for fc = 0:0.005:1

            count_row = 1;

            for c = 0:0.01:0.7

                % average payoff for D

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
        plot(xx, yyy_n,'LineWidth',2,'Color',color(mm,:))
        
        mm = mm + 1;
        
    end


end

set(gca,'Xtick',1:40:201,'XTickLabel',{'0','0.2','0.4','0.6','0.8','1.0'},'FontSize',16) 
set(gca,'Ytick',[1 11 21 31 41 51 61 71],'YTickLabel',{'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7'},'FontSize',16) 

axis([1 201 0 72])

%xlabel('fc','FontSize',16,'FontWeight','bold')
%ylabel('c/b','FontSize',16,'FontWeight','bold')
title('n = 8 ','FontSize',16,'FontWeight','bold')
box on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n = 9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('n = 9')
left= 0.2;
bottom=0.06;
width=0.18;
height=0.18;
axes('Position',[left bottom width height])


tau = 0;

color = jet(5);

mm = 1;

for n = 9:9
    
    for tau = 0:0.25:1

        % average payoff for D

        diff_WC_WD = zeros(71,201);

        count_col = 1;

        for fc = 0:0.005:1

            count_row = 1;

            for c = 0:0.01:0.7

                % average payoff for D

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
        plot(xx, yyy_n,'LineWidth',2,'Color',color(mm,:))
        
        mm = mm + 1;
        
    end


end

set(gca,'Xtick',1:40:201,'XTickLabel',{'0','0.2','0.4','0.6','0.8','1.0'},'FontSize',16) 
set(gca,'Ytick',[1 11 21 31 41 51 61 71],'YTickLabel',{'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7'},'FontSize',16) 

axis([1 201 0 72])

xlabel('fc','FontSize',16,'FontWeight','bold')
ylabel('c/b','FontSize',16,'FontWeight','bold')
title('n = 9 ','FontSize',16,'FontWeight','bold')
box on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n = 10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('n = 10')
left= 0.46;
bottom=0.06;
width=0.18;
height=0.18;

axes('Position',[left bottom width height])

tau = 0;

color = jet(5);

mm = 1;

for n = 10:10
    
    for tau = 0:0.25:1

        % average payoff for D

        diff_WC_WD = zeros(71,201);

        count_col = 1;

        for fc = 0:0.005:1

            count_row = 1;

            for c = 0:0.01:0.7

                % average payoff for D

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
        plot(xx, yyy_n,'LineWidth',2,'Color',color(mm,:))
        
        mm = mm + 1;
        
    end


end

set(gca,'Xtick',1:40:201,'XTickLabel',{'0','0.2','0.4','0.6','0.8','1.0'},'FontSize',16) 
set(gca,'Ytick',[1 11 21 31 41 51 61 71],'YTickLabel',{'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7'},'FontSize',16) 

axis([1 201 0 72])

xlabel('fc','FontSize',16,'FontWeight','bold')
%ylabel('c/b','FontSize',16,'FontWeight','bold')
title('n = 10 ','FontSize',16,'FontWeight','bold')
box on




