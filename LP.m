clear all
data=xlsread('MONETARY DATASET.xlsx','ALL DATA','B2:K223');

%y = [diff(log(data(:,[1 8])))*100 diff(log(data(1:end,5)))*100 data(2:end,4)];
y = [diff(log(data(1:end,[1 8 5])))*100  data(2:end,4)];
z = [data(2:end,3)];
dummy18 = data(2:end,10);
ndummy18 = 1-dummy18;

var_names = {'PRICE LEVEL (cum. inflation)', 'RESERVES (log level, cum.)', 'REER (log level, cum.)', 'MMR'};


max_p = 12;

AIC = zeros(1,max_p);

for p = 1:max_p
    [beta,res,T,N] = OLS(y,p,0,0);
    Sigma = res.'*res/(T-p);
    AIC(p) = log(det(Sigma)) + (2*N^2*p)/T;  
end


[best_AIC, optimal_p] = min(AIC);
fprintf('Optimal lag for VAR: p = %d (AIC = %.4f)\n', optimal_p, best_AIC);


[A_best, res_best] = OLS(y,optimal_p,0,0);


T = size(y,1);
N = size(y,2);



s=40;

Theta = zeros(s+1, N);   
Beta_0 = zeros(s+1, N);
Beta_1 = zeros(s+1, N);



for i = 0:s
    ylpi = i;
    xlpi = T-i;

    Y = y(optimal_p+1+i:T,:); 
    X = [];
    for j = 1:optimal_p
        X = [X y(optimal_p+1-j:T-j-i,1:end)];
    end

    X_B = [ones(size(Y,1),1) z(optimal_p+1:T-i,:) X];

    shocks = (X_B.'*X_B)\(X_B.'*Y);
    theta_i = shocks(2,:);   
    Theta(i+1,:) = theta_i;
   
    
end

Theta_level = Theta;
Theta_level(:,1:3) = cumsum(Theta(:,1:3),1); 


h = 0:s;


Boot = 1000;
s = 40;
T = size(y,1);
N = size(y,2);

IRF_boot = zeros(s+1, N, Boot);

[~, res_full] = OLS(y, optimal_p, 0, 0);
Teff = size(res_full, 1);

for b = 1:Boot
    
    
    v = (rand(Teff,1) > 0.5)*2 - 1;
    res_boot = res_full .* v;

  
    YB = y;
    YB(optimal_p+1:end,:) = YB(optimal_p+1:end,:) + res_boot;

    
    for h = 0:s
        
        Yh = YB(optimal_p+1+h:T, :);
        
        
        Xlags = [];
        for j = 1:optimal_p
            Xlags = [Xlags YB(optimal_p+1-j:T-j-h, :)];
        end
        
        
        X = [ones(size(Yh,1),1), z(optimal_p+1:T-h,:), Xlags];
        
       
        bh = (X.'*X) \ (X.'*Yh);
        
        
        IRF_boot(h+1, :, b) = bh(2, :);
    end
end


cum_idx = 1:3;  
IRF_boot(:, cum_idx, :) = cumsum(IRF_boot(:, cum_idx, :), 1);


IRF_lo = prctile(IRF_boot, 16, 3);  
IRF_hi = prctile(IRF_boot, 84, 3); 


t = (0:s)';

figure('Color', 'w');
for j = 1:N
    subplot(2, 2, j);
    hold on;
    
    
    if j == 1
        fill_color = [0.2 0.6 0.8];
    elseif j == 2
        fill_color = [0.8 0.4 0.2]; 
    elseif j == 3
        fill_color = [0.3 0.7 0.4]; 
    else
        fill_color = [0.5 0.5 0.5]; 
    end
    
    
    patch([t; flipud(t)], [IRF_lo(:,j); flipud(IRF_hi(:,j))], fill_color, 'FaceAlpha', 0.25,'EdgeColor', 'none');
    
    
    plot(t, Theta_level(:,j), 'Color', fill_color * 0.7, 'LineWidth', 2.5);
    
    yline(0, 'k--');
    grid on;
    title(var_names{j}, 'FontSize', 11);
    xlabel('Horizon', 'FontSize', 9);
    ylabel('Response', 'FontSize', 9);
    
   
    if j == 1
        legend('68% CI', 'IRF', 'Location', 'best');
    end
    
    hold off;
end

sgtitle('IRF Oil Kanzig Shock - LP (2007-2025) 68% cf', 'FontSize', 14, 'FontWeight', 'bold');