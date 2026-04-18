clear all

data=xlsread('MONETARY DATASET.xlsx','ALL DATA','B2:K223');
y = [diff(log(data(1:end,[1 8 5])))*100 data(2:end,4)];
z = data(2:end,3);
var_names = {'PRICE LEVEL (cum. inflation)', 'RESERVES (log level, cum.)', 'REER (log level, cum.)', 'MMR'};
max_p = 13;
AIC = zeros(1,max_p);
for p = 1:max_p
    [beta,res,T,N,X] = OLS(y,p,z);
    T_eff = size(res,1);
    Sigma = res.'*res/T_eff;
    AIC(p) = log(det(Sigma)) + 2*(size(X,2)*N)/T_eff;  
end
[best_AIC, optimal_p] = min(AIC);
fprintf('Optimal lag for VAR: p = %d (AIC = %.4f)\n', optimal_p, best_AIC);
 
[A_best, res_best, T,N,X_best] = OLS(y,optimal_p,z);
T_eff = size(y,1)-optimal_p;
C = companion(A_best,optimal_p,N);
H = 40;
J = [eye(N), zeros(N, N*(optimal_p-1))];
C_h = zeros(N, N, H+1);
C_h(:,:,1) = eye(N);
for h = 1:H
    Fh = C^h;                       
    C_h(:,:,h+1) = J * Fh * J.';
end
IRF = zeros(N, H+1);
z_ind = 1+ (N*optimal_p)+1;
B = A_best(z_ind, :)'; 

for h = 0:H
    IRF(:,h+1) = C_h(:,:,h+1) * B;
end
% --- CUMULATE IRF for variables in differences (1:3) ---
cum_idx = 1:3;                 % INFLATION, RESERVES, REER
IRF(cum_idx,:) = cumsum(IRF(cum_idx,:), 2);


Boot = 1000;
cl= 0.68;
IRF_boot = zeros(N, H+1, Boot);
%make boot starting
for b = 1:Boot
    YB(1:optimal_p,:) = y(1:optimal_p,:);
    for t = optimal_p+1:T_eff+optimal_p
        resind = randi(T_eff,1,1);
        boot_res = res_best(resind,:);
        X_boot = [];
        for lag = 1:optimal_p
            X_boot = [X_boot YB(t-lag, :)];   
        end
        X_boot = [1 X_boot z(t, :)];  
        YB(t,:)= X_boot *A_best + boot_res;
    end
    [A_boot] = OLS(YB,optimal_p,z);
    C_boot = companion(A_boot, optimal_p, N);
    B_boot = A_boot(z_ind,:).';
    IRF_boot(:,1,b) = B_boot;  
    for h = 1:H
        Psi_h = J*(C_boot^h)*J.';
        IRF_boot(:,h+1,b) = Psi_h*B_boot;  
    end
end

% --- CUMULATE BOOTSTRAP IRFs for variables in differences (1:3) ---
cum_idx = 1:3;
IRF_boot(cum_idx,:,:) = cumsum(IRF_boot(cum_idx,:,:), 2);

% irf=cumsum(irf([1 2 3],:),2);
% irfboot=cumsum(IRF_boot([1 2 3],:,:),2)

IRF_lower = prctile(IRF_boot,16,3);
IRF_upper = prctile(IRF_boot,84,3);


figure('Color', 'w');
t_axis = 0:H; 

for i = 1:N
    subplot(2, 2, i);
    hold on;
    
    boot_series = squeeze(IRF_boot(i, :, :)); 
    low = squeeze(IRF_lower(i, :))';           
    up  = squeeze(IRF_upper(i, :))';           
    
    is_inside = all(boot_series >= low & boot_series <= up, 1);
    filtered_boot = boot_series(:, is_inside);
    
    
    if ~isempty(filtered_boot)
        plot(t_axis, filtered_boot, 'Color', [0.7 0.7 0.7 0.1], 'LineWidth', 0.3);
    end

    

    t_axis_row = t_axis(:)'; 
    low_row = low(:)';         
    up_row = up(:)';          
    
    
    x_fill = [t_axis_row, fliplr(t_axis_row)];
    y_fill = [low_row, fliplr(up_row)];
    
    
   
    if i == 1
        fill_color = [0.2 0.6 0.8]; 
    elseif i == 2
        fill_color = [0.8 0.4 0.2]; 
    elseif i == 3
        fill_color = [0.3 0.7 0.4]; 
    else
        fill_color = [0.5 0.5 0.5]; 
    end
    
    
    fill(x_fill, y_fill, fill_color, 'FaceAlpha', 0.25, 'EdgeColor', 'none');
    
    
    plot(t_axis_row, IRF(i, :), 'Color', fill_color * 0.7, 'LineWidth', 2.5);
    
    title(var_names{i}, 'FontSize', 11);
    xlabel('Horizon', 'FontSize', 9);
    ylabel('Response', 'FontSize', 9);
    grid on;
    line([0 H], [0 0], 'Color', 'k', 'LineWidth', 0.5);
    
    hold off;
end
sgtitle('IRF Oil Kanzig Shock - VARX (2007-2025) 68% cf', 'FontSize', 14, 'FontWeight', 'bold');