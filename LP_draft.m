clear all
data=xlsread('MONETARY DATASET.xlsx','ALL DATA','B2:K223');

%y = [diff(log(data(:,[1 8])))*100 diff(log(data(1:end,5)))*100 data(2:end,4)];
y = [log(data(2:end,[1 8 5]))*100 data(2:end,4)];
z = [data(2:end,3)];
dummy18 = data(2:end,10);
ndummy18 = 1-dummy18;

var_names = {'INFLATION', 'RESERVES', 'REER', 'MMR'};


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

    X_D = [ones(size(Y,1),1), dummy18(optimal_p+1:T-i,:), z(optimal_p+1:T-i,:) .* ndummy18(optimal_p+1:T-i,:), z(optimal_p+1:T-i,:) .* dummy18(optimal_p+1:T-i,:),  X, dummy18(optimal_p+1:T-i,:) .* X];
    shocks_D = (X_D.'*X_D) \ (X_D.'*Y);
    beta_0 = shocks_D(3,:); %pre D=0
    beta_1 = shocks_D(4,:); %post D=1
    Beta_0(i+1,:) = beta_0;
    Beta_1(i+1,:) = beta_1;
    
   
   

    
end





h = 0:s;

figure;
for j = 1:N
    subplot(2,2,j)
    plot(h, Theta(:,j), 'LineWidth', 1.8)
    yline(0,'k--')
    title(var_names{j})
    xlabel('Horizon')
    grid on
end

sgtitle('Local Projection IRF – Oil Shock Kanzig 2007-2025');

figure;
for j = 1:N
    subplot(2, 2, j);
    
    % Plot delle due IRF sullo stesso grafico
    plot(h, Beta_0(:,j), 'b-', 'LineWidth', 2, 'DisplayName', 'Pre 2018 (D=0)');
    hold on;
    plot(h, Beta_1(:,j), 'r-', 'LineWidth', 2, 'DisplayName', 'Post 2018 (D=1)');
    yline(0, 'k--');
    
    title(var_names{j});
    xlabel('Horizon');
    ylabel('Response');
    grid on;
    
    % Legenda solo nel primo subplot (o in tutti se vuoi)
    if j == 1
        legend('show', 'Location', 'best');
    end
    
    hold off;
end

sgtitle('Local Projection IRF – Shock z (by regime)');






% =========================
% 1) BOOTSTRAP LP (crea IRF0_boot e IRF1_boot)
% =========================
Boot = 1000;

IRF0_boot = zeros(s+1, N, Boot);   % (horizon x var x boot)
IRF1_boot = zeros(s+1, N, Boot);

% residui dal VAR(p) di base (una volta sola, non dentro il loop)
[~, res_full] = OLS(y, optimal_p, 0, 0);
Teff = size(res_full,1);

for b = 1:Boot
    
    % resample residui (i.i.d. bootstrap)
    v = (rand(Teff,1) > 0.5)*2 - 1;
    res_boot = res_full .* v;

    % genera YB (versione semplice: aggiungi residui ricampionati)
    YB = y;
    YB(optimal_p+1:end,:) = YB(optimal_p+1:end,:) + res_boot;

    % rifai LP regime-dependent su YB
    for h = 0:s
        
        Yh = YB(optimal_p+1+h:T,:);
        
        Xlags = [];
        for j = 1:optimal_p
            Xlags = [Xlags YB(optimal_p+1-j:T-j-h,:)];
        end
        
        X_D = [ones(size(Yh,1),1), ...
               dummy18(optimal_p+1:T-h,:), ...
               z(optimal_p+1:T-h,:) .* ndummy18(optimal_p+1:T-h,:), ...
               z(optimal_p+1:T-h,:) .* dummy18(optimal_p+1:T-h,:), ...
               Xlags, ...
               dummy18(optimal_p+1:T-h,:) .* Xlags];
        
        bh = (X_D.'*X_D)\(X_D.'*Yh);

        IRF0_boot(h+1,:,b) = bh(3,:);   % pre (z * (1-D))
        IRF1_boot(h+1,:,b) = bh(4,:);   % post (z * D)
    end
end

% =========================
% 2) DIFFERENZA: point estimate + bootstrap bands
% =========================
dIRF = Beta_1 - Beta_0;                 % (s+1) x N
dIRF_boot = IRF1_boot - IRF0_boot;      % (s+1) x N x Boot

d_lo68 = prctile(dIRF_boot,16,3);       % (s+1) x N
d_hi68 = prctile(dIRF_boot,84,3);

% =========================
% 3) PLOT differenza (corretto in dimensioni)
% =========================
t = (0:s)';

figure;
for j = 1:N
    subplot(N,1,j); hold on;

    lo = d_lo68(:,j);
    hi = d_hi68(:,j);

    patch([t; flipud(t)], [lo; flipud(hi)], ...
          [0.7 0.7 0.7], 'FaceAlpha',0.35, 'EdgeColor','none');

    plot(t, dIRF(:,j), 'k', 'LineWidth', 1.8);

    yline(0,'k--');
    grid on;
    title(['\Delta IRF (Post - Pre): ' var_names{j}]);
    xlabel('Horizon');
end

sgtitle('\Delta IRF (LP) with 68% bootstrap bands');


% =========================
% BOOTSTRAP BANDS for IRF0 and IRF1 (percentile bands)
% =========================

% Point estimates already computed:
% Beta_0  (s+1 x N)  pre
% Beta_1  (s+1 x N)  post
% Bootstrap draws:
% IRF0_boot (s+1 x N x Boot)
% IRF1_boot (s+1 x N x Boot)

% Choose confidence level (68% like your delta plot)
% lo_q = 16;
% hi_q = 84;
% 
% IRF0_lo = prctile(IRF0_boot, lo_q, 3);   % (s+1 x N)
% IRF0_hi = prctile(IRF0_boot, hi_q, 3);
% 
% IRF1_lo = prctile(IRF1_boot, lo_q, 3);
% IRF1_hi = prctile(IRF1_boot, hi_q, 3);
% 
% t = (0:s)';

% =========================
% PLOT: IRF0 vs IRF1 with bands (blue vs red)
% =========================
% figure;
% for j = 1:N
%     subplot(2,2,j); hold on;
% 
%     % --- BLUE band (Pre 2018) ---
%     patch([t; flipud(t)], [IRF0_lo(:,j); flipud(IRF0_hi(:,j))],'b', 'FaceAlpha', 0.18, 'EdgeColor', 'none');
% 
%     % --- RED band (Post 2018) ---
%     patch([t; flipud(t)], ...
%           [IRF1_lo(:,j); flipud(IRF1_hi(:,j))], ...
%           'r', 'FaceAlpha', 0.18, 'EdgeColor', 'none');
% 
%     % --- Lines ---
%     plot(t, Beta_0(:,j), 'b-', 'LineWidth', 2, 'DisplayName', 'Pre 2018 (D=0)');
%     plot(t, Beta_1(:,j), 'r-', 'LineWidth', 2, 'DisplayName', 'Post 2018 (D=1)');
% 
%     yline(0,'k--');
%     grid on;
%     title(var_names{j});
%     xlabel('Horizon');
%     ylabel('Response');
% 
%     if j == 1
%         legend('show','Location','best');
%     end
% end
% sgtitle('LP IRFs by regime with 68% bootstrap bands');




%%%---FINAL---%%%



pLP = min(optimal_p,4);   % consigliato 2 o 4
s = 40;

T = size(y,1);
N = size(y,2);

Beta_full = zeros(s+1,N);

for h = 0:s
    for j = 1:N
        
        Yh = y(optimal_p+1+h:T, j);
        
        % solo lag della variabile j
        Xlags = [];
        for L = 1:optimal_p
            Xlags = [Xlags y(optimal_p+1-L:T-L-h, j)];
        end
        
        Xh = [ones(size(Yh,1),1), ...
              z(optimal_p+1:T-h,:), ...
              Xlags];
        
        bh = (Xh'*Xh)\(Xh'*Yh);
        
        Beta_full(h+1,j) = bh(2);  
    end
end




Boot = 1000;
IRF_boot = zeros(s+1,N,Boot);

for h = 0:s
    for j = 1:N
        
        Yh = y(optimal_p+1+h:T, j);
        
        % regressori
        Xlags = [];
        for L = 1:optimal_p
            Xlags = [Xlags y(optimal_p+1-L:T-L-h, j)];
        end
        
        Xh = [ones(size(Yh,1),1), ...
              z(optimal_p+1:T-h,:), ...
              Xlags];
        
        % stima originale
        bh = (Xh'*Xh)\(Xh'*Yh);
        uhat = Yh - Xh*bh;
        
        T_eff = length(Yh);
        
        for b = 1:Boot
            
            % wild bootstrap Rademacher
            v = (rand(T_eff,1) > 0.5)*2 - 1;
            Y_star = Xh*bh + uhat .* v;
            
            bh_star = (Xh'*Xh)\(Xh'*Y_star);
            
            IRF_boot(h+1,j,b) = bh_star(2);
        end
    end
end

% percentile bands 68%
IRF_lo = prctile(IRF_boot,16,3);
IRF_hi = prctile(IRF_boot,84,3);


hgrid = 0:s;
t = hgrid';

figure;
for j = 1:N
    subplot(2,2,j); hold on;
    
    patch([t; flipud(t)], ...
          [IRF_lo(:,j); flipud(IRF_hi(:,j))], ...
          [0.7 0.7 0.7], ...
          'FaceAlpha',0.35,'EdgeColor','none');
    
    plot(hgrid, Beta_full(:,j), 'k','LineWidth',2);
    
    yline(0,'k--');
    grid on;
    title(var_names{j});
    xlabel('Horizon');
end

sgtitle('Local Projection IRFs – Oil Shock (Full Sample) with 68% Wild Bootstrap Bands');