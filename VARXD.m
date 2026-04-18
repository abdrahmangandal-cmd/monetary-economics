clear all

data=xlsread('MONETARY DATASET.xlsx','ALL DATA','B2:K223');

y = [diff(log(data(1:end,[1 8 5])))*100 data(2:end,4)];
z = data(2:end,3);
dummy = data(2:end,10);

var_names = {'PRICE LEVEL (cum. inflation)', 'RESERVES (log level, cum.)', 'REER (log level, cum.)', 'MMR'};

max_p = 12;

AIC = zeros(1,max_p);

for p = 1:max_p
    [beta,res,T,N,X] = OLS(y,p,z,dummy);
    T_eff = size(res,1);
    Sigma = res.'*res/T_eff;
    AIC(p) = log(det(Sigma)) + 2*(size(X,2)*N)/T_eff;  
end


[best_AIC, optimal_p] = min(AIC);
fprintf('Optimal lag for VAR: p = %d (AIC = %.4f)\n', optimal_p, best_AIC);
 

[A_best, res_best, T,N,X] = OLS(y,optimal_p,z,dummy);




Np = N*optimal_p;


A0 = A_best(2:1+Np, :);  

C0 = zeros(Np,Np);
C0(1:N,1:Np) = reshape(A0.', N, Np);
C0(N+1:end,1:end-N) = eye(Np-N);


Aint = A_best(2+Np:1+2*Np, :);  
A1 = A0 + Aint;

C1 = zeros(Np,Np);
C1(1:N,1:Np) = reshape(A1.', N, Np);
C1(N+1:end,1:end-N) = eye(Np-N);



H = 40;

J = [eye(N), zeros(N, N*(optimal_p-1))];


C0_h = zeros(N, N, H+1);
C1_h = zeros(N, N, H+1);

C0_h(:,:,1) = eye(N);
C1_h(:,:,1) = eye(N);

for h = 1:H
    C0_h(:,:,h+1) = J * (C0^h) * J.';
    C1_h(:,:,h+1) = J * (C1^h) * J.';
end




Np = N*optimal_p;
k_z  = 2*Np + 3;  
k_zD = 2*Np + 4;  

B0 = A_best(k_z , :).';                   
B1 = A_best(k_z , :).'+ A_best(k_zD,:).'; 

IRF0 = zeros(N, H+1);
IRF1 = zeros(N, H+1);

for h = 0:H
    IRF0(:,h+1) = C0_h(:,:,h+1) * B0;
    IRF1(:,h+1) = C1_h(:,:,h+1) * B1;
end


cum_idx = 1:3;                
IRF0(cum_idx,:) = cumsum(IRF0(cum_idx,:), 2);
IRF1(cum_idx,:) = cumsum(IRF1(cum_idx,:), 2);


Boot = 1000;

IRF0_boot = zeros(N, H+1, Boot);
IRF1_boot = zeros(N, H+1, Boot);

T_total = size(y,1);
T_eff = T_total - optimal_p;

for b = 1:Boot
    
    
    YB = zeros(T_total,N);
    YB(1:optimal_p,:) = y(1:optimal_p,:);
    
    for t = optimal_p+1:T_total
        
        res_ind = randi(T_eff);
        eps_b = res_best(res_ind,:);
        
        
        Xb = [];
        for lag = 1:optimal_p
            Xb = [Xb YB(t-lag,:)];
        end
        
        d_t = dummy(t);
        z_t = z(t);
        
        Xb = [1 Xb Xb.*d_t d_t z_t z_t*d_t];
        
        YB(t,:) = Xb * A_best + eps_b;
    end
    
    
    [A_boot, ~, ~, ~, ~] = OLS(YB,optimal_p,z,dummy);
    
    
    Np = N*optimal_p;
    
    A0_boot = A_boot(2:1+Np,:);
    Aint_boot = A_boot(2+Np:1+2*Np,:);
    A1_boot = A0_boot + Aint_boot;
    
    
    C0_boot = zeros(Np,Np);
    C0_boot(1:N,1:Np) = reshape(A0_boot.',N,Np);
    C0_boot(N+1:end,1:end-N) = eye(Np-N);
    
    
    C1_boot = zeros(Np,Np);
    C1_boot(1:N,1:Np) = reshape(A1_boot.',N,Np);
    C1_boot(N+1:end,1:end-N) = eye(Np-N);
    
    
    k_z  = 2*Np + 3;
    k_zD = 2*Np + 4;
    
    B0_boot = A_boot(k_z , :).';
    B1_boot = A_boot(k_z , :).'+ A_boot(k_zD,:).';
    
    
    J = [eye(N), zeros(N, N*(optimal_p-1))];
    
    for h = 0:H
        
        Psi0 = J*(C0_boot^h)*J.';
        Psi1 = J*(C1_boot^h)*J.';
        
        IRF0_boot(:,h+1,b) = Psi0 * B0_boot;
        IRF1_boot(:,h+1,b) = Psi1 * B1_boot;
    end
    
end

cum_idx = 1:3;
IRF0_boot(cum_idx,:,:) = cumsum(IRF0_boot(cum_idx,:,:), 2);
IRF1_boot(cum_idx,:,:) = cumsum(IRF1_boot(cum_idx,:,:), 2);


IRF0_low = prctile(IRF0_boot,16,3);
IRF0_up  = prctile(IRF0_boot,84,3);

IRF1_low = prctile(IRF1_boot,16,3);
IRF1_up  = prctile(IRF1_boot,84,3);

t_axis = (0:H)';

figure;

for i = 1:N
    subplot(2,2,i); hold on;

    
    lo0 = squeeze(IRF0_low(i,:))';
    hi0 = squeeze(IRF0_up(i,:))';
    
    patch([t_axis; flipud(t_axis)], [lo0; flipud(hi0)], [0.6 0.8 1], ...
        'FaceAlpha',0.35, 'EdgeColor','none', 'HandleVisibility','off');

    plot(t_axis, IRF0(i,:)', 'Color',[0 0.2 0.8],'LineWidth',1.8);

    
    lo1 = squeeze(IRF1_low(i,:))';
    hi1 = squeeze(IRF1_up(i,:))';

    patch([t_axis; flipud(t_axis)], [lo1; flipud(hi1)], [1 0.7 0.7], ...
        'FaceAlpha',0.35, 'EdgeColor','none', 'HandleVisibility','off');

    plot(t_axis, IRF1(i,:)', 'Color',[0.8 0 0], 'LineWidth',1.8);

    yline(0,'k--', 'HandleVisibility','off')  % Anche la linea zero esclusa
    grid on
    title(var_names{i})

    if i==1
        legend('Pre-2018','Post-2018','Location','best')
    end
end

sgtitle('IRF Oil Kanzig Shock – Regime-dependent VARXD');





