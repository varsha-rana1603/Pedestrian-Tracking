clear; clc; close all;
 
stilltime     = 5;    
smoothing       = 0.10;    
thr_mag     = 0.8;     
thr_std     = 0.5;     

p_noise    = 0.5;    
R_zupt      = 1e-3;    

T  = readtable('CST05_left.csv');
sr = T{:,1};                 

T  = T(2001:end, :);
T  = T(T{:,1} > 1, :);
t_ = T.time;
T = T([true; diff(t_) > 0], :); 

t   = T.time - T.time(1);
A   = [T.ax, T.ay, T.az] * 9.81;    
N   = numel(t);
dtv = [median(diff(t)); diff(t)];  
fs  = 1/median(diff(t));
fprintf('Samples: %d, duration %.1f s, ~%.1f Hz\n', N, t(end), fs);

idx0   = t < 16;
a_bias = mean(A(idx0,:), 1);      
a_lin  = A - a_bias;

v_ = cumtrapz(t, a_lin);        
p_ = cumtrapz(t, v_) ;

mag    = vecnorm(A, 2, 2);
w      = max(3, round(smoothing * fs));
mag_s  = movmean(mag, w);
mag_sd = movstd(mag, w);
stationary = (abs(mag_s - 9.81) < thr_mag) & (mag_sd < thr_std);
fprintf('Stationary samples: %.1f %%\n', 100*mean(stationary));

p_kf = zeros(N, 3);  v_kf = zeros(N, 3);
for ax = 1:3
    [p_kf(:,ax), v_kf(:,ax)] = kf_axis(a_lin(:,ax), stationary, dtv, p_noise, R_zupt);
end
p_kf = p_kf * 1000;               

fprintf('Final position without KF (mm): [%.0f %.0f %.0f]\n', p_(end,:));
fprintf('Final position with KF    (mm): [%.0f %.0f %.0f]\n', p_kf(end,:));

labels = {'x','y','z'};

figure('Name','Stationary detection');
plot(t, mag, 'Color', [0.7 0.7 0.7]); hold on;
plot(t, 9.81*stationary, 'r'); plot([0 t(end)], [9.81 9.81], 'k--');
xlabel('Time (s)'); ylabel('|a| (m/s^2)'); grid on;
legend('|a|','stationary flag (scaled)','g'); title('Detected stance phases');

figure('Name','Part (a)');
for ax = 1:3
    subplot(3,1,ax); plot(t, p_(:,ax), 'LineWidth', 1.2); grid on;
    ylabel(sprintf('%s (mm)', labels{ax}));
    if ax == 1, title('Part (a): position without Kalman filtering'); end
end
xlabel('Time (s)');

figure('Name','Part (b)');
for ax = 1:3
    subplot(3,1,ax); plot(t, p_kf(:,ax), 'LineWidth', 1.2); grid on;
    ylabel(sprintf('%s (mm)', labels{ax}));
    if ax == 1, title('Part (b): position with Kalman filter (ZUPT)'); end
end
xlabel('Time (s)');

figure('Name','Comparison');
for ax = 1:3
    subplot(3,1,ax);
    plot(t, p_(:,ax), 'r', t, p_kf(:,ax), 'b', 'LineWidth', 1.2); grid on;
    ylabel(sprintf('%s (mm)', labels{ax}));
    if ax == 1
        title('Raw integration vs. Kalman filter');
        legend('Without KF','With KF (ZUPT)','Location','northwest');
    end
end
xlabel('Time (s)');

figure('Name','Velocity');
for ax = 1:3
    subplot(3,1,ax);
    plot(t, v_(:,ax), 'r', t, v_kf(:,ax), 'b', 'LineWidth', 1.0); grid on;
    ylabel(sprintf('v_%s (m/s)', labels{ax}));
    if ax == 1
        title('Velocity: drifts without KF, reset at each footfall with KF');
        legend('Without KF','With KF (ZUPT)','Location','northwest');
    end
end
xlabel('Time (s)');

figure('Name','XY path');
plot(p_kf(:,1), p_kf(:,2), 'b', 'LineWidth', 1.2); hold on;
plot(0, 0, 'ko', 'MarkerFaceColor', 'k');
axis equal; grid on; xlabel('x (mm)'); ylabel('y (mm)');
title('Estimated horizontal path with Kalman filter');

function [p, v] = kf_axis(acc, stationary, dtv, sigma_a, R)

    N = numel(acc);
    H = [0 1];
    x = [0; 0];  P = zeros(2);
    p = zeros(N,1);  v = zeros(N,1);
    for k = 1:N
        dt = dtv(k);
        F = [1 dt; 0 1];
        B = [0.5*dt^2; dt];
        Q = sigma_a^2 * [dt^4/4, dt^3/2; dt^3/2, dt^2];
        x = F*x + B*acc(k);       
        P = F*P*F' + Q;
        if stationary(k)          
            S = H*P*H' + R;
            K = P*H' / S;
            x = x + K*(0 - H*x);
            P = (eye(2) - K*H) * P;
        end
        p(k) = x(1);  v(k) = x(2);
    end
end
