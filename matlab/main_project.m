function R = main_project(file, outdir, opts)
%MAIN_PROJECT  Problem 1 (a)(b)(c) for one IMU file. Saves figures + returns metrics.
%   R = main_project('CST03_left.csv', 'results')
if exist('OCTAVE_VERSION','builtin'), pkg load signal; end
if nargin < 2, outdir = 'results'; end
if nargin < 3, opts = struct(); end
if ~exist(outdir,'dir'), mkdir(outdir); end
[~,tag] = fileparts(file); g = 9.80665;

%% 0) Load + stance detection + initial still window
d = load_imu(file);  t = d.t;  fs = d.fs;
[stance, ~] = detect_stance(d.acc, d.gyr, fs);
init = find_init_window(stance, fs);
fprintf('\n=== %s (%s): %.0f s @ %d Hz, stance %.0f%%, init window %.2f-%.2f s\n', ...
    tag, d.fmt, t(end), fs, 100*mean(stance), t(init(1)), t(init(end)));

%% a) Naive double integration (sensor frame, constant gravity removal)
nv = ins_naive(d.acc, t, init);

%% b) Kalman filter + ZUPT on the same accelerometer-only pipeline
a_lin = d.acc - mean(d.acc(init,:),1);
kf = kf_zupt(a_lin, fs, stance);

%% c) Approach 1 = accelerometer only (= b); Approach 2 = accelerometer + gyroscope
bias_idx = still_mask(stance, fs, 1.0);          % pool gyro bias over all standing periods
fu = ins_eskf(d.acc, d.gyr, fs, stance, init, struct('bias_idx', bias_idx));

%% Metrics
S1 = analyze_strides(kf.p, kf.v, stance, fs);          % accel-only
S2 = analyze_strides(fu.p, fu.v, stance, fs);          % fusion
Sn = analyze_strides(nv.p, nv.v, stance, fs);
x = stance(:); dd = diff([0;x;0]); st0 = find(dd==1);
R.file = tag; R.fmt = d.fmt; R.dur = t(end); R.fs = fs; R.stance_frac = mean(stance);
R.naive_max_pos = max(sqrt(sum(nv.p.^2,2)));   R.naive_end_pos = norm(nv.p(end,:));
R.naive_v_at_stance_rms = sqrt(mean(sum(nv.v(stance,:).^2,2)));
R.kf_v_pre_rms   = sqrt(mean(sum(kf.v_pre(st0,:).^2,2)));   % velocity drift accumulated per stride
R.fu_v_pre_rms   = sqrt(mean(sum(fu.v_pre(st0,:).^2,2)));
R.kf_max_pos = max(sqrt(sum(kf.p.^2,2)));      R.fu_max_pos = max(sqrt(sum(fu.p.^2,2)));
R.kf_stride_mean = mean(S1.len);  R.kf_stride_std = std(S1.len);
R.fu_stride_mean = mean(S2.len);  R.fu_stride_std = std(S2.len);
R.fu_stride_med  = median(S2.len);
R.n_strides = S2.n_strides;
R.fu_path = sum(S2.len3);  R.kf_path = sum(S1.len3);
R.fu_mean_abs_dz = mean(abs(S2.dz));  R.fu_max_abs_z = max(abs(fu.p(:,3)));
R.fu_net_dz = sum(S2.dz);
R.fu_bouts = S2.bouts;  R.kf_bouts = S1.bouts;
R.gyro_bias_dps = fu.gyro_bias*180/pi;  R.info = d.info;
R.stride_time_med = median(diff(S2.t_mid));
R.fu_end_to_start = norm(fu.p(end,1:2) - fu.p(1,1:2));   % loop-closure error if the path is closed
R.fu_end_z = fu.p(end,3);  R.kf_end_to_start = norm(kf.p(end,1:2));

%% -------- Figures
set(0,'DefaultAxesFontSize',9);
col = struct('naive',[0.85 0.33 0.10],'kf',[0 0.45 0.74],'fu',[0.47 0.67 0.19]);
% window used for zooms: first 25 s of walking
mv = find(movmean(sqrt(sum(d.gyr.^2,2)),fs) > 1.0, 1);  if isempty(mv), mv = 1; end
z0 = max(1, mv - round(2*fs)); z1 = min(numel(t), z0 + round(25*fs)); zi = z0:z1;
lab = {'x','y','z'};

% ---- Fig 1: part (a)
f = figure('Visible','off','Position',[50 50 1200 800]);
subplot(3,2,1); plot(t, d.acc/g); grid on; xlabel('t [s]'); ylabel('a [g]'); title('Raw accelerometer (sensor frame)'); legend(lab,'Location','northeast');
subplot(3,2,2); plot(t, d.gyr*180/pi); grid on; xlabel('t [s]'); ylabel('\omega [deg/s]'); title('Raw gyroscope (not used in (a),(b))');
subplot(3,2,3); plot(t, nv.p); grid on; xlabel('t [s]'); ylabel('position [m]'); title('(a) Double integration: full record'); legend(lab,'Location','northwest');
subplot(3,2,4); plot(t(zi), nv.p(zi,:)-nv.p(zi(1),:)); grid on; xlabel('t [s]'); ylabel('position [m]'); title('(a) zoom: first 25 s of walking');
subplot(3,2,5); plot(t, nv.v); grid on; xlabel('t [s]'); ylabel('velocity [m/s]'); title('(a) velocity: never returns to 0 at stance'); hold on; plot(t(stance), zeros(nnz(stance),1),'k.','MarkerSize',2);
subplot(3,2,6); plot3(nv.p(zi,1)-nv.p(zi(1),1), nv.p(zi,2)-nv.p(zi(1),2), nv.p(zi,3)-nv.p(zi(1),3),'Color',col.naive); grid on; xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]'); title('(a) 3-D trajectory (25 s zoom)'); view(35,25);
print(f,'-dpng','-r100',fullfile(outdir,[tag '_a_naive.png'])); close(f);

% ---- Fig 2: part (b)
f = figure('Visible','off','Position',[50 50 1200 800]);
for ax = 1:3
    subplot(3,3,ax); plot(t(zi), nv.p(zi,ax)-nv.p(zi(1),ax),'Color',col.naive); hold on; plot(t(zi), kf.p(zi,ax)-kf.p(zi(1),ax),'Color',col.kf,'LineWidth',1.2);
    grid on; xlabel('t [s]'); ylabel([lab{ax} ' [m]']); title(['(b) position ' lab{ax} ' (25 s zoom)']); if ax==1, legend('no KF','KF+ZUPT','Location','best'); end
end
for ax = 1:3
    subplot(3,3,3+ax); plot(t(zi), nv.v(zi,ax),'Color',col.naive); hold on; plot(t(zi), kf.v(zi,ax),'Color',col.kf);
    grid on; xlabel('t [s]'); ylabel(['v_' lab{ax} ' [m/s]']); title(['velocity ' lab{ax}]);
end
subplot(3,3,7); plot(t, sqrt(sum(nv.p.^2,2)),'Color',col.naive); hold on; plot(t, sqrt(sum(kf.p.^2,2)),'Color',col.kf); set(gca,'YScale','log'); grid on; xlabel('t [s]'); ylabel('|p| [m] (log)'); title('position magnitude, full record');
subplot(3,3,8); plot(t(zi), stance(zi)*1,'k'); hold on; plot(t(zi), sqrt(sum(d.gyr(zi,:).^2,2))/20,'Color',[.6 .6 .6]); grid on; xlabel('t [s]'); title('stance flag (ZUPT active) and |\omega|/20'); ylim([-0.1 1.2]);
subplot(3,3,9); plot(kf.p(zi,1)-kf.p(zi(1),1), kf.p(zi,2)-kf.p(zi(1),2),'Color',col.kf); grid on; axis equal; xlabel('x [m]'); ylabel('y [m]'); title('(b) KF+ZUPT x-y path (sensor-frame axes)');
print(f,'-dpng','-r100',fullfile(outdir,[tag '_b_kalman.png'])); close(f);

% ---- Fig 3: part (c)
f = figure('Visible','off','Position',[50 50 1200 800]);
subplot(2,3,1); plot(kf.p(:,1), kf.p(:,2),'Color',col.kf); grid on; axis equal; xlabel('x [m]'); ylabel('y [m]'); title('Approach 1: accel only (KF+ZUPT)');
subplot(2,3,2); plot(fu.p(:,1), fu.p(:,2),'Color',col.fu); hold on; plot(fu.p(1,1),fu.p(1,2),'ko','MarkerFaceColor','k'); plot(fu.p(end,1),fu.p(end,2),'rs','MarkerFaceColor','r'); grid on; axis equal; xlabel('x [m]'); ylabel('y [m]'); title('Approach 2: accel + gyro (ESKF + ZUPT)'); legend('path','start','end','Location','best');
subplot(2,3,3); plot3(fu.p(:,1), fu.p(:,2), fu.p(:,3),'Color',col.fu); grid on; xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]'); title('Approach 2: 3-D trajectory'); view(35,25);
subplot(2,3,4); plot(t, kf.p(:,3),'Color',col.kf); hold on; plot(t, fu.p(:,3),'Color',col.fu); grid on; xlabel('t [s]'); ylabel('z [m]'); title('Vertical coordinate'); legend('accel only','fusion','Location','best');
subplot(2,3,5); plot(t, sqrt(sum(kf.p(:,1:2).^2,2)),'Color',col.kf); hold on; plot(t, sqrt(sum(fu.p(:,1:2).^2,2)),'Color',col.fu); grid on; xlabel('t [s]'); ylabel('range from start [m]'); title('Horizontal distance from start'); legend('accel only','fusion','Location','best');
subplot(2,3,6); plot(t, rad2deg(fu.yaw),'Color',col.fu); grid on; xlabel('t [s]'); ylabel('\psi [deg]'); title('Fusion: integrated yaw rate (heading drifts, no magnetometer)');
print(f,'-dpng','-r100',fullfile(outdir,[tag '_c_fusion.png'])); close(f);

% ---- Fig 4: scenario identification
f = figure('Visible','off','Position',[50 50 1200 800]);
subplot(2,3,1); plot(t, sqrt(sum(d.gyr.^2,2))*180/pi,'Color',[.3 .3 .3]); hold on; plot(t, 300*stance,'r'); grid on; xlabel('t [s]'); ylabel('|\omega| [deg/s]'); title('Gait activity and detected stance (red)');
subplot(2,3,2); hist(S2.len, 20); grid on; xlabel('stride length [m]'); ylabel('count'); title(sprintf('Stride length  %.2f \\pm %.2f m (median %.2f)', mean(S2.len), std(S2.len), median(S2.len)));
subplot(2,3,3); plot(S2.t_mid, S2.dz,'o-','MarkerSize',3); grid on; xlabel('t [s]'); ylabel('\Deltaz per stride [m]'); title('Height change per stride (stairs \Rightarrow consistent sign)');
subplot(2,3,4); plot(S2.t_mid, rad2deg(S2.heading),'.-'); grid on; xlabel('t [s]'); ylabel('stride heading [deg]'); title('Heading of successive strides (turns \Rightarrow steps)');
subplot(2,3,5); plot(S2.t_mid, S2.len,'.-'); grid on; xlabel('t [s]'); ylabel('stride length [m]'); title('Stride length vs time');
subplot(2,3,6); plot(t, fu.yaw_rate*180/pi); grid on; xlabel('t [s]'); ylabel('deg/s'); title('Yaw rate about vertical');
print(f,'-dpng','-r100',fullfile(outdir,[tag '_d_scenario.png'])); close(f);

R.t = t; R.nv = nv; R.kf = kf; R.fu = fu; R.stance = stance; R.S2 = S2;
save('-mat', fullfile(outdir,[tag '_results.mat']), 'R');
print_summary(R);
end

function print_summary(R)
fprintf('  (a) naive: max |p| = %.0f m, final |p| = %.0f m, RMS v at stance = %.2f m/s\n', R.naive_max_pos, R.naive_end_pos, R.naive_v_at_stance_rms);
fprintf('  (b) KF+ZUPT accel-only: max |p| = %.0f m, v drift per stride (RMS before ZUPT) = %.2f m/s, mean stride %.2f +- %.2f m\n', R.kf_max_pos, R.kf_v_pre_rms, R.kf_stride_mean, R.kf_stride_std);
fprintf('  (c) fusion: strides = %d, stride = %.2f +- %.2f m (median %.2f), v drift/stride = %.2f m/s, path = %.1f m, max|z| = %.2f m, mean|dz|/stride = %.3f m\n', ...
    R.n_strides, R.fu_stride_mean, R.fu_stride_std, R.fu_stride_med, R.fu_v_pre_rms, R.fu_path, R.fu_max_abs_z, R.fu_mean_abs_dz);
disp('  fusion bouts [t0 t1 nStr path net_disp dHeading_deg dz]:'); disp(round(100*R.fu_bouts)/100);
end
