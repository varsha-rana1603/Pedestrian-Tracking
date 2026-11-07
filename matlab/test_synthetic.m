% Validate the whole pipeline on a synthetic walk with known ground truth.
if exist('OCTAVE_VERSION','builtin'), pkg load signal; end
s = sim_walk(100, 24, 1.3, 0);
[stance,T] = detect_stance(s.acc, s.gyr, s.fs);
init = find_init_window(stance, s.fs);
nv = ins_naive(s.acc, s.t, init);
kf = kf_zupt(s.acc - mean(s.acc(init,:),1), s.fs, stance);
fu = ins_fusion(s.acc, s.gyr, s.fs, stance, init);
true_end = s.pos_true(end,:);
fprintf('TRUE end position        : [%.2f %.2f %.2f] m (path %.1f m)\n', true_end, s.nstrides*s.L);
fprintf('Naive double integration : [%.1f %.1f %.1f] m\n', nv.p(end,:));
fprintf('Accel-only KF+ZUPT       : [%.2f %.2f %.2f] m\n', kf.p(end,:));
fprintf('Acc+Gyro fusion + ZUPT   : [%.2f %.2f %.2f] m\n', fu.p(end,:));
fprintf('stance fraction %.2f\n', mean(stance));
