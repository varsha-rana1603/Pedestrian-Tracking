function [stance, T] = detect_stance(acc, gyr, fs, p)
%DETECT_STANCE  SHOE-style stance (zero-velocity) detector.
%   T(k) = mean over a window of  |w|^2/sig_g^2 + |a - g*a_mean/|a_mean||^2/sig_a^2
%   stance(k) = T(k) < gamma   (foot flat on the ground), short blips removed.
if nargin < 4, p = struct(); end
sig_a = gf(p,'sig_a',0.6);   % m/s^2
sig_g = gf(p,'sig_g',0.35);  % rad/s
gamma = gf(p,'gamma', 30);
W     = max(3, round(gf(p,'win',0.1)*fs));
minst = round(gf(p,'min_stance',0.06)*fs);   % shortest stance accepted [s]
g = 9.80665;
am  = movmean(acc, W, 1);
res = acc - g * am ./ max(sqrt(sum(am.^2,2)), eps);
T   = (sum(res.^2,2)/sig_a^2 + sum(gyr.^2,2)/sig_g^2);
T   = movmean(T, W);
stance = T < gamma;
stance = clean_runs(stance, minst);          % drop too-short stance runs
stance = ~clean_runs(~stance, round(0.02*fs)); % fill tiny gaps
end

function x = clean_runs(x, nmin)
x = x(:); d = diff([0; x; 0]); s = find(d==1); e = find(d==-1)-1;
for i = 1:numel(s)
    if (e(i)-s(i)+1) < nmin, x(s(i):e(i)) = false; end
end
end

function v = gf(p, f, dflt)
if isfield(p, f), v = p.(f); else v = dflt; end
end
