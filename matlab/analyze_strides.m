function S = analyze_strides(p, v, stance, fs, yaw_rate)
%ANALYZE_STRIDES  Stride-wise statistics from a position track and stance flags.
%   Stance run mid-points are the stride anchors. Two stances belong to the same walking
%   bout if the swing between them is shorter than 1 s.
x = stance(:); d = diff([0; x; 0]); s = find(d==1); e = find(d==-1)-1;
n = numel(s); mid = zeros(n,3);
for i = 1:n, mid(i,:) = mean(p(s(i):e(i),:),1); end
gap = (s(2:end) - e(1:end-1)) / fs;                % swing duration between stances [s]
dur = (e - s + 1) / fs;                            % stance durations [s]
walk = gap < 1.0 & dur(1:end-1) < 1.2 & dur(2:end) < 1.2;  % long stance = standing, breaks a bout
dp = mid(2:end,:) - mid(1:end-1,:);
S.dp = dp(walk,:);  S.len = sqrt(sum(S.dp(:,1:2).^2,2));  S.dz = S.dp(:,3);
S.len3 = sqrt(sum(S.dp.^2,2));
S.t_mid = (s(2:end)+e(2:end))/2/fs;  S.t_mid = S.t_mid(walk);
S.heading = unwrap(atan2(S.dp(:,2), S.dp(:,1)));
S.n_strides = numel(S.len);
% bouts
b_start = [1; find(~walk)+1]; b_end = [find(~walk); n-1];
S.bouts = [];
for k = 1:numel(b_start)
    ii = b_start(k):b_end(k);                       % stride indices (in 'dp' indexing)
    ii = ii(walk(ii));
    if numel(ii) < 3, continue; end
    pp = [mid(ii(1),:); mid(ii+1,:)];
    h  = unwrap(atan2(diff(pp(:,2)), diff(pp(:,1))));
    S.bouts = [S.bouts; (s(ii(1)))/fs, e(ii(end)+1)/fs, numel(ii), ...
        sum(sqrt(sum(diff(pp).^2,2))), norm(pp(end,1:2)-pp(1,1:2)), ...
        rad2deg(h(end)-h(1)), pp(end,3)-pp(1,3)]; %#ok<AGROW>
end
% columns: t_start t_end n_strides path_len net_horiz_disp heading_change_deg net_dz
S.mid = mid;
end
