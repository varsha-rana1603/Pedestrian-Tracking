function idx = find_init_window(stance, fs, min_len)
%FIND_INIT_WINDOW  First long stationary stretch (used for gravity + gyro bias).
if nargin < 3, min_len = 0.8; end
x = stance(:); d = diff([0; x; 0]); s = find(d==1); e = find(d==-1)-1;
L = e - s + 1;
k = find(L >= round(min_len*fs), 1, 'first');
if isempty(k), [~,k] = max(L); end        % fall back to the longest run
% trim 10% at both ends so transients are excluded
n = L(k); a = s(k) + round(0.1*n); b = e(k) - round(0.1*n);
idx = a:b;
end
