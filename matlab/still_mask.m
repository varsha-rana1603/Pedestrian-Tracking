function m = still_mask(stance, fs, min_len)
%STILL_MASK  Samples that belong to long standing periods (stance runs >= min_len s),
%            trimmed by 0.15 s at both ends. Used to pool gyro-bias estimates.
if nargin < 3, min_len = 1.0; end
x = stance(:); d = diff([0; x; 0]); s = find(d==1); e = find(d==-1)-1;
m = false(size(x)); tr = round(0.15*fs);
for i = 1:numel(s)
    if (e(i)-s(i)+1) >= round(min_len*fs) && (e(i)-s(i)+1) > 2*tr
        m(s(i)+tr : e(i)-tr) = true;
    end
end
end
