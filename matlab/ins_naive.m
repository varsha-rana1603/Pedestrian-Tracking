function out = ins_naive(acc, t, init)
%INS_NAIVE  Part (a): remove a constant gravity offset, integrate twice (sensor frame).
%   No orientation tracking, no filtering: a_lin = a - mean(a_static).
a0   = mean(acc(init,:), 1);
alin = acc - a0;
v = cumtrapz(t, alin);
p = cumtrapz(t, v);
out.a = alin; out.v = v; out.p = p;
end
