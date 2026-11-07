function out = kf_zupt(a, fs, stance, prm)
%KF_ZUPT  Linear Kalman filter, state x=[p(3); v(3)], acceleration as control input,
%         zero-velocity pseudo-measurement z=0 during stance (ZUPT).
%   a      : N x 3 linear acceleration in the frame where p,v are estimated [m/s^2]
%   stance : N x 1 logical
%   prm    : .sig_a (process, m/s^2)  .sig_v (ZUPT meas. noise, m/s)
if nargin < 4, prm = struct(); end
sig_a = gf(prm,'sig_a',0.3);  sig_v = gf(prm,'sig_v',0.01);
N = size(a,1); dt = 1/fs; I3 = eye(3);
F = [I3, dt*I3; zeros(3), I3];
B = [0.5*dt^2*I3; dt*I3];
Q = blkdiag((0.5*dt^2*sig_a)^2*I3, (dt*sig_a)^2*I3);
H = [zeros(3), I3];  R = sig_v^2*I3;
x = zeros(6,1); P = 1e-6*eye(6);
X = zeros(N,6); Vpre = zeros(N,3); Pd = zeros(N,6);
for k = 1:N
    x = F*x + B*a(k,:)';
    P = F*P*F' + Q;
    Vpre(k,:) = x(4:6)';
    if stance(k)
        S = H*P*H' + R;  K = (P*H')/S;
        x = x + K*(0 - H*x);
        IKH = eye(6) - K*H;  P = IKH*P*IKH' + K*R*K';   % Joseph form
    end
    X(k,:) = x';  Pd(k,:) = diag(P)';
end
out.p = X(:,1:3); out.v = X(:,4:6); out.v_pre = Vpre; out.Pdiag = Pd;
end
function v = gf(p, f, d)
if isfield(p,f), v = p.(f); else v = d; end
end
