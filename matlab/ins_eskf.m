function out = ins_eskf(acc, gyr, fs, stance, init, prm)
%INS_ESKF  Part (c) Approach 2 (full): strapdown INS + error-state Kalman filter with ZUPT.
%   Nominal states : orientation quaternion q (sensor->world, z up), velocity v, position p
%   Error state    : dx = [dp; dv; dphi] (9x1), world-frame orientation error dphi
%       dp' = dv ;   dv' = -[f_w]x dphi + noise ;   dphi' = noise
%       with f_w = R*f (specific force in the world frame)
%   ZUPT: z = 0 - v (pseudo-measurement) while stance(k) is true.
%   The cross-covariance between dv and dphi lets each ZUPT also correct tilt
%   (roll/pitch). Heading (yaw) stays unobservable (no magnetometer) and drifts.
if nargin < 6, prm = struct(); end
sa  = gf(prm,'sig_a',0.5);                 % accel noise density     [m/s^2/sqrt(Hz)]
sg  = gf(prm,'sig_g',0.05);                % gyro  noise density     [rad/s/sqrt(Hz)]
sv  = gf(prm,'sig_v',0.03);                % ZUPT measurement std    [m/s]
N = size(acc,1); dt = 1/fs;
g = mean(sqrt(sum(acc(init,:).^2,2)));     % local |g| from the still window
bg = mean(gyr(init,:),1);                  % gyro bias: still window ...
if isfield(prm,'bias_idx') && nnz(prm.bias_idx) > fs   % ... or pooled over all standing periods
    bg = mean(gyr(prm.bias_idx,:),1);
end
w = gyr - bg;
% --- initial attitude: align mean specific force with world +z
a0 = mean(acc(init,:),1)'; a0 = a0/norm(a0); z = [0;0;1];
ax = cross(a0,z); s = norm(ax); c = dot(a0,z);
if s < 1e-9, q = [1;0;0;0]; else, ang = atan2(s,c); ax = ax/s; q = [cos(ang/2); ax*sin(ang/2)]; end
v = zeros(3,1); p = zeros(3,1);
P = blkdiag(1e-6*eye(3), 1e-6*eye(3), (0.5*pi/180)^2*eye(3));
Qd = blkdiag(zeros(3), sa^2*dt*eye(3), sg^2*dt*eye(3));
H = [zeros(3), eye(3), zeros(3)];  Rm = sv^2*eye(3);
Pp = zeros(N,3); Vv = zeros(N,3); Vpre = zeros(N,3); Q4 = zeros(N,4); Aw = zeros(N,3);
yawrate = zeros(N,1); Rall = zeros(3,3,N);
for k = 1:N
    R  = quat2rot(q);  fw = R*acc(k,:)';  aw = fw - [0;0;g];
    Rall(:,:,k) = R;  Aw(k,:) = aw';  yawrate(k) = R(3,:)*w(k,:)';
    % ---- nominal propagation (position/velocity with current a, then attitude)
    p = p + v*dt + 0.5*aw*dt^2;  v = v + aw*dt;
    q = q + 0.5*dt*qmult(q,[0; w(k,:)']);  q = q/norm(q);
    % ---- error-state covariance propagation
    F = eye(9); F(1:3,4:6) = dt*eye(3); F(4:6,7:9) = -dt*skew(fw);
    P = F*P*F' + Qd;
    Vpre(k,:) = v';
    % ---- ZUPT update
    if stance(k)
        S = H*P*H' + Rm;  K = (P*H')/S;
        dx = K*(0 - v);
        p = p + dx(1:3);  v = v + dx(4:6);
        dphi = dx(7:9);   dq = [1; 0.5*dphi]; dq = dq/norm(dq);
        q = qmult(dq, q); q = q/norm(q);
        IKH = eye(9) - K*H;  P = IKH*P*IKH' + K*Rm*K';  P = (P+P')/2;
    end
    Pp(k,:) = p'; Vv(k,:) = v'; Q4(k,:) = q';
end
out.p = Pp; out.v = Vv; out.v_pre = Vpre; out.q = Q4; out.alin = Aw; out.R = Rall;
out.gyro_bias = bg; out.yaw_rate = yawrate; out.yaw = cumsum(yawrate)*dt;
end

function S = skew(a), S = [0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0]; end
function R = quat2rot(q)
w=q(1); x=q(2); y=q(3); z=q(4);
R = [1-2*(y^2+z^2), 2*(x*y-z*w),   2*(x*z+y*w);
     2*(x*y+z*w),   1-2*(x^2+z^2), 2*(y*z-x*w);
     2*(x*z-y*w),   2*(y*z+x*w),   1-2*(x^2+y^2)];
end
function r = qmult(a,b)
r = [a(1)*b(1)-a(2)*b(2)-a(3)*b(3)-a(4)*b(4);
     a(1)*b(2)+a(2)*b(1)+a(3)*b(4)-a(4)*b(3);
     a(1)*b(3)-a(2)*b(4)+a(3)*b(1)+a(4)*b(2);
     a(1)*b(4)+a(2)*b(3)-a(3)*b(2)+a(4)*b(1)];
end
function v = gf(p,f,d), if isfield(p,f), v = p.(f); else v = d; end, end
