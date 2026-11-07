function out = ins_fusion(acc, gyr, fs, stance, init, prm)
%INS_FUSION  Part (c) Approach 2: accelerometer + gyroscope fusion with ZUPT.
%   1) gyro bias from the initial still window
%   2) Mahony complementary filter: gyro integration (quaternion), accelerometer
%      tilt correction applied ONLY during stance (when |a| ~ g is trustworthy)
%   3) rotate specific force to the world frame (z up), subtract gravity
%   4) Kalman filter + ZUPT (kf_zupt) on the world-frame acceleration
%   Heading (yaw) is unobservable without a magnetometer -> starts at 0 and drifts.
if nargin < 6, prm = struct(); end
Kp = gf(prm,'Kp_stance',8);
g = mean(sqrt(sum(acc(init,:).^2,2)));      % local |g| from the still window (removes scale error)
N = size(acc,1); dt = 1/fs;
bg = mean(gyr(init,:),1);  w = gyr - bg;                  % bias-corrected gyro
% initial orientation: align mean specific force with world +z
a0 = mean(acc(init,:),1)'; a0 = a0/norm(a0); z = [0;0;1];
ax = cross(a0,z); s = norm(ax); c = dot(a0,z);
if s < 1e-9, q = [1 0 0 0]'; else, ang = atan2(s,c); ax = ax/s; q = [cos(ang/2); ax*sin(ang/2)]; end
alin = zeros(N,3); Rall = zeros(3,3,N); Q4 = zeros(N,4);
for k = 1:N
    R = quat2rot(q);  Rall(:,:,k) = R;  Q4(k,:) = q';
    alin(k,:) = (R*acc(k,:)' - [0;0;g])';
    % ---- Mahony update
    om = w(k,:)';
    if stance(k)
        an = acc(k,:)'/norm(acc(k,:));  u = R'*z;
        om = om + Kp*cross(an,u);
    end
    q = q + 0.5*dt*qmult(q,[0; om]);  q = q/norm(q);
end
kf = kf_zupt(alin, fs, stance, prm);
out.alin = alin; out.p = kf.p; out.v = kf.v; out.v_pre = kf.v_pre; out.q = Q4;
out.R = Rall; out.gyro_bias = bg;
% yaw rate about the world vertical, and heading of the sensor x/y axis
wz = zeros(N,1);
for k = 1:N, wz(k) = Rall(3,:,k)*gyr(k,:)' - Rall(3,:,k)*bg'; end
out.yaw_rate = wz; out.yaw = cumsum(wz)*dt;
end

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
function v = gf(p,f,d)
if isfield(p,f), v = p.(f); else v = d; end
end
