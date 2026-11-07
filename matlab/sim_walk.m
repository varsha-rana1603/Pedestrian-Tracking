function s = sim_walk(fs, nstrides, L, turn_deg)
%SIM_WALK  Synthetic foot-mounted IMU with KNOWN ground truth (for validating the code).
%   Stance 0.55 s, swing 0.45 s. Foot moves L [m] per stride in swing (min-jerk),
%   lifts 0.12 m, pitches +-35 deg about the lateral axis. Optional turn per stride.
%   Sensor axes (like CST data): x = vertical when foot flat, y = forward, z = lateral.
if nargin<1, fs=100; end, if nargin<2, nstrides=24; end
if nargin<3, L=1.3; end, if nargin<4, turn_deg=0; end
g = 9.80665; dt = 1/fs; Tst = 0.55; Tsw = 0.45; Tcyc = Tst+Tsw;
lead = 3; N = round((lead + nstrides*Tcyc + 2)*fs);
t = (0:N-1)'*dt;
pos = zeros(N,3); pitch = zeros(N,1); yaw = zeros(N,1);
pcur = [0 0 0]; yawcur = 0; hd = 0;
for k = 1:N
    tk = t(k) - lead;
    if tk < 0 || tk >= nstrides*Tcyc, ph = -1; else, ph = mod(tk,Tcyc); end
    if ph >= Tst
        u = (ph-Tst)/Tsw;  m = 10*u^3-15*u^4+6*u^5;           % min-jerk 0..1
        stride_i = floor(tk/Tcyc);
        hd = deg2rad(turn_deg)*stride_i;
        pos(k,:) = pcur + L*m*[cos(hd) sin(hd) 0] + [0 0 0.12*sin(pi*u)];
        pitch(k) = deg2rad(35)*sin(2*pi*u);  yaw(k) = hd;
    else
        pos(k,:) = pcur; pitch(k)=0; yaw(k)=hd;
        if ph >= 0 && ph+dt >= Tst   % about to leave stance: commit new base position
            m1 = 1; pcur_next = pcur + L*[cos(hd+deg2rad(turn_deg)) sin(hd+deg2rad(turn_deg)) 0];
        end
    end
    % advance base position at end of each swing
    if ph >= Tst && ph+dt >= Tcyc
        hd2 = deg2rad(turn_deg)*floor(tk/Tcyc);
        pcur = pcur + L*[cos(hd2) sin(hd2) 0];
    end
end
v = gradient(pos',dt)'; a_w = gradient(v',dt)';
% orientation: world frame (X fwd0, Y left, Z up). Sensor: x=up, y=fwd, z=lateral.
% R_ws maps sensor->world: base (flat, yaw=0): sx->Z, sy->X, sz->Y
Rbase = [0 1 0; 0 0 1; 1 0 0];   % columns = sensor axes expressed in world
acc = zeros(N,3); gyr = zeros(N,3); Rs = zeros(3,3,N);
for k = 1:N
    Rz = [cos(yaw(k)) -sin(yaw(k)) 0; sin(yaw(k)) cos(yaw(k)) 0; 0 0 1];
    Ry = [cos(pitch(k)) 0 sin(pitch(k)); 0 1 0; -sin(pitch(k)) 0 cos(pitch(k))]; % pitch about world Y (lateral)
    Rs(:,:,k) = Rz*Ry*Rbase;
    acc(k,:) = (Rs(:,:,k)' * (a_w(k,:)' + [0;0;g]))';
end
for k = 2:N-1
    dR = Rs(:,:,k-1)' * Rs(:,:,k+1);               % body-frame rotation over 2 dt
    S = (dR - dR')/2;  gyr(k,:) = [S(3,2) S(1,3) S(2,1)]/(2*dt);
end
gyr(1,:) = gyr(2,:); gyr(N,:) = gyr(N-1,:);
rng(1);
acc = acc + 0.05*randn(N,3) + [0.08 -0.05 0.06];     % noise + constant bias
gyr = gyr + 0.01*randn(N,3) + [0.004 -0.003 0.002];
s.t=t; s.acc=acc; s.gyr=gyr; s.fs=fs; s.pos_true=pos; s.L=L; s.nstrides=nstrides;
end
