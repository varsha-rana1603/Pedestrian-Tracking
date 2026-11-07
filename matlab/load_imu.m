function d = load_imu(file, fs_out)
%LOAD_IMU  Read a CST*/WIPAD* (instructor) or phone CSV into one common struct.
%   d = load_imu(file)            -> uniform sampling (CST: 100 Hz, WIPAD: 200 Hz)
%   d = load_imu(file, fs_out)    -> uniform sampling at fs_out [Hz]
%
%   Output (SI units):  d.t [s], d.acc [m/s^2] (specific force, gravity included),
%                       d.gyr [rad/s], d.fs, d.fmt, d.info (data-quality notes)
%
%   CST   : t [s], ax..az [g], gx..gz [deg/s], mx..mz [uT], ClassifierState
%   WIPAD : currTime [ms] (wraps every 20000 ms!), accel [mg], gyro [deg/s]
if nargin < 2, fs_out = []; end
g = 9.80665;
fid = fopen(file); hdr = lower(strtrim(fgetl(fid))); fclose(fid);
M = dlmread(file, ',', 1, 0);
d = struct(); d.info = struct();

if ~isempty(strfind(hdr, 'currtime'))                  % ---------------- WIPAD
    d.fmt = 'WIPAD';
    t_ms = M(:,1);
    wraps = [0; cumsum(diff(t_ms) < -10000)];          % 20 s counter wrap
    t    = (t_ms + 20000*wraps) / 1000;                % stitched, continuous time
    acc  = M(:,2:4) / 1000 * g;                        % mg -> m/s^2
    gyr  = M(:,5:7) * pi/180;                          % deg/s -> rad/s
    d.info.n_blocks   = max(wraps) + 1;
    d.info.n_clipped  = sum(abs(M(:,2:4)) >= 3990, 1); % +-4 g accelerometer range
    d.info.eff_rate   = numel(t) / (t(end)-t(1));
    if isempty(fs_out), fs_out = 200; end
elseif ~isempty(strfind(hdr, 'classifierstate'))       % ------------------- CST
    d.fmt = 'CST';
    t   = M(:,2);
    acc = M(:,3:5) * g;                                % g -> m/s^2
    gyr = M(:,6:8) * pi/180;                           % deg/s -> rad/s
    d.mag   = M(:,9:11);  d.label = M(:,12);
    bad = all(M(:,3:8) == 0, 2);                       % e.g. all-zero first row of CST05
    t(bad) = NaN;
    d.info.n_bad_rows = sum(bad);
    d.info.eff_rate   = 1 / median(diff(t(~isnan(t))));
    if isempty(fs_out), fs_out = 100; end
else                                                   % ------ generic: t,ax,ay,az,gx,gy,gz
    d.fmt = 'GENERIC';
    t = M(:,1); acc = M(:,2:4); gyr = M(:,5:7);
    if median(sqrt(sum(acc.^2,2))) < 3, acc = acc * g; end   % looks like g
    if max(abs(gyr(:))) > 50, gyr = gyr * pi/180; end        % looks like deg/s
    if isempty(fs_out), fs_out = round(1/median(diff(t))); end
end

% ---- make time strictly increasing, then resample on a uniform grid
ok = ~isnan(t);  t = t(ok); acc = acc(ok,:); gyr = gyr(ok,:);
if isfield(d,'label'), d.label = d.label(ok); end
if isfield(d,'mag'),   d.mag   = d.mag(ok,:);  end
[t, iu] = unique(t, 'last');  acc = acc(iu,:); gyr = gyr(iu,:);
if isfield(d,'label'), d.label = d.label(iu); end
if isfield(d,'mag'),   d.mag   = d.mag(iu,:);  end
tu = (t(1) : 1/fs_out : t(end))';
d.acc = interp1(t, acc, tu, 'linear');
d.gyr = interp1(t, gyr, tu, 'linear');
if isfield(d,'label'), d.label = interp1(t, d.label, tu, 'nearest'); end
if isfield(d,'mag'),   d.mag   = interp1(t, d.mag,   tu, 'linear');  end
d.t  = tu - tu(1);
d.fs = fs_out;
d.name = file;
end
