function out_csv = load_phone(acc_file, gyr_file, out_csv)
%LOAD_PHONE  Convert a phone-app export (separate accelerometer + gyroscope CSVs) to the
%            generic 7-column CSV  t,ax,ay,az,gx,gy,gz  that main_project/load_imu reads.
%   Works with Sensor Logger (Accelerometer.csv / Gyroscope.csv: time, seconds_elapsed, z, y, x)
%   and phyphox raw exports (Time (s), X, Y, Z). Columns are found by NAME, case-insensitive.
%   Units are auto-detected later: |a|~1 -> g, |a|~9.8 -> m/s^2; gyro > 50 -> deg/s else rad/s.
%   Example:  load_phone('Accelerometer.csv','Gyroscope.csv','my_walk.csv'); main_project('my_walk.csv','results_phone')
if nargin < 3, out_csv = 'phone_walk.csv'; end
[ta, A] = read_xyz(acc_file);  [tg, G] = read_xyz(gyr_file);
t0 = max(ta(1), tg(1));  t1 = min(ta(end), tg(end));
fs = round(1/median(diff(ta)));  fs = min(max(fs,50),200);       % target 50..200 Hz
t  = (t0 : 1/fs : t1)';
A  = interp1(ta, A, t, 'linear');  G = interp1(tg, G, t, 'linear');
fid = fopen(out_csv,'w'); fprintf(fid,'t,ax,ay,az,gx,gy,gz\n'); fclose(fid);
dlmwrite(out_csv, [t - t(1), A, G], '-append', 'precision', '%.6f');
fprintf('wrote %s : %.1f s @ %d Hz\n', out_csv, t(end)-t(1), fs);
end

function [t, XYZ] = read_xyz(f)
fid = fopen(f); hdr = fgetl(fid); fclose(fid);
h = strsplit(hdr, ',');  h = lower(strtrim(strrep(h, '"', '')));
M = dlmread(f, ',', 1, 0);
ct = find(strcmp(h,'seconds_elapsed'), 1);
if isempty(ct), ct = find(~cellfun(@isempty, regexp(h,'^time')), 1); end
if isempty(ct), ct = 1; end
t = M(:,ct);  if max(t) > 1e9, t = (t - t(1))/1e9; end            % ns epoch -> s
idx = zeros(1,3);  nm = {'x','y','z'};
for k = 1:3
    j = find(strcmp(h, nm{k}), 1);
    if isempty(j), j = find(~cellfun(@isempty, regexp(h, ['(^|\s)' nm{k} '(\s|\(|$)'])), 1); end
    if isempty(j), error('load_phone: column %s not found in %s (header: %s)', nm{k}, f, hdr); end
    idx(k) = j;
end
XYZ = M(:, idx);
[t, iu] = unique(t); XYZ = XYZ(iu,:);
end
