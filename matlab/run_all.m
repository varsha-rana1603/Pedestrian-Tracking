%RUN_ALL  Process every instructor file (CST*, WIPAD*) and write results/summary.csv
%   Edit DATA_DIR if your files live elsewhere.
if exist('OCTAVE_VERSION','builtin'), pkg load signal; end
DATA_DIR = '../data';   OUT = '../results';
files = [dir(fullfile(DATA_DIR,'CST*.csv')); dir(fullfile(DATA_DIR,'WIPAD*.CSV'))];
fid = fopen(fullfile(OUT,'summary.csv'),'w');
fprintf(fid,'file,fs,dur_s,stance_frac,naive_end_pos_m,naive_v_stance_rms,kf_max_pos_m,kf_v_pre_rms,fu_max_pos_m,fu_v_pre_rms,n_strides,stride_mean_m,stride_std_m,fu_path_m,fu_max_abs_z_m,fu_mean_abs_dz_m,stride_time_s,fu_end_to_start_m,fu_end_z_m\n');
for i = 1:numel(files)
    R = main_project(fullfile(DATA_DIR,files(i).name), OUT);
    fprintf(fid,'%s,%d,%.1f,%.3f,%.1f,%.2f,%.1f,%.3f,%.1f,%.3f,%d,%.3f,%.3f,%.1f,%.2f,%.3f,%.3f,%.2f,%.2f\n', R.file,R.fs,R.dur,R.stance_frac,R.naive_end_pos,R.naive_v_at_stance_rms,R.kf_max_pos,R.kf_v_pre_rms,R.fu_max_pos,R.fu_v_pre_rms,R.n_strides,R.fu_stride_mean,R.fu_stride_std,R.fu_path,R.fu_max_abs_z,R.fu_mean_abs_dz,R.stride_time_med,R.fu_end_to_start,R.fu_end_z);
end
fclose(fid);
