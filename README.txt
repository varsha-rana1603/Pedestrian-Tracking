Problem 1 - Part (c): accelerometer-only vs accelerometer+gyro tracking
Run (MATLAB or Octave), from the matlab/ folder:
  1. Put your IMU file(s) in ../data/   (CST*.csv or WIPAD*.CSV)
  2. main_project('../data/WIPAD09.CSV', '../results')     % one file
     run_all                                                % every file + summary.csv
Files: load_imu (loader), detect_stance (SHOE), find_init_window, still_mask,
ins_naive (a), kf_zupt (b / Approach 1), ins_eskf (Approach 2), analyze_strides,
main_project (figures+metrics), sim_walk/test_synthetic (validation), load_phone.
