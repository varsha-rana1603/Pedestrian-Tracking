PEDESTRIAN TRACKING

Setup Description: You are provided with an accelerometer (7Semi KXTJ3-1057 or 7Semi BMA400)
and a microcontroller (7Semi RP2040 Microcontroller).

Experimental Procedure: Attach the sensor system securely to your foot and perform a walk-
ing experiment along a predefined path. Record the accelerometer data using the microcontroller
and save the measurements in a CSV file for post-processing. Record a video of the walking ex-
periment to document the motion and support qualitative comparison with the estimated trajectory.

Objective:
(a) Pedestrian Tracking Without Kalman Filtering: Integrate the measured acceleration data
twice to estimate the foot’s position along the (x), (y), and (z) axes. Plot the estimated posi-
tion trajectories and compare them qualitatively with the observed walking motion. Does the
estimated trajectory accurately represent the motion? Discuss the sources of error and explain
why position estimates obtained through double integration may drift over time.

(b) Pedestrian Tracking Using Kalman Filtering: Apply a Kalman filter to the acceleration-based
position estimation pipeline developed in Part (a). Plot the estimated position trajectories and
compare the results with those obtained without Kalman filtering. Does filtering improve the
estimates? Discuss the observed improvements and limitations.

(c) Pedestrian Tracking Through Sensor Fusion: Since the provided sensor system contains only
an accelerometer and does not constitute a complete IMU, each team will receive a designated
IMU dataset from the instructor containing accelerometer and gyroscope measurements. Per-
form pedestrian tracking using two approaches:

Approach 1: Use only the accelerometer measurements.
Approach 2: Fuse the accelerometer and gyroscope measurements to estimate the pedestrian’s
motion.

Plot and compare the resulting trajectories. Identify the walking scenario represented by the
assigned dataset and discuss how sensor fusion influences the accuracy and reliability of pedes-
trian tracking.
