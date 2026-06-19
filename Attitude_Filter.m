% Attitude_Filter.m
% Complementary filter for roll/pitch (and yaw by integration) using IMU data.
% Data from MPU9250 (FLU frame) is converted to FRD (Forward-Right-Down).
% The filter matches the C code implementation provided earlier.

clear;
close all;

%% ------------------------------------------------------------------------
%  READ DATA FILES (adjust path if needed)
%  ------------------------------------------------------------------------
MPU9250 = readmatrix("Datasets\Dynamic\DS1\IMU.csv");

%% ------------------------------------------------------------------------
%  INITIALIZATIONS AND DEFINITIONS
%  ------------------------------------------------------------------------
imu_starting = 755;               % first IMU sample index
imu_obssec_s = MPU9250(imu_starting:end, 1);

% IMU data: accelerometer and gyro
% Original MPU9250 axes: X=Forward, Y=Left, Z=Up (FLU)
% We convert to FRD: X=Forward (same), Y=Right (negate Y), Z=Down (negate Z)
raw_f_ib_b = [ MPU9250(imu_starting:end, 2), ...          % X: forward
              -MPU9250(imu_starting:end, 3), ...          % Y: left -> right (negate)
              -MPU9250(imu_starting:end, 4)];             % Z: up   -> down  (negate)

raw_w_ib_b = [ deg2rad(MPU9250(imu_starting:end, 5)), ... % X gyro
              -deg2rad(MPU9250(imu_starting:end, 6)), ... % Y gyro (negate)
              -deg2rad(MPU9250(imu_starting:end, 7))];    % Z gyro (negate)

%% ------------------------------------------------------------------------
%  BIAS ESTIMATION (first 6000 samples assumed static)
%  ------------------------------------------------------------------------
avg_gyro = mean(raw_w_ib_b(1:6000, 1:3), 1);   % gyro bias (rad/s)
avg_accel = mean(raw_f_ib_b(1:6000, 1:3), 1);  % average acceleration (m/s^2)

%% ------------------------------------------------------------------------
%  INITIAL ATTITUDE (using the same formulas as the complementary filter)
%  ------------------------------------------------------------------------
% Normalize average acceleration
ax0 = avg_accel(1); ay0 = avg_accel(2); az0 = avg_accel(3);
acc_norm0 = sqrt(ax0^2 + ay0^2 + az0^2);
ax0 = ax0 / acc_norm0;
ay0 = ay0 / acc_norm0;
az0 = az0 / acc_norm0;

% Roll and pitch from accelerometer (in degrees)
roll_deg  = rad2deg( atan( ay0 / sqrt(ax0^2 + az0^2) ) );
pitch_deg = rad2deg( -atan( ax0 / sqrt(ay0^2 + az0^2) ) );

% Yaw – we set an initial value from your original code (156°)
% In a real application you might compute it from magnetometer or GNSS.
yaw_deg = 156.0;

%% ------------------------------------------------------------------------
%  FILTER PARAMETERS
%  ------------------------------------------------------------------------
COMP_FILTER_ALPHA = 0.995;   % trust gyro heavily

%% ------------------------------------------------------------------------
%  MAIN LOOP
%  ------------------------------------------------------------------------
loop_count = length(imu_obssec_s);

% Preallocate history (in degrees)
roll_hist  = zeros(loop_count, 1);
pitch_hist = zeros(loop_count, 1);
yaw_hist   = zeros(loop_count, 1);

for i = 2:loop_count
    % Time step (convert seconds to milliseconds)
    dt_ms = (imu_obssec_s(i) - imu_obssec_s(i-1)) * 1000;

    % Accelerometer readings (already in FRD)
    ax = raw_f_ib_b(i, 1) - 0.15;
    ay = raw_f_ib_b(i, 2) + 0.1;
    az = raw_f_ib_b(i, 3) + 0.45;

    % Gyro rates after bias removal, convert rad/s -> deg/s
    gx_deg = rad2deg(raw_w_ib_b(i, 1) - avg_gyro(1));
    gy_deg = rad2deg(raw_w_ib_b(i, 2) - avg_gyro(2));
    gz_deg = rad2deg(raw_w_ib_b(i, 3) - avg_gyro(3));

    % Update roll and pitch using complementary filter
    [roll_deg, pitch_deg] = computeAttitude(ax, ay, az, gx_deg, gy_deg, ...
                                            dt_ms, roll_deg, pitch_deg, ...
                                            COMP_FILTER_ALPHA);

    % Update yaw by pure integration (no correction)
    yaw_deg = yaw_deg + gz_deg * (dt_ms / 1000);

    % Store history
    roll_hist(i)  = roll_deg;
    pitch_hist(i) = pitch_deg;
    yaw_hist(i)   = yaw_deg;
end

%% ------------------------------------------------------------------------
%  PLOT RESULTS (decimated for clarity)
%  ------------------------------------------------------------------------
figure('Name', 'Attitude');
plot(roll_hist(1:100:end), 'LineWidth', 1.5); hold on;
plot(pitch_hist(1:100:end), 'LineWidth', 1.5);
% plot(yaw_hist(1:100:end), 'LineWidth', 1.5);
grid on;
xlabel('Sample index (decimated by 100)');
ylabel('Angle (degrees)');
legend('Roll', 'Pitch','Location', 'best');
title('Complementary Filter Attitude Estimates');

%% ------------------------------------------------------------------------
%  LOCAL FUNCTION: COMPLEMENTARY FILTER FOR ROLL/PITCH
%  ------------------------------------------------------------------------
function [phi, theta] = computeAttitude(ax, ay, az, gx_deg, gy_deg, dt_ms, phi, theta, alpha)
% COMPUTEATTITUDE  Complementary filter for roll and pitch (angles in degrees)
%   Inputs:
%     ax, ay, az   - accelerometer readings in m/s^2 (body frame)
%     gx_deg, gy_deg - gyro rates in deg/s (roll and pitch rates)
%     dt_ms        - time step in milliseconds
%     phi, theta   - previous roll and pitch angles in degrees
%     alpha        - filter weight (0<alpha<1), higher = trust gyro more
%   Outputs:
%     phi, theta   - updated angles in degrees

    % Convert time step to seconds
    dt = dt_ms / 1000.0;

    % Normalize accelerometer vector
    acc_norm = sqrt(ax^2 + ay^2 + az^2);
    if acc_norm == 0
        return;  % avoid division by zero
    end
    ax = ax / acc_norm;
    ay = ay / acc_norm;
    az = az / acc_norm;

    % Accelerometer-based roll and pitch (radians -> degrees)
    accel_roll_rad  = atan( ay / sqrt(ax^2 + az^2) );
    accel_pitch_rad = -atan( ax / sqrt(ay^2 + az^2) );
    accel_roll_deg  = rad2deg(accel_roll_rad);
    accel_pitch_deg = rad2deg(accel_pitch_rad);

    % Gyro integration (gyro rates are in deg/s, dt in s -> deg)
    gyro_roll_deg  = phi  + gx_deg * dt;
    gyro_pitch_deg = theta + gy_deg * dt;

    % Complementary filter blending
    phi   = alpha * gyro_roll_deg  + (1 - alpha) * accel_roll_deg;
    theta = alpha * gyro_pitch_deg + (1 - alpha) * accel_pitch_deg;
end