% LCA_Simulation_Dynamic.m
% Loosely Coupled INS/GNSS Integration with 5 simulated GNSS outages (20 s each)
% Outage periods are highlighted in the plots and error statistics are computed.
% Corrected: outage error now means the drift at the very end of each outage.

clear; 
% close all;

%% ------------------------------------------------------------------------
%  READ DATA FILES (adjust paths as needed)
%  ------------------------------------------------------------------------
RLG_JZ105 = readmatrix("C:\Users\HP\Desktop\Abdullah Wasim\Datasets\Dynamic\DS1\IMU.csv");
PocketSDR_PVT = readmatrix("C:\Users\HP\Desktop\Abdullah Wasim\Datasets\Dynamic\DS1\GPS.csv");

%% ------------------------------------------------------------------------
%  INITIALIZATIONS AND DEFINITIONS
%  ------------------------------------------------------------------------
% For DS 1 Static
imu_starting = 755;          % first IMU sample index
neo_starting = 49;            % first GNSS sample index

% For DS 2 Static
% imu_starting = 1561;          % first IMU sample index
% neo_starting = 61;            % first GNSS sample index

micro_g_to_meters_per_second_squared = 9.80665E-6;

% Time vectors
imu_obssec_s     = RLG_JZ105(imu_starting:end, 1);
pocket_pvt_obssec_s = PocketSDR_PVT(neo_starting:end, 1);

% IMU data (accel, gyro)
raw_f_ib_b = RLG_JZ105(imu_starting:end, 2:4);
raw_w_ib_b = deg2rad(RLG_JZ105(imu_starting:end, 5:7));

% GNSS data (position, velocity)
pocket_lat = PocketSDR_PVT(neo_starting:end, 2);
pocket_lng = PocketSDR_PVT(neo_starting:end, 3);
pocket_h_b = PocketSDR_PVT(neo_starting:end, 4);
pocket_vn  = PocketSDR_PVT(neo_starting:end, 5);
pocket_ve  = PocketSDR_PVT(neo_starting:end, 6);

% max(pocket_vn)
% max(pocket_ve)

%% ------------------------------------------------------------------------
%  INITIAL ATTITUDE, POSITION, VELOCITY
%  ------------------------------------------------------------------------
roll  = 2.27962 * pi/180;
pitch = 1.43812 * pi/180;
yaw   = (156) * pi/180;

% Biases Calculations
avg_gyro = mean(raw_w_ib_b(1:6000,1:3), 1);

% Initialize Attitude
avg_accel = mean(raw_f_ib_b(1:6000,1:3), 1);
ax = avg_accel(1); ay = avg_accel(2); az = avg_accel(3);

norm_g = sqrt(ax^2 + ay^2 + az^2);
ax = ax / norm_g; ay = ay / norm_g; az = az / norm_g;

pitch = atan(ax / sqrt(ay^2 + az^2));
roll = atan2( -ay, -az);

wx = avg_gyro(1); wy = avg_gyro(2); wz = avg_gyro(3);

sin_psi = -wy*cos(roll) + wz*sin(roll);
cos_psi = wx*cos(pitch) + ...
    wy*sin(roll)*sin(pitch) + ...
    wz*cos(roll)*sin(pitch);
psi_nb = atan2(sin_psi, cos_psi);
% yaw = psi_nb;

init_C = Euler2DCM(roll, pitch, yaw);

init_lat  = deg2rad(mean(pocket_lat(1)));
init_lng  = deg2rad(mean(pocket_lng(1)));
init_h_b  = mean(pocket_h_b(1));

init_vn = mean(pocket_vn(1));
init_ve = mean(pocket_ve(1));
init_vd = 0;
init_v  = [init_vn; init_ve; init_vd];

[init_r_eb_e, init_v_eb_e, init_C_b_e] = NED2ECEF(init_lat, init_lng, init_h_b, init_v, init_C);

%% ------------------------------------------------------------------------
%  SIDE VARIABLES FOR THE MAIN LOOP
%  ------------------------------------------------------------------------
loop_count  = length(imu_obssec_s);
pocketCount = 61;          % current GNSS index
lastMatch   = -19;

%% ------------------------------------------------------------------------
%  DEFINE GNSS OUTAGES (5 intervals of 20 seconds each)
%  ------------------------------------------------------------------------
outage_intervals = [
    %  For DS 1 Static
    % 588000, 588020;
    % 588100, 588120;
    % 588200, 588220;
    % 588400, 588420;
    %  For DS 2 Static | Range (595292.204 to 597226.001)
    % 595800, 595820;
    % 596000, 596020;
    % 596700, 596720;
    % 596900, 596920;

    % For DS 1 Dynamic | Range (43389.300 to 44477.600)

    % 43630.00, 43650.00;
    % 43700.00, 43720.00;
    43800.00, 43820.00;
    43900.00, 43920.00;
    44000.00, 44020.00;
    44100.00, 44120.00;
    44230.00, 44250.00;
    44300.00, 44320.00;
    44400.00, 44420.00;
    
    ];

% Create logical mask: true if GNSS measurement is available (not in outage)
gnss_valid = true(length(pocket_pvt_obssec_s), 1);
for k = 1:size(outage_intervals, 1)
    in_interval = (pocket_pvt_obssec_s >= outage_intervals(k,1)) & ...
        (pocket_pvt_obssec_s <= outage_intervals(k,2));
    gnss_valid(in_interval) = false;
end

%% ------------------------------------------------------------------------
%  MAIN STATE VARIABLES
%  ------------------------------------------------------------------------
meas_w_ib_b = raw_w_ib_b(1,:)';
meas_f_ib_b = raw_f_ib_b(1,:)';

est_r_eb_e = init_r_eb_e;
est_v_eb_e = init_v_eb_e;
est_C_b_e  = init_C_b_e;

est_L_b     = init_lat;
est_lambda_b= init_lng;
est_h_b     = init_h_b;

output_profile = zeros(10, loop_count);   % [time, lat, lon, h, vn, ve, vd, roll, pitch, yaw]
output_bias = zeros(6, loop_count);   % [time, lat, lon, h, vn, ve, vd, roll, pitch, yaw]

kalman_performance = zeros(1, length(pocket_pvt_obssec_s));
performance_index = 0;

%% ------------------------------------------------------------------------
%  KALMAN FILTER CONFIGURATION
%  ------------------------------------------------------------------------
LC_KF_config.init_att_unc = deg2rad(1);
LC_KF_config.init_vel_unc = 0.1;
LC_KF_config.init_pos_unc = 10;
LC_KF_config.init_b_a_unc = 20 * micro_g_to_meters_per_second_squared;
LC_KF_config.init_b_g_unc = deg2rad(0.1) / 3600;

LC_KF_config.gyro_noise_PSD   = (deg2rad(0.02) / 60)^2;
LC_KF_config.accel_noise_PSD  = (200 * micro_g_to_meters_per_second_squared)^2;
LC_KF_config.accel_bias_PSD   = 1.0E-7;
LC_KF_config.gyro_bias_PSD    = 2.0E-12;

% LC_KF_config.init_att_unc = deg2rad(10);
% LC_KF_config.init_vel_unc = 0.5;
% LC_KF_config.init_pos_unc = 20;
% LC_KF_config.init_b_a_unc = 50000 * micro_g_to_meters_per_second_squared;
% LC_KF_config.init_b_g_unc = deg2rad(1) / 3600;
% 
% LC_KF_config.gyro_noise_PSD   = (deg2rad(0.5) / 60)^2;       % ~ 2.1e-8 (rad/s)^2/Hz
% LC_KF_config.accel_noise_PSD  = (0.1 * 9.80665e-3)^2;       % ~ 9.6e-7 (m/s^2)^2/Hz
% LC_KF_config.accel_bias_PSD   = 1.0e-5;                     % (m/s^2)^2/Hz
% LC_KF_config.gyro_bias_PSD    = 1.0e-9;                     % (rad/s)^2/Hz

% LC_KF_config.gyro_noise_PSD   = (deg2rad(0.02) / 60)^2 * 100;   % was *1
% LC_KF_config.accel_noise_PSD  = (200 * micro_g_to_meters_per_second_squared)^2 * 100;
% LC_KF_config.accel_bias_PSD   = 1.0E-7 * 100;
% LC_KF_config.gyro_bias_PSD    = 2.0E-12 * 100;

LC_KF_config.pos_meas_SD = 2.5;
LC_KF_config.vel_meas_SD = 0.1;

% LC_KF_config.pos_meas_SD = 0.5;   % was 1.5 m
% LC_KF_config.vel_meas_SD = 0.05;  % was 0.1 m/s

P_matrix = Initialize_P_Matrix(LC_KF_config);
est_IMU_bias = zeros(6,1);

% Lever Arm
l_ba_b = [0, 0, 0]';

%% ------------------------------------------------------------------------
%  MAIN LOOP (INS Mechanization + KF update with outage simulation)
%  ------------------------------------------------------------------------
for i = imu_starting:loop_count

    % Store old values
    old_w_ib_b = meas_w_ib_b;
    old_f_ib_b = meas_f_ib_b;
    old_r_eb_e = est_r_eb_e;
    old_v_eb_e = est_v_eb_e;
    old_C_b_e  = est_C_b_e;
    old_L_b    = est_L_b;
    old_lambda_b = est_lambda_b;
    old_h_b    = est_h_b;

    % INS mechanization
    dt_imu = imu_obssec_s(i) - imu_obssec_s(i-1);
    meas_w_ib_b = raw_w_ib_b(i,:)' - est_IMU_bias(4:6);
    meas_f_ib_b = raw_f_ib_b(i,:)' - est_IMU_bias(1:3);

    [est_r_eb_e, est_v_eb_e, est_C_b_e, P_matrix] = INS_Mechanization(dt_imu, old_r_eb_e, ...
        old_v_eb_e, old_C_b_e, P_matrix, meas_f_ib_b, meas_w_ib_b, est_L_b, LC_KF_config);
    % [est_r_eb_e, est_v_eb_e, est_C_b_e, ~] = INS_Mechanization(dt_imu, old_r_eb_e, ...
    %     old_v_eb_e, old_C_b_e, P_matrix, meas_f_ib_b, meas_w_ib_b, est_L_b, LC_KF_config);

    % Check if a GNSS measurement is due
    error = imu_obssec_s(i) - pocket_pvt_obssec_s(pocketCount);
    if (abs(error) < 0.0035) && ((i - lastMatch) >= 20) && (pocketCount < length(pocket_pvt_obssec_s))

        % fprintf("IMU: %0.3f | GPS: %0.3f\n", imu_obssec_s(i), pocket_pvt_obssec_s(pocketCount));

        if ~gnss_valid(pocketCount)
            % ----- GNSS OUTAGE: skip this measurement -----
            pocketCount = pocketCount + 1;   % move to next GNSS epoch
            % Do NOT update lastMatch -> next available GNSS will be used
        else
            % ----- Normal processing (valid GNSS) -----
            dt_gnss = pocket_pvt_obssec_s(pocketCount) - pocket_pvt_obssec_s(pocketCount-1);

            % Convert GNSS position/velocity to ECEF
            [GNSS_r_eb_e, GNSS_v_eb_e] = pv_NED2ECEF(deg2rad(pocket_lat(pocketCount)), ...
                deg2rad(pocket_lng(pocketCount)), pocket_h_b(pocketCount), ...
                [pocket_vn(pocketCount); pocket_ve(pocketCount); 0]);

            [est_C_b_e, est_v_eb_e, est_r_eb_e, est_IMU_bias, P_matrix, innov] = ...
                LC_KF_Epoch(GNSS_r_eb_e, GNSS_v_eb_e, dt_gnss, est_C_b_e, est_v_eb_e, ...
                est_r_eb_e, est_IMU_bias, P_matrix, meas_f_ib_b, meas_w_ib_b, est_L_b, LC_KF_config, l_ba_b);

            performance_index = performance_index + 1;

            kalman_performance(performance_index) = innov;

            % Advance to next GNSS epoch and record the match
            pocketCount = pocketCount + 1;
            lastMatch = i;
        end

    end

    % Convert state back to NED for logging
    [est_L_b, est_lambda_b, est_h_b, v_eb_n, C_b_n] = ECEF2NED(est_r_eb_e, est_v_eb_e, est_C_b_e);

    % Save output profile
    output_profile(1, i) = imu_obssec_s(i);
    output_profile(2:4, i) = [rad2deg(est_L_b), rad2deg(est_lambda_b), est_h_b];
    output_profile(5:7, i) = v_eb_n;
    output_profile(8:10, i) = rad2deg(DCM2Euler(C_b_n));

    output_bias(1:3, i) = est_IMU_bias(1:3);
    output_bias(4:6, i) = est_IMU_bias(4:6);
end

%% ------------------------------------------------------------------------
%  POST-PROCESSING: Compute errors at every GNSS epoch
%  ------------------------------------------------------------------------
num_gnss = length(pocket_pvt_obssec_s);
gnss_errors = zeros(num_gnss, 5);   % [time, horz_error, h_error, vn_error, ve_error]

for j = 1:num_gnss
    [~, idx_imu] = min(abs(output_profile(1,:) - pocket_pvt_obssec_s(j)));

    est_lat = output_profile(2, idx_imu);
    est_lon = output_profile(3, idx_imu);
    est_h   = output_profile(4, idx_imu);
    est_vn  = output_profile(5, idx_imu);
    est_ve  = output_profile(6, idx_imu);

    horz_error = haversine(est_lat, est_lon, pocket_lat(j), pocket_lng(j), 'm');
    h_error    = est_h - pocket_h_b(j);
    vn_error   = est_vn - pocket_vn(j);
    ve_error   = est_ve - pocket_ve(j);

    gnss_errors(j,:) = [pocket_pvt_obssec_s(j), horz_error, h_error, vn_error, ve_error];
end

% Logical mask for outage periods (based on GNSS time)
in_outage = false(num_gnss, 1);
for k = 1:size(outage_intervals, 1)
    in_outage = in_outage | (pocket_pvt_obssec_s >= outage_intervals(k,1) & ...
        pocket_pvt_obssec_s <= outage_intervals(k,2));
end

%% ------------------------------------------------------------------------
%  COMPUTE ERROR AT THE END OF EACH OUTAGE (drift)
%  ------------------------------------------------------------------------
num_outages = size(outage_intervals, 1);
end_outage_errors = zeros(num_outages, 5); % [time, horz_error, h_error, vn_error, ve_error]

for k = 1:num_outages
    % Indices of GNSS epochs inside this outage interval
    in_interval = (pocket_pvt_obssec_s >= outage_intervals(k,1)) & ...
        (pocket_pvt_obssec_s <= outage_intervals(k,2));
    idx_in_outage = find(in_interval);

    if ~isempty(idx_in_outage)
        % Take the last index (end of outage)
        last_idx = idx_in_outage(end);
        end_outage_errors(k, :) = gnss_errors(last_idx, :);
    else
        % No GNSS epoch inside this interval (should not happen for 20 s outages)
        end_outage_errors(k, :) = NaN;
    end
end

% Remove NaN rows if any
end_outage_errors = rmmissing(end_outage_errors);

%% ------------------------------------------------------------------------
%  STATISTICS
%  ------------------------------------------------------------------------
outage_err_all    = gnss_errors(in_outage, :);
nonoutage_err     = gnss_errors(~in_outage, :);

fprintf('\nNon‑outage periods (%d samples)\n', size(nonoutage_err,1));
fprintf('  Horizontal RMS : %.2f m\n', rms(nonoutage_err(:,2)));
fprintf('  Horizontal Max  : %.2f m\n', max(nonoutage_err(:,2)));
fprintf('  Height RMS      : %.2f m\n', rms(nonoutage_err(:,3)));
fprintf('  Vn RMS          : %.2f m/s\n', rms(nonoutage_err(:,4)));
fprintf('  Ve RMS          : %.2f m/s\n', rms(nonoutage_err(:,5)));

fprintf('\n=== End‑of‑Outage Drift Statistics (error at the very end of each outage) ===\n');
fprintf('  Horizontal RMS : %.2f m\n', rms(end_outage_errors(:,2)));
fprintf('  Horizontal Max  : %.2f m\n', max(end_outage_errors(:,2)));
fprintf('  Height RMS      : %.2f m\n', rms(end_outage_errors(:,3)));
fprintf('  Vn RMS          : %.2f m/s\n', rms(end_outage_errors(:,4)));
fprintf('  Ve RMS          : %.2f m/s\n', rms(end_outage_errors(:,5)));

%% ------------------------------------------------------------------------
%  PLOTTING
%  ------------------------------------------------------------------------

% 1. Horizontal error time series with outage shading
figure('Name', 'Horizontal Error');
plot(gnss_errors(:,1), gnss_errors(:,2), 'b-', 'LineWidth', 1); hold on;
% Mark end-of-outage drift points
plot(end_outage_errors(:,1), end_outage_errors(:,2), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
yl = ylim;
for k = 1:size(outage_intervals,1)
    x = [outage_intervals(k,1), outage_intervals(k,2), outage_intervals(k,2), outage_intervals(k,1)];
    y = [yl(1), yl(1), yl(2), yl(2)];
    patch(x, y, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
end
xlabel('GPS time (s)'); ylabel('Horizontal position error (m)');
title('Horizontal Error with Outage Periods Highlighted');
legend('All errors','End‑of‑outage drift','Location','best');
grid on;

% 2. Velocity errors with shading
figure('Name', 'Velocity Errors');
subplot(2,1,1);
plot(gnss_errors(:,1), gnss_errors(:,4), 'b-'); hold on;
yl = ylim;
for k = 1:size(outage_intervals,1)
    x = [outage_intervals(k,1), outage_intervals(k,2), outage_intervals(k,2), outage_intervals(k,1)];
    y = [yl(1), yl(1), yl(2), yl(2)];
    patch(x, y, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
end
ylabel('North velocity error (m/s)'); grid on;

subplot(2,1,2);
plot(gnss_errors(:,1), gnss_errors(:,5), 'b-'); hold on;
yl = ylim;
for k = 1:size(outage_intervals,1)
    x = [outage_intervals(k,1), outage_intervals(k,2), outage_intervals(k,2), outage_intervals(k,1)];
    y = [yl(1), yl(1), yl(2), yl(2)];
    patch(x, y, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
end
xlabel('GPS time (s)'); ylabel('East velocity error (m/s)'); grid on;

% 3. Trajectory (estimated line + valid GNSS reference dots)
figure('Name', 'Trajectory');
line(output_profile(2, 822:end), output_profile(3, 822:end), ...
     'LineStyle', 'none', 'Marker', '.', 'DisplayName', 'Estimated');
hold on;

% Plot only GNSS points that are NOT in an outage
valid_idx = find(gnss_valid);
line(pocket_lat(valid_idx), pocket_lng(valid_idx), ...
     'LineStyle', 'none', 'Marker', '.', 'DisplayName', 'GNSS (valid)');
hold off;

legend('Location', 'best');
xlabel('Latitude (deg)');
ylabel('Longitude (deg)');
title('Trajectory (Outage‑free GNSS shown)');
grid on;

% 4. Plot Biases
figure('Name',  'Accelerometer Biases')
plot(output_bias(1, 822:100:end));
hold on;
plot(output_bias(2, 822:100:end));
hold on;
plot(output_bias(3, 822:100:end));
hold off;
legend('X','Y','Z')

figure('Name',  'Gyroscope Biases')
plot(rad2deg(output_bias(4, 822:100:end)));
hold on;
plot(rad2deg(output_bias(5, 822:100:end)));
hold on;
plot(rad2deg(output_bias(6, 822:100:end)));
hold off;
legend('X','Y','Z')

figure('Name',  'Kalman Innovation');
plot(gnss_errors(100:end,1), kalman_performance(100:end));
% ylim([0 100]);

%% ------------------------------------------------------------------------
%  STORE 1 Hz ESTIMATED TRAJECTORY (regular time grid)
%  ------------------------------------------------------------------------
% Define a 1‑Hz time grid covering the whole dataset
t_start = ceil(min(output_profile(1, imu_starting:end)));
t_end   = floor(max(output_profile(1, imu_starting:end)));
t_1hz   = (t_start:t_end)';   % integer seconds

num_1hz = length(t_1hz);
output_1hz = zeros(num_1hz, 6);   % [time, est_lat, est_lon, est_h, est_vn, est_ve]

for j = 1:num_1hz
    [~, idx_imu] = min(abs(output_profile(1,:) - t_1hz(j)));
    output_1hz(j,1) = t_1hz(j);
    output_1hz(j,2) = output_profile(2, idx_imu);   % est_lat
    output_1hz(j,3) = output_profile(3, idx_imu);   % est_lon
    output_1hz(j,4) = output_profile(4, idx_imu);   % est_h
    output_1hz(j,5) = output_profile(5, idx_imu);   % est_vn
    output_1hz(j,6) = output_profile(6, idx_imu);   % est_ve
end

% Write raw matrix (no headers)
csv_filename = 'C:\Users\HP\Desktop\Abdullah Wasim\Datasets\estimated_trajectory_1Hz.csv';
writematrix(output_1hz, csv_filename);
fprintf('1 Hz estimated trajectory saved to %s (%d samples)\n', csv_filename, num_1hz);

% Write table with headers (more readable)
T = array2table(output_1hz, 'VariableNames', ...
    {'GPStime', 'lat_deg', 'lon_deg', 'h_m', 'vn_mps', 've_mps'});
writetable(T, 'C:\Users\HP\Desktop\Abdullah Wasim\Datasets\estimated_trajectory_1Hz_with_headers.csv');
fprintf('1 Hz estimated trajectory with headers saved.\n');

%% ------------------------------------------------------------------------
%  STORE GNSS REFERENCE DATA WITH OUTAGES REMOVED
%  ------------------------------------------------------------------------
% Only keep GNSS epochs that are NOT in the simulated outage intervals
valid_indices = find(gnss_valid);  % indices where gnss_valid is true
gnss_ref_valid = [pocket_pvt_obssec_s(valid_indices), ...
    pocket_lat(valid_indices), ...
    pocket_lng(valid_indices), ...
    pocket_h_b(valid_indices), ...
    pocket_vn(valid_indices), ...
    pocket_ve(valid_indices)];

% Write raw matrix
csv_ref_filename = 'C:\Users\HP\Desktop\Abdullah Wasim\Datasets\gnss_reference_no_outages.csv';
writematrix(gnss_ref_valid, csv_ref_filename);
fprintf('GNSS reference data (outages removed) saved to %s (%d samples)\n', ...
    csv_ref_filename, length(valid_indices));

% Write table with headers
T_ref = array2table(gnss_ref_valid, 'VariableNames', ...
    {'GPStime', 'lat_deg', 'lon_deg', 'h_m', 'vn_mps', 've_mps'});
writetable(T_ref, 'C:\Users\HP\Desktop\Abdullah Wasim\Datasets\gnss_reference_no_outages_with_headers.csv');
fprintf('GNSS reference with headers saved.\n');

%% ========================================================================
%  FUNCTIONS (unchanged from original)
%  ========================================================================

function [r_eb_e, v_eb_e, C_b_e] = NED2ECEF(L_b, lambda_b, h_b, v_eb_n, C_b_n)
R_0 = 6378137;
e = 0.0818191908425;
R_E = R_0 / sqrt(1 - (e * sin(L_b))^2);
cos_lat = cos(L_b); sin_lat = sin(L_b);
cos_long = cos(lambda_b); sin_long = sin(lambda_b);
r_eb_e = [(R_E + h_b) * cos_lat * cos_long;...
    (R_E + h_b) * cos_lat * sin_long;...
    ((1 - e^2) * R_E + h_b) * sin_lat];
C_e_n = [-sin_lat * cos_long, -sin_lat * sin_long,  cos_lat;...
    -sin_long,            cos_long,        0;...
    -cos_lat * cos_long, -cos_lat * sin_long, -sin_lat];
v_eb_e = C_e_n' * v_eb_n;
C_b_e = C_e_n' * C_b_n;
end

function [r_eb_e, v_eb_e] = pv_NED2ECEF(L_b, lambda_b, h_b, v_eb_n)
R_0 = 6378137;
e = 0.0818191908425;
R_E = R_0 / sqrt(1 - (e * sin(L_b))^2);
cos_lat = cos(L_b); sin_lat = sin(L_b);
cos_long = cos(lambda_b); sin_long = sin(lambda_b);
r_eb_e = [(R_E + h_b) * cos_lat * cos_long;...
    (R_E + h_b) * cos_lat * sin_long;...
    ((1 - e^2) * R_E + h_b) * sin_lat];
C_e_n = [-sin_lat * cos_long, -sin_lat * sin_long,  cos_lat;...
    -sin_long,            cos_long,        0;...
    -cos_lat * cos_long, -cos_lat * sin_long, -sin_lat];
v_eb_e = C_e_n' * v_eb_n;
end

function [L_b, lambda_b, h_b, v_eb_n, C_b_n] = ECEF2NED(r_eb_e, v_eb_e, C_b_e)
R_0 = 6378137;
e = 0.0818191908425;
lambda_b = atan2(r_eb_e(2), r_eb_e(1));
k1 = sqrt(1 - e^2) * abs(r_eb_e(3));
k2 = e^2 * R_0;
beta = sqrt(r_eb_e(1)^2 + r_eb_e(2)^2);
E = (k1 - k2) / beta;
F = (k1 + k2) / beta;
P = 4/3 * (E*F + 1);
Q = 2 * (E^2 - F^2);
D = P^3 + Q^2;
V = (sqrt(D) - Q)^(1/3) - (sqrt(D) + Q)^(1/3);
G = 0.5 * (sqrt(E^2 + V) + E);
T = sqrt(G^2 + (F - V * G) / (2 * G - E)) - G;
L_b = sign(r_eb_e(3)) * atan((1 - T^2) / (2 * T * sqrt(1 - e^2)));
h_b = (beta - R_0 * T) * cos(L_b) + ...
    (r_eb_e(3) - sign(r_eb_e(3)) * R_0 * sqrt(1 - e^2)) * sin(L_b);
cos_lat = cos(L_b); sin_lat = sin(L_b);
cos_long = cos(lambda_b); sin_long = sin(lambda_b);
C_e_n = [-sin_lat * cos_long, -sin_lat * sin_long,  cos_lat;...
    -sin_long,            cos_long,        0;...
    -cos_lat * cos_long, -cos_lat * sin_long, -sin_lat];
v_eb_n = C_e_n * v_eb_e;
C_b_n = C_e_n * C_b_e;
end

function [r_eb_e, v_eb_e, C_b_e, P_matrix_propagated] = INS_Mechanization(tor_i, old_r_eb_e, ...
    old_v_eb_e, old_C_b_e, P_matrix_old, f_ib_b, omega_ib_b, est_L_b_old, LC_KF_config)
omega_ie = 7.292115E-5;
alpha_ie = omega_ie * tor_i;
C_Earth = [cos(alpha_ie), sin(alpha_ie), 0;...
    -sin(alpha_ie), cos(alpha_ie), 0;...
    0,             0,             1];
alpha_ib_b = omega_ib_b * tor_i;
mag_alpha = sqrt(alpha_ib_b' * alpha_ib_b);
Alpha_ib_b = Skew_symmetric(alpha_ib_b);
if mag_alpha > 1.E-8
    C_new_old = eye(3) + sin(mag_alpha)/mag_alpha * Alpha_ib_b + ...
        (1 - cos(mag_alpha))/mag_alpha^2 * Alpha_ib_b * Alpha_ib_b;
    ave_C_b_e = old_C_b_e * (eye(3) + (1 - cos(mag_alpha))/mag_alpha^2 * Alpha_ib_b + ...
        (1 - sin(mag_alpha)/mag_alpha)/mag_alpha^2 * Alpha_ib_b * Alpha_ib_b) - ...
        0.5 * Skew_symmetric([0;0;alpha_ie]) * old_C_b_e * tor_i;
else
    C_new_old = eye(3) + (1-(mag_alpha^2)/6)*Alpha_ib_b + (0.5-(mag_alpha^2)/24)*Alpha_ib_b^2;
    ave_C_b_e = old_C_b_e * (eye(3) + (1-(mag_alpha^2)/6)*Alpha_ib_b + ...
        (0.5-(mag_alpha^2)/24)*Alpha_ib_b^2) - ...
        0.5 * Skew_symmetric([0;0;alpha_ie]) * old_C_b_e * tor_i;
end
C_b_e = C_Earth * old_C_b_e * C_new_old;
f_ib_e = ave_C_b_e * f_ib_b;
v_eb_e_pred = old_v_eb_e + tor_i * (f_ib_e + Gravity_ECEF(old_r_eb_e) - ...
    2 * Skew_symmetric([0;0;omega_ie]) * old_v_eb_e);
v_eb_e = old_v_eb_e + tor_i * (f_ib_e + Gravity_ECEF(old_r_eb_e) - ...
    0.5 * 2 * Skew_symmetric([0;0;omega_ie]) * (old_v_eb_e + v_eb_e_pred));
r_eb_e = old_r_eb_e + (v_eb_e + old_v_eb_e) * 0.5 * tor_i;

R_0 = 6378137;
e = 0.0818191908425;

Omega_ie = Skew_symmetric([0;0;omega_ie]);
% Transition matrix (first order)
Phi_matrix = eye(15);
Phi_matrix(1:3,1:3) = Phi_matrix(1:3,1:3) - Omega_ie * tor_i;
Phi_matrix(1:3,13:15) = C_b_e * tor_i;
Phi_matrix(4:6,1:3) = -tor_i * Skew_symmetric(C_b_e * f_ib_b);
Phi_matrix(4:6,4:6) = Phi_matrix(4:6,4:6) - 2 * Omega_ie * tor_i;
geocentric_radius = R_0 / sqrt(1 - (e * sin(est_L_b_old))^2) * ...
    sqrt(cos(est_L_b_old)^2 + (1 - e^2)^2 * sin(est_L_b_old)^2);
Phi_matrix(4:6,7:9) = -tor_i * 2 * Gravity_ECEF(r_eb_e) / geocentric_radius * ...
    r_eb_e' / sqrt(r_eb_e' * r_eb_e);
Phi_matrix(4:6,10:12) = C_b_e * tor_i;
Phi_matrix(7:9,4:6) = eye(3) * tor_i;
% System noise covariance
Q_prime_matrix = zeros(15);
Q_prime_matrix(1:3,1:3) = eye(3) * LC_KF_config.gyro_noise_PSD * tor_i;
Q_prime_matrix(4:6,4:6) = eye(3) * LC_KF_config.accel_noise_PSD * tor_i;
Q_prime_matrix(10:12,10:12) = eye(3) * LC_KF_config.accel_bias_PSD * tor_i;
Q_prime_matrix(13:15,13:15) = eye(3) * LC_KF_config.gyro_bias_PSD * tor_i;
% Propagate state (all zeros due to closed loop)\
P_matrix_propagated = Phi_matrix * (P_matrix_old + 0.5 * Q_prime_matrix) * ...
    Phi_matrix' + 0.5 * Q_prime_matrix;
end

function g = Gravity_ECEF(r_eb_e)
R_0 = 6378137;
mu = 3.986004418E14;
J_2 = 1.082627E-3;
omega_ie = 7.292115E-5;
mag_r = sqrt(r_eb_e' * r_eb_e);
if mag_r == 0
    g = [0;0;0];
else
    z_scale = 5 * (r_eb_e(3) / mag_r)^2;
    gamma = -mu / mag_r^3 * (r_eb_e + 1.5 * J_2 * (R_0 / mag_r)^2 * ...
        [(1 - z_scale) * r_eb_e(1); (1 - z_scale) * r_eb_e(2); (3 - z_scale) * r_eb_e(3)]);
    g = [gamma(1:2) + omega_ie^2 * r_eb_e(1:2); gamma(3)];
end
end

function [est_C_b_e_new, est_v_eb_e_new, est_r_eb_e_new, est_IMU_bias_new, P_matrix_new, innov] = ...
    LC_KF_Epoch(GNSS_r_eb_e, GNSS_v_eb_e, tor_s, est_C_b_e_old, est_v_eb_e_old, ...
    est_r_eb_e_old, est_IMU_bias_old, P_matrix_old, meas_f_ib_b, meas_omega_ib_b, est_L_b_old, LC_KF_config, l_ba_b)
c = 299792458;
omega_ie = 7.292115E-5;
R_0 = 6378137;
e = 0.0818191908425;
Omega_ie = Skew_symmetric([0;0;omega_ie]);
% Transition matrix (first order)
Phi_matrix = eye(15);
Phi_matrix(1:3,1:3) = Phi_matrix(1:3,1:3) - Omega_ie * tor_s;
Phi_matrix(1:3,13:15) = est_C_b_e_old * tor_s;
Phi_matrix(4:6,1:3) = -tor_s * Skew_symmetric(est_C_b_e_old * meas_f_ib_b);
Phi_matrix(4:6,4:6) = Phi_matrix(4:6,4:6) - 2 * Omega_ie * tor_s;
geocentric_radius = R_0 / sqrt(1 - (e * sin(est_L_b_old))^2) * ...
    sqrt(cos(est_L_b_old)^2 + (1 - e^2)^2 * sin(est_L_b_old)^2);
Phi_matrix(4:6,7:9) = -tor_s * 2 * Gravity_ECEF(est_r_eb_e_old) / geocentric_radius * ...
    est_r_eb_e_old' / sqrt(est_r_eb_e_old' * est_r_eb_e_old);
Phi_matrix(4:6,10:12) = est_C_b_e_old * tor_s;
Phi_matrix(7:9,4:6) = eye(3) * tor_s;
% System noise covariance
Q_prime_matrix = zeros(15);
Q_prime_matrix(1:3,1:3) = eye(3) * LC_KF_config.gyro_noise_PSD * tor_s;
Q_prime_matrix(4:6,4:6) = eye(3) * LC_KF_config.accel_noise_PSD * tor_s;
Q_prime_matrix(10:12,10:12) = eye(3) * LC_KF_config.accel_bias_PSD * tor_s;
Q_prime_matrix(13:15,13:15) = eye(3) * LC_KF_config.gyro_bias_PSD * tor_s;
% Propagate state (all zeros due to closed loop)
x_est_propagated = zeros(15,1);
P_matrix_propagated = Phi_matrix * (P_matrix_old + 0.5 * Q_prime_matrix) * ...
    Phi_matrix' + 0.5 * Q_prime_matrix;
% P_matrix_propagated = P_matrix_old;
% Measurement matrix
H_matrix = zeros(6,15);
H_matrix(1:3,7:9) = -eye(3);
H_matrix(4:6,4:6) = -eye(3);
% Measurement noise
R_matrix = zeros(6,6);
R_matrix(1:3,1:3) = eye(3) * LC_KF_config.pos_meas_SD^2;
R_matrix(4:6,4:6) = eye(3) * LC_KF_config.vel_meas_SD^2;

delta_z = [GNSS_r_eb_e - est_r_eb_e_old - est_C_b_e_old * l_ba_b; GNSS_v_eb_e - est_v_eb_e_old - est_C_b_e_old * cross(meas_omega_ib_b, l_ba_b) + Omega_ie * est_C_b_e_old * l_ba_b];

C = H_matrix * P_matrix_propagated * H_matrix' + R_matrix;
innov = delta_z' / C * delta_z;

% Kalman gain
K_matrix = P_matrix_propagated * H_matrix' * inv(H_matrix * P_matrix_propagated * H_matrix' + R_matrix);
% State update
x_est_new = x_est_propagated + K_matrix * delta_z;
P_matrix_new = (eye(15) - K_matrix * H_matrix) * P_matrix_propagated;
% Closed-loop correction
est_C_b_e_new = (eye(3) - Skew_symmetric(x_est_new(1:3))) * est_C_b_e_old;
est_v_eb_e_new = est_v_eb_e_old - x_est_new(4:6);
est_r_eb_e_new = est_r_eb_e_old - x_est_new(7:9);
est_IMU_bias_new = est_IMU_bias_old + x_est_new(10:15);
end

function P_matrix = Initialize_P_Matrix(LC_KF_config)
P_matrix = zeros(15);
P_matrix(1:3,1:3) = eye(3) * LC_KF_config.init_att_unc^2;
P_matrix(4:6,4:6) = eye(3) * LC_KF_config.init_vel_unc^2;
P_matrix(7:9,7:9) = eye(3) * LC_KF_config.init_pos_unc^2;
P_matrix(10:12,10:12) = eye(3) * LC_KF_config.init_b_a_unc^2;
P_matrix(13:15,13:15) = eye(3) * LC_KF_config.init_b_g_unc^2;
end

function A = Skew_symmetric(a)
A = [0, -a(3), a(2); a(3), 0, -a(1); -a(2), a(1), 0];
end

function distance = haversine(lat1, lon1, lat2, lon2, varargin)
if nargin < 5
    unit = 'km';
else
    unit = lower(varargin{1});
end
R = 6371000; % meters
lat1 = deg2rad(lat1); lon1 = deg2rad(lon1);
lat2 = deg2rad(lat2); lon2 = deg2rad(lon2);
dlat = lat2 - lat1;
dlon = lon2 - lon1;
a = sin(dlat/2).^2 + cos(lat1) .* cos(lat2) .* sin(dlon/2).^2;
c = 2 * atan2(sqrt(a), sqrt(1-a));
distance_m = R * c;
switch unit
    case {'km','kilometers'}, distance = distance_m / 1000;
    case {'miles','mi'},      distance = distance_m / 1609.344;
    case {'m','meters'},      distance = distance_m;
    case {'nm','nauticalmiles','nmi'}, distance = distance_m / 1852;
    otherwise,                 distance = distance_m / 1000;
end
end

function C = Euler2DCM(phi, theta, psi)
sin_phi = sin(phi); cos_phi = cos(phi);
sin_theta = sin(theta); cos_theta = cos(theta);
sin_psi = sin(psi); cos_psi = cos(psi);
C = [cos_theta*cos_psi, cos_theta*sin_psi, -sin_theta;
    -cos_phi*sin_psi + sin_phi*sin_theta*cos_psi, cos_phi*cos_psi + sin_phi*sin_theta*sin_psi, sin_phi*cos_theta;
    sin_phi*sin_psi + cos_phi*sin_theta*cos_psi, -sin_phi*cos_psi + cos_phi*sin_theta*sin_psi, cos_phi*cos_theta];
end

function eul = DCM2Euler(C)
eul = [atan2(C(2,3), C(3,3)); -asin(C(1,3)); atan2(C(1,2), C(1,1))];
end