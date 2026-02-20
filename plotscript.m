%% tank 1, no disturbance
time = out.h1sim(:, 1);
h1meas = load('h1meas.mat', 'h1'); % load measurement data
h1meas = h1meas.h1; % rename

h1sim = out.h1sim(:, 2); % load simulated data

h1stat = out.h1stat(:, 2); % load stationary kalman filter estimation

h1nstat = out.h1nstat(:, 2); % load non-stationary kalman filter estimation

h1ext = out.h1ext(:, 2); % load extended kalman filter estimation

%%%%%%%%%%%%% plots for tank 1 %%%%%%%%%%%%%%%%
figure;
plot(time, h1meas, 'k--', 'LineWidth', 2)
hold on
plot(time, h1sim, 'r-', 'LineWidth', 2)
plot(time, h1stat, 'm-', 'LineWidth', 2)
plot(time, h1nstat, 'b-', 'LineWidth', 2)
plot(time, h1ext, 'g-', 'LineWidth', 2)
grid on
legend('Measurement data', 'Simulated data', 'SKF', 'nSKF', 'EKF', 'Location', 'SouthEast')
title('Measured, simulated and estimated water level in upper tank')
xlabel('Time (s)')
ylabel('Water level (cm)')

%%%%%%%%%%%%%% error plots tank 1 %%%%%%%%%%%%%%
error_h1sim = h1meas - h1sim;
error_h1stat = h1meas - h1stat;
error_h1nstat = h1meas - h1nstat;
error_h1ext = h1meas - h1ext;

figure;
hold on
plot(time, error_h1sim, 'r-', 'LineWidth', 2)
plot(time, error_h1stat, 'm-', 'LineWidth', 2)
plot(time, error_h1nstat, 'b--', 'LineWidth', 2)
plot(time, error_h1ext, 'g--', 'LineWidth', 2)
grid on
legend('Simulated data', 'SKF', 'nSKF', 'EKF', 'Location', 'SouthEast')
title('Error plot upper tank')
xlabel('Time (s)')
ylabel('Model offset (cm)')

% rms error
rms_h1sim = rmse(h1meas, h1sim)
rms_h1stat = rmse(h1meas, h1stat)
rms_h1nstat = rmse(h1meas, h1nstat)
rms_h1ext = rmse(h1meas, h1ext)

%% tank 2, no disturbance

h2meas = load('h2meas.mat', 'h2'); % measurement
h2meas = h2meas.h2; % rename

h2sim = out.h2sim(:, 2); % simulation data

h2stat = out.h2stat(:, 2); % SKF

h2nstat = out.h2nstat(:, 2); %nSKF

h2ext = out.h2ext(:, 2); % EKF

%%%%%%%%%%%%% plots for tank 2 %%%%%%%%%%%%%%%%
figure;
plot(time, h2meas, 'k--', 'LineWidth', 2)
hold on
plot(time, h2sim, 'r-', 'LineWidth', 2)
plot(time, h2stat, 'm-', 'LineWidth', 2)
plot(time, h2nstat, 'b--', 'LineWidth', 2)
plot(time, h2ext, 'g--', 'LineWidth', 2)
grid on
legend('Measurement data', 'Simulated data', 'SKF', 'nSKF', 'EKF', 'Location', 'SouthEast')
title('Measured, simulated and estimated water level in lower tank')
xlim([0, 600])
ylim([-3, 15])
xlabel('Time (s)')
ylabel('Water level (cm)')

%%%%%%%%%%%%%% error plots tank 1 %%%%%%%%%%%%%%
error_h2sim = h2meas - h2sim;
error_h2stat = h2meas - h2stat;
error_h2nstat_d = h2meas - h2nstat;
error_h2ext = h2meas - h2ext;

figure;
hold on
plot(time, error_h2sim, 'r-', 'LineWidth', 2)
plot(time, error_h2stat, 'm-', 'LineWidth', 2)
plot(time, error_h2nstat_d, 'b--', 'LineWidth', 2)
plot(time, error_h2ext, 'g--', 'LineWidth', 2)
grid on
legend('Simulated data', 'SKF', 'nSKF', 'EKF', 'Location', 'SouthEast')
title('Error plot lower tank')
xlim([0, 600])
ylim([-2, 2])
xlabel('Time (s)')
ylabel('Model offset (cm)')

% rms error
rms_h2sim = rmse(h2meas, h2sim)
rms_h2stat = rmse(h2meas, h2stat)
rms_h2nstat = rmse(h2meas, h2nstat)
rms_h2ext = rmse(h2meas, h2ext)

%% tank 1, disturbance

h1meas_d = load('h1meas_d.mat', 'h1');
h1meas_d = h1meas_d.h1; % rename
time = out.h1sim_d(:, 1);

h1sim_d = out.h1sim_d(:, 2); % load simulated data

h1stat_d = out.h1stat_d(:, 2); % load stationary kalman filter estimation

h1nstat_d = out.h1nstat_d(:, 2); % load non-stationary kalman filter estimation

h1ext_d = out.h1ext_d(:, 2); % load extended kalman filter estimation

%%%%%%%%%%%%% plots for tank 1 %%%%%%%%%%%%%%%%
figure;
plot(time, h1meas_d, 'k--', 'LineWidth', 2)
hold on
plot(time, h1sim_d, 'r-', 'LineWidth', 2)
plot(time, h1stat_d, 'm-', 'LineWidth', 2)
plot(time, h1nstat_d, 'b-', 'LineWidth', 2)
plot(time, h1ext_d, 'g-', 'LineWidth', 2)
grid on
legend('Measurement data', 'Simulated data', 'SKF', 'nSKF', 'EKF', 'Location', 'SouthEast')
title('Measured, simulated and estimated water level in upper tank')
xlabel('Time (s)')
ylabel('Water level (cm)')

%%%%%%%%%%%%%% error plots tank 1 %%%%%%%%%%%%%%
error_h1sim_d = h1meas_d - h1sim_d;
error_h1stat_d = h1meas_d - h1stat_d;
error_h1nstat_d = h1meas_d - h1nstat_d;
error_h1ext_d = h1meas_d - h1ext_d;

figure;
hold on
plot(time, error_h1sim_d, 'r-', 'LineWidth', 2)
plot(time, error_h1stat_d, 'm-', 'LineWidth', 2)
plot(time, error_h1nstat_d, 'b--', 'LineWidth', 2)
plot(time, error_h1ext_d, 'g--', 'LineWidth', 2)
grid on
legend('Simulated data', 'SKF', 'nSKF', 'EKF', 'Location', 'SouthEast')
title('Error plot upper tank')
xlabel('Time (s)')
ylabel('Model offset (cm)')

% rms error
rms_h1sim_d = rmse(h1meas_d, h1sim_d)
rms_h1stat_d = rmse(h1meas_d, h1stat_d)
rms_h1nstat_d = rmse(h1meas_d, h1nstat_d)
rms_h1ext_d = rmse(h1meas_d, h1ext_d)


%% tank 2, disturbance

h2meas_d = load('h2meas_d.mat', 'h2');
h2meas_d = h2meas_d.h2; % rename

h2sim_d = out.h2sim_d(:, 2); % simulation data

h2stat_d = out.h2stat_d(:, 2); % SKF

h2nstat_d = out.h2nstat_d(:, 2); %nSKF

h2ext_d = out.h2ext_d(:, 2); % EKF

%%%%%%%%%%%%% plots for tank 2 %%%%%%%%%%%%%%%%
figure;
plot(time, h2meas_d, 'k--', 'LineWidth', 2)
hold on
plot(time, h2sim_d, 'r-', 'LineWidth', 2)
plot(time, h2stat_d, 'm-', 'LineWidth', 2)
plot(time, h2nstat_d, 'b--', 'LineWidth', 2)
plot(time, h2ext_d, 'g--', 'LineWidth', 2)
grid on
legend('Measurement data', 'Simulated data', 'SKF', 'nSKF', 'EKF', 'Location', 'SouthEast')
title('Measured, simulated and estimated water level in lower tank')
xlim([0, 600])
ylim([-3, 15])
xlabel('Time (s)')
ylabel('Water level (cm)')

%%%%%%%%%%%%%% error plots tank 2 %%%%%%%%%%%%%%
error_h2sim_d = h2meas_d - h2sim_d;
error_h2stat_d = h2meas_d - h2stat_d;
error_h2nstat_d = h2meas_d - h2nstat_d;
error_h2ext_d = h2meas_d - h2ext_d;

figure;
hold on
plot(time, error_h2sim_d, 'r-', 'LineWidth', 2)
plot(time, error_h2stat_d, 'm-', 'LineWidth', 2)
plot(time, error_h2nstat_d, 'b--', 'LineWidth', 2)
plot(time, error_h2ext_d, 'g--', 'LineWidth', 2)
grid on
legend('Simulated data', 'SKF', 'nSKF', 'EKF', 'Location', 'SouthEast')
title('Error plot lower tank')
xlim([0, 600])
ylim([-2, 2])
xlabel('Time (s)')
ylabel('Model offset (cm)')

% rms error
rms_h2sim_d = rmse(h2meas_d, h2sim_d)
rms_h2stat_d = rmse(h2meas_d, h2stat_d)
rms_h2nstat_d = rmse(h2meas_d, h2nstat_d)
rms_h2ext_d = rmse(h2meas_d, h2ext_d)

%% Valve opening

valve_d = out.valve_d(:, 2); % simulation data

valvestat_d = out.valvestat_d(:, 2); % SKF

valvenstat_d = out.valvenstat_d(:, 2); %nSKF

valveext_d = out.valveext_d(:, 2); % EKF

%%%%%%%%%%%%% plots for valve %%%%%%%%%%%%%%%%
figure;
hold on
plot(time, valve_d, 'r-', 'LineWidth', 2)
plot(time, valvestat_d, 'm-', 'LineWidth', 2)
plot(time, valvenstat_d, 'g-', 'LineWidth', 2)
plot(time, valveext_d, 'b--', 'LineWidth', 2)
grid on
legend('Valve opening', 'SKF', 'nSKF', 'EKF', 'Location', 'SouthEast')
title('Actuated and estimated valve opening')
xlim([0, 600])
ylim([-0.2, 0.2])
xlabel('Time (s)')
ylabel('Valve opening area (cm^2)')

%%%%%%%%%%%%%% error plots valve %%%%%%%%%%%%%%

error_valvestat_d = valve_d - valvestat_d;
error_valvenstat_d = valve_d - valvenstat_d;
error_valveext_d = valve_d - valveext_d;

figure;
hold on
plot(time, error_valvestat_d, 'm-', 'LineWidth', 2)
plot(time, error_valvenstat_d, 'g-', 'LineWidth', 2)
plot(time, error_valveext_d, 'b--', 'LineWidth', 2)
grid on
legend('SKF', 'nSKF', 'EKF', 'Location', 'SouthEast')
title('Error plot valve opening')
xlim([0, 600])
ylim([-0.2, 0.2])
xlabel('Time (s)')
ylabel('Model offset (cm)')

% rms error
rms_valvestat_d = rmse(valve_d, valvestat_d)
rms_valvenstat_d = rmse(valve_d, valvenstat_d)
rms_valveext_d = rmse(valve_d, valveext_d)

%% Reference voltage

t = out.vref(:, 1);
voltage = out.vref(:, 2);

figure;
hold on
grid on
plot(t, voltage, 'k-', 'LineWidth', 2)
legend('Reference signal', 'Location', 'SouthEast')
title('Reference voltage')
xlim([0, 600])
ylim([-0.1, 5])
xlabel('Time (s)')
ylabel('Voltage (V)')