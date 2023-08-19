% Dummy distance values from 0 to 1000 m
%distance = linspace(-500, 0, 1000);

offset = 300;

idx = N/2-offset:N/2;

distance = x(idx) - x(idx(end));

% Generate dummy values for each quantity
%pressure = rand(1, 1000)*1e5;    % Pressure values in Pa
%energy = rand(1, 1000)*1e6;      % Energy values in Joules
%velocity = rand(1, 1000)*100;    % Velocity values in m/s
%temperature = rand(1, 1000)*300; % Temperature values in K


% pressure = pressure_plot(idx,[1,2,4,6,8]);
% energy = e_plot(idx,[1,2,4,6,8]);
% rho = rho_plot(idx,[1,2,4,6,8]);
% velocity = u_plot(idx,[1,2,4,6,8]);
% temperature = T_plot(idx,[1,2,4,6,8]);
% 
% time_labels = ["Time = 0s", "Time = 1s", "Time = 3s", "Time = 5s", "Time = 7s"];

% pressure = [pressure_plot(idx,end), pressure_plot_save(idx,end)];
% energy = [e_plot(idx,end), e_plot_save(idx,end)];
% rho = [rho_plot(idx,end), rho_plot_save(idx,end)];
% velocity = [u_plot(idx,end), u_plot_save(idx,end)];
% temperature = [T_plot(idx,end), T_plot_save(idx,end)];
% 
% time_labels = ["No friction", "Friction"];

pressure = [pressure_plot(idx,8), pressure_plot(idx,end), pressure_plot_save(idx,8), pressure_plot_save(idx,end)];
energy = [e_plot(idx,8), e_plot(idx,end), e_plot_save(idx,8), e_plot_save(idx,end)];
rho = [rho_plot(idx,8), rho_plot(idx,end), rho_plot_save(idx,8), rho_plot_save(idx,end)];
velocity = [u_plot(idx,8), u_plot(idx,end), u_plot_save(idx,8), u_plot_save(idx,end)];
temperature = [T_plot(idx,8), T_plot(idx,end), T_plot_save(idx,8), T_plot_save(idx,end)];

time_labels = ["No friction t=7", "No friction t=60", "Friction t=7", "Friction t=60"];


% Plot

% Pressure vs distance
subplot(2, 2, 1);
plot(distance, pressure);
title('Pressure');
xlabel('Position (m)');
ylabel('Pressure (Pa)');
legend(time_labels, 'Location', 'best');
grid on;

% Energy vs distance
% subplot(2, 2, 2);
% plot(distance, energy);
% title('Energy');
% xlabel('Position (m)');
% ylabel('Energy (J)');
% legend(time_labels, 'Location', 'best');

subplot(2, 2, 2);
plot(distance, rho);
title('Density');
xlabel('Position (m)');
ylabel('Density (kg/m3)');
legend(time_labels, 'Location', 'best');

grid on;

% Velocity vs distance
subplot(2, 2, 3);
plot(distance, velocity);
title('Velocity');
xlabel('Position (m)');
ylabel('Velocity (m/s)');
legend(time_labels, 'Location', 'best');

grid on;

% Temperature vs distance
subplot(2, 2, 4);
plot(distance, temperature);
title('Temperature');
xlabel('Position (m)');
ylabel('Temperature (K)');
legend(time_labels, 'Location', 'best');

grid on;

% Create a new set of axes that cover the entire figure
%hLegend = axes('visible', 'off');

% Create the legend in the new axes
%legend(hLegend, time_labels, 'Location', 'northoutside', 'Orientation', 'horizontal');