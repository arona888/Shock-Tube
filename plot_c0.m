% Constants
gamma = 1.31;
R = 8.3145;
M = 16.04e-3;

% Temperature range (0C to 100C, in Kelvin)
T = linspace(200, 300, 1000);

% Compute speed of sound
v = sqrt(gamma .* R .* T ./ M);

% Plot
figure;
plot(T - 273, v);
xlabel('Temperature (C)');
ylabel('Speed of Sound (m/s)');
title('Speed of Sound in Methane vs. Temperature');
grid on;