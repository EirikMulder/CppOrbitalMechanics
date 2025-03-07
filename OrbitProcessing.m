%% Orbital Display

orbit = csvread('orbit.csv', 1);

X = orbit(:, 1);
Y = orbit(:, 2);
Z = orbit(:, 3);

plot3(X, Y, Z, LineWidth=2)