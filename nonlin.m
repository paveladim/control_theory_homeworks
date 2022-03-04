clc;
clear;
g  = 9.81;  % коэф. свободного падения

plot_nonlin();

function A = retA(phi)
    m  = 0.127;
    M  = 1.206; 
    I  = 0.001; 
    l  = 0.178; 
    A = [m + M, -m * l * cos(phi); -m * l * cos(phi), I + m * l^2];
end

function B = retB()
    Bc = 5.4;
    Bp = 0.002;
    B = diag([Bc, Bp]);
end

function C = retC(phi, dphi)
    m  = 0.127;
    l  = 0.178;
    g = 9.81;

    C = [m * l * dphi^2 * sin(phi); -m * g * l * sin(phi)];
end

function F = retF(X)
    theta = [0, 0, 0, 0];
    F = theta * X';
end

function ddX = retddX(t, X)
    A = retA(X(2));
    B = retB();
    C = retC(X(2), X(4));
    F = retF(X);

    bigA = [zeros(2, 2), diag([1, 1]); zeros(2, 2), -inv(A) * B];
    tempA = inv(A) * [1; 0];
    bigF = [0; 0; -inv(A) * C + tempA * F];
    ddX = (bigA * X' + bigF)';
end

function plot_nonlin()
    % задаем время отрисовки графиков
    TIME = 10.0;
    % задаем точки, в которых нужно вычислять решение системы
    ticks = 0 : 0.001 : TIME;
    % задаем параметры метода решения системы дифференциальных уравнений
    options = odeset('RelTol', 1e-5, 'AbsTol', 1e-5 * ones(1, 4));
    % задаем начальные условия 
    X0 = zeros(1, 4);

    X0(1) = 0.1;
    X0(2) = -0.1;
    
    retddX(0.1, X0);
    % находим решение замкнутой линейной системы
    [TL, YL] = ode45(@(t, X) retddX(t, X), ticks, X0, options);

    fhandle = figure;
    subplot(2, 1, 1)
        plot(TL, YL(:, 1), 'b', TL, YL(:, 3), 'r', 'LineWidth', 2.0)
        grid on;
        xlabel('t', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('x(t)', 'FontSize', 12, 'FontWeight', 'bold');
    subplot(2, 1, 2)
        plot(TL, YL(:, 2), 'b', TL, YL(:, 4), 'r', 'LineWidth', 2.0)
        grid on;
        xlabel('t', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('\xi(t)', 'FontSize', 12, 'FontWeight', 'bold');
end