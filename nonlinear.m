nonlinear_plot();

function A = give_matrix(phi)
    m  = 0.127;
    M  = 1.206;
    L  = 0.178;
    I  = 0.001;

    A = [m + M, -m * L * cos(phi); -m * L * cos(phi), I + m * L^2];
end

function F = give_right_part(phi, dphi)
    m  = 0.127;
    g  = 9.81;
    L  = 0.178;
    F = [-m * L * dphi * dphi * sin(phi); m * g * L * sin(phi)];
end

function f = rp(t, X)
    th1 =  19.6330;
    th2 = -64.6908;
    th3 =  21.4010;
    th4 =  -9.2987;

    Bc = 5.4;
    Bp = 0.002;

    B = [Bc - th3, -th4; 0, Bp];
    C = [-th1, -th2; 0, 0];
    A = give_matrix(X(2));
    F = give_right_part(X(2), X(4));

    bigA = [zeros(2, 2), eye(2, 2); -inv(A) * C, -inv(A) * B];
    bigF = [0; 0; F];

    f = bigA * X + bigF;
end

function nonlinear_plot
    % задаем время отрисовки графиков
    TIME = 20.0;
    % задаем точки, в которых нужно вычислять решение системы
    ticks = 0 : 0.001 : TIME;
    % задаем параметры метода решения системы дифференциальных уравнений
    options = odeset('RelTol', 1e-5, 'AbsTol', 1e-5 * ones(1, 4));
    % задаем начальные условия 
    X0 = zeros(1, 4);   

    X0(1) = 0.9;
    X0(2) = -0.9;

    % находим решение замкнутой линейной системы
    [TL, YL] = ode45(@(t, X)rp(ticks, X), ticks, X0, options);
    
    fhandle = figure;
    subplot(2, 2, 1)
        plot(TL, YL(:, 1), 'b', 'LineWidth', 2.0)
        grid on;
        xlabel('t', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('x(t)', 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('x_1^0 = %0.3f', X0(1)));
    subplot(2, 2, 2)
        plot(TL, YL(:, 2), 'b', 'LineWidth', 2.0)
        grid on;
        xlabel('t', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('x(t)', 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('x_1^0 = %0.3f', X0(2)));
    subplot(2, 2, 3)
        plot(TL, YL(:, 3), 'b', 'LineWidth', 2.0)
        grid on;
        xlabel('t', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('x(t)', 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('x_1^0 = %0.3f', X0(3)));
    subplot(2, 2, 4)
        plot(TL, YL(:, 4), 'b', 'LineWidth', 2.0)
        grid on;
        xlabel('t', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('x(t)', 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('x_1^0 = %0.3f', X0(4)));
end