plotnonlin();

function Y = dxdt(t, X)
    m  = 0.127; % масса маятника
    M  = 1.206; % масса тележки
    I  = 0.001; % момент инерции маятника относительно центра масс
    l  = 0.178; % расстояние от точки крепления до центра масс
    Bc = 5.4;   % коэф. вязкого трения между кареткой и направляющей
    Bp = 0.002; % коэф. вязкого трения в точке крепления
    g  = 9.81;  % коэф. свободного падения

    theta = [0, 0, 0, 0]; % регулятор
    theta = [19.6330, -64.6908, 21.4010, -9.2987];

    x    = X(1);
    phi  = X(2);
    dx   = X(3);
    dphi = X(4);

    A = [
         m + M           ,-m * l * cos(phi); 
        -m * l * cos(phi), I + m * l^2
        ];

    F = theta * X;

    b = [
          F - Bc * dx - m * l * dphi^2 * sin(phi);
         -Bp * phi + m * g * l * sin(phi)
        ];

    Y = [dx; dphi; inv(A) * b];
end

function plotnonlin
    nx = 4; % количество строчек

    % задаем время отрисовки графиков
    TIME = 10.0;
    % задаем точки, в которых нужно вычислять решение системы
    ticks = 0 : 0.001 : TIME;
    % задаем параметры метода решения системы дифференциальных уравнений
    options = odeset('RelTol', 1e-5, 'AbsTol', 1e-5 * ones(1, nx));
    % задаем начальные условия 
    X0 = zeros(1, nx);   

    X0(1) = 0.5;
    X0(2) = -0.4;

    % находим решение замкнутой линейной системы
    % [TL, YL] = ode45(@dxdt, ticks, X0, options);

    t = 0;
    dt = 0.01;
    TL = [];
    YL = [];
    TL = [TL; t];
    YL = [YL; X0'];
    Y = X0';
    while (t < TIME)
        t = t + dt;
        Y = Y + dt * dxdt(t, Y);
        TL = [TL; t];
        YL = [YL, Y];
    end

    % --------------------------------------------------------------------
	% СТРОИМ ГРАФИКИ РЕШЕНИЙ
    % --------------------------------------------------------------------  

    fhandle = figure;
    subplot(2, 1, 1)
        plot(TL, YL(1, :), 'b', 'LineWidth', 2.0)
        grid on;
        xlabel('t', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('x(t)', 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('x_1^0 = %0.3f', X0(1)));
    subplot(2, 1, 2)
        plot(TL, YL(2, :), 'b', 'LineWidth', 2.0)
        grid on;
        xlabel('t', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('\xi(t)', 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('x_2^0 = %0.3f', X0(2)));
end