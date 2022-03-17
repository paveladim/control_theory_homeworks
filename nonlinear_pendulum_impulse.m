nonlin_IMP();

function Y = dxdt(t, X, theta)
    m  = 0.127; % масса маятника
    M  = 1.206; % масса тележки
    I  = 0.001; % момент инерции маятника относительно центра масс
    l  = 0.178; % расстояние от точки крепления до центра масс
    Bc = 5.4;   % коэф. вязкого трения между кареткой и направляющей
    Bp = 0.002; % коэф. вязкого трения в точке крепления
    g  = 9.81;  % коэф. свободного падения

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
         -Bp * dphi + m * g * l * sin(phi)
        ];

    Y = [dx; dphi; inv(A) * b];
end

function nonlin_IMP
    m  = 0.127; % масса маятника
    M  = 1.206; % масса тележки
    I  = 0.001; % момент инерции маятника относительно центра масс
    l  = 0.178; % расстояние от точки крепления до центра масс
    Bc = 5.4;   % коэф. вязкого трения между кареткой и направляющей
    Bp = 0.002; % коэф. вязкого трения в точке крепления
    g  = 9.81;  % коэф. свободного падения

    % определим матрицы системы
    A0 = [m + M, -m * l; 
          -m * l, I + m * l^2];

    A1 = diag([Bc, Bp]);
    A2 = diag([0, -m * g * l]);
    B = [0; 0; inv(A0) * [1; 0]];

    % основная матрица системы
    A = [zeros(2, 2), eye(2); -inv(A0) * A2, -inv(A0) * A1];

    h = 0.1;
    Ad = expm(A * h);

    Bd = expm(A * h) * B;
    
    theta = -place(Ad, Bd, [-0.1, -0.2, 0.2, 0.1]);

    nx = size(Ad, 1); % количество строчек
    % задаем время отрисовки графиков
    TIME = 1.0;
    % задаем параметры метода решения системы дифференциальных уравнений
    options = odeset('RelTol', 1e-5, 'AbsTol', 1e-5 * ones(1, nx));
    % задаем начальные условия 
    x0 = zeros(1, nx);   

    x0(1) = 0.02;
    x0(2) = -0.03;

    discr = 0 : h : TIME;
    tlst = [];
    xlst = [;;;];
    tulst = [];
    ulst = [];

    for i = 1 : length(discr) - 1
        ticks = discr(i) : 0.001 : discr(i + 1);
        [TL, YL] = ode45(@(TL, YL)dxdt(TL, YL, theta), ticks, x0, options);
        tlst = [tlst; TL];
        xlst = [xlst; YL];
        if (i == 1)
            ulst = [ulst, 0];
        else
            ulst = [ulst, theta * x0'];
        end
        tulst = [tulst, h * (i - 1)];
        x0 = YL(end, :);
        x0 = (B * theta + eye(size(A, 1), size(A, 2))) * x0';
        x0 = x0';
    end

    fhandle = figure;
    subplot(3, 1, 1)
        plot(tlst, xlst(:, 1), 'b', 'LineWidth', 2.0)
        grid on;
        xlabel('t', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('x(t)', 'FontSize', 12, 'FontWeight', 'bold');
    subplot(3, 1, 2)
        plot(tlst, xlst(:, 2), 'b', 'LineWidth', 2.0)
        grid on;
        xlabel('t', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('\xi(t)', 'FontSize', 12, 'FontWeight', 'bold');
    subplot(3, 1, 3)
        plot(tulst, ulst(1, :), 'ob', 'LineWidth', 2.0)
        grid on;
        xlabel('t', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('u(t)', 'FontSize', 12, 'FontWeight', 'bold');
end