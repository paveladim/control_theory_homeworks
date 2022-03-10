function pendulum_Deltafixed_control
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
    delta = h * 0.3;
    Ad = expm(A * h);
    Bd = expm(A * h) * B * delta;

    theta = -place(Ad, Bd, [-0.2, -0.1, 0.1, 0.2]);

    nx = size(Ad, 1); % количество строчек
    % задаем время отрисовки графиков
    TIME = 1.0;
    % задаем параметры метода решения системы дифференциальных уравнений
    options = odeset('RelTol', 1e-5, 'AbsTol', 1e-5 * ones(1, nx));
    % задаем начальные условия 
    x0 = zeros(1, nx);   

    x0(1) = 0.0;
    x0(2) = -0.3;

    discr = 0 : h : TIME;
    tlst = [];
    xlst = [;;;];
    ulst = [];

    ticks = discr(1) : 0.001 : discr(2);
    [TL, YL] = ode45(@(t, X)(A * X), ticks, x0, options);
    tlst = [tlst; TL];
    xlst = [xlst; YL];
    ulst = [ulst, 0 * ones(1, size(TL, 1))];
    x0 = YL(end, :);

    for i = 2 : size(discr, 2) - 1
        ticks = discr(i) : 0.001 : discr(i + 1);
        [TL, YL] = ode45(@(t, X)(A * X + B * theta * x0' * (heaviside(t - discr(i)) - heaviside(t - discr(i) - delta))), ticks, x0, options);
        tlst = [tlst; TL];
        xlst = [xlst; YL];
        ulst = [ulst, theta * x0' * ones(1, size(TL, 1))];
        x0 = YL(end, :);
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
        plot(tlst, ulst(1, :), 'b', 'LineWidth', 2.0)
        grid on;
        xlabel('t', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('u(t)', 'FontSize', 12, 'FontWeight', 'bold');
end