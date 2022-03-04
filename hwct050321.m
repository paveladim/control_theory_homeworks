function hwct050321
    m  = 0.127; % масса маятника
    M  = 1.206; % масса тележки
    I  = 0.001; % момент инерции маятника относительно центра масс
    l  = 0.178; % расстояние от точки крепления до центра масс
    Bc = 5.4;   % коэф. вязкого трения между кареткой и направляющей
    Bp = 0.002; % коэф. вязкого трения в точке крепления
    g  = 9.81;  % коэф. свободного падения

    % определим матрицы системы
    A0 = [m + M , -m * l; 
          -m * l, I + m * l^2];

    A1 = diag([Bc, Bp]);
    A2 = diag([0, -m * g * l]);
    B = [0; 0; inv(A0) * [1; 0]];

    % основная матрица системы
    A = [zeros(2, 2), eye(2); -inv(A0) * A2, -inv(A0) * A1];

    h = 0.1;
    Ad = exp(A * h);
    f = @(s)(exp(A * s) * B);

    Bd = [0; 0; 0; 0];
    step = h / 10000;
    for i = 1:10000
        Bd = Bd + f(step * (i - 1)) * step;
    end
    
    theta = -place(A, B, [-0.1, -0.2, 0.2, 0.1]);
    
    nx = size(A, 1); % количество строчек
    % задаем время отрисовки графиков
    TIME = 10.0;
    % задаем параметры метода решения системы дифференциальных уравнений
    options = odeset('RelTol', 1e-5, 'AbsTol', 1e-5 * ones(1, nx));
    % задаем начальные условия 
    x0 = zeros(nx, 1);   

    x0(1) = 0.5;
    x0(2) = -0.4;

    discr = 0 : h : TIME;
    tlst = [];
    xlst = [;;;];

    for i = 1 : size(discr, 2) - 1
        ticks = discr(i) : 0.001 : discr(i + 1);
        [TL, YL] = ode45(@(t, X)(A * X + B * theta * x0), ticks, x0', options);
        tlst = [tlst; TL];
        xlst = [xlst; YL];
        x0 = YL(end, :)';
    end

    fhandle = figure;
    subplot(2, 1, 1)
        plot(tlst, xlst(:, 1), 'b', 'LineWidth', 2.0)
        grid on;
        xlabel('t', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('x(t)', 'FontSize', 12, 'FontWeight', 'bold');
    subplot(2, 1, 2)
        plot(tlst, xlst(:, 2), 'b', 'LineWidth', 2.0)
        grid on;
        xlabel('t', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('\xi(t)', 'FontSize', 12, 'FontWeight', 'bold');
end