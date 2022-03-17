nonlin_DF();

function Y = dxdt(t, X, theta, X0)
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

    F = theta * X0;

    b = [
          F - Bc * dx - m * l * dphi^2 * sin(phi);
         -Bp * dphi + m * g * l * sin(phi)
        ];

    Y = [dx; dphi; inv(A) * b];
end

function nonlin_DF
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
    B  = [0; 0; inv(A0) * [1; 0]];

    % основная матрица системы
    A = [zeros(2, 2), eye(2); -inv(A0) * A2, -inv(A0) * A1];

    h     = 0.1;
    delta = h * 0.99;
    Ad    = expm(A * h);
    f = @(s)(expm(A * s) * B);
    Bd = integral(f, 0, delta, "ArrayValued", true);

    theta = -place(Ad, Bd, [-0.2, -0.1, 0.1, 0.2]);

    nx = size(Ad, 1); % количество строчек
    % задаем время отрисовки графиков
    TIME = 1.0;
    % задаем начальные условия 
    x0 = zeros(1, nx);   

    x0(1) = 0.02;
    x0(2) = -0.03;

    discr = 0 : h : TIME;
    tlst  = [];
    ylst  = [;;;];
    ulst  = [];
 
    y = x0';
    dt = 0.001;
    for i = 1:(length(discr) - 1)
        t = discr(i);
        tlst = [tlst; t];
        ylst = [ylst; y'];
        while (t < discr(i) + delta)
            t = t + dt;
            ulst = [ulst; theta * x0'];
            y = y + dt * dxdt(t, y, theta, x0');
            tlst = [tlst; t];
            ylst = [ylst; y'];
        end
 
        while (t < discr(i + 1))
            t = t + dt;
            ulst = [ulst; 0];
            y = y + dt * dxdt(t, y, [0, 0, 0, 0], x0');
            tlst = [tlst; t];
            ylst = [ylst; y'];
        end

        ulst = [ulst; 0];
        x0 = ylst(end, :);
    end

    fhandle = figure;
    subplot(3, 1, 1)
        plot(tlst, ylst(:, 1), 'b', 'LineWidth', 2.0)
        grid on;
        xlabel('t', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('x(t)', 'FontSize', 12, 'FontWeight', 'bold');
    subplot(3, 1, 2)
        plot(tlst, ylst(:, 2), 'b', 'LineWidth', 2.0)
        grid on;
        xlabel('t', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('\xi(t)', 'FontSize', 12, 'FontWeight', 'bold');
    subplot(3, 1, 3)
        plot(tlst, ulst(:, 1), 'b', 'LineWidth', 2.0)
        grid on;
        xlabel('t', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('u(t)', 'FontSize', 12, 'FontWeight', 'bold');
end