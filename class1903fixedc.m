function class1903fixedc
    A = [0, 1; 1, 0];
    B = [0; 1];
    h = 0.1;
    Ad = expm(A * h);
    f = @(s)(expm(A * s) * B);
    Bd = integral(f, 0, h, "ArrayValued", true);

    theta = -place(Ad, Bd, [-0.5, 0.5]);

    nx = size(Ad, 1); % количество строчек
    % задаем время отрисовки графиков
    TIME = 1.0;
    % задаем параметры метода решения системы дифференциальных уравнений
    options = odeset('RelTol', 1e-5, 'AbsTol', 1e-5 * ones(1, nx));
    % задаем начальные условия 
    x0 = zeros(1, nx);   

    x0(1) = 0.0;
    x0(2) = -0.6;

    discr = 0 : h : TIME;
    tlst = [];
    xlst = [];
    ulst = [];

    x0lst = [];

    for i = 1 : size(discr, 2) - 1
        ticks = discr(i) : 0.001 : discr(i + 1);
        [TL, YL] = ode45(@(t, X)(A * X + B * theta * x0'), ticks, x0, options);
        tlst = [tlst; TL];
        xlst = [xlst; YL];
        ulst = [ulst, theta * x0' * ones(1, size(TL, 1))];
        x0 = YL(end, :);
        x0lst = [x0lst; x0];
    end

    fhandle = figure;
    subplot(2, 2, 1)
        plot(tlst, xlst(:, 1), 'g', 'LineWidth', 2.0)
        grid on;
        xlabel('t', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('x(t)', 'FontSize', 12, 'FontWeight', 'bold');
    subplot(2, 2, 2)
        plot(tlst, xlst(:, 2), 'g', 'LineWidth', 2.0)
        grid on;
        xlabel('t', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('y(t)', 'FontSize', 12, 'FontWeight', 'bold');
    subplot(2, 2, 3)
        plot(xlst(:, 1), xlst(:, 2), 'r', x0lst(:, 1), x0lst(:, 2), 'ob', 'LineWidth', 1.0)
        grid on;
        xlabel('x(t)', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('y(t)', 'FontSize', 12, 'FontWeight', 'bold');
    subplot(2, 2, 4)
        plot(tlst, ulst(1, :), 'r', 'LineWidth', 1.0)
        grid on;
        xlabel('t', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('u(t)', 'FontSize', 12, 'FontWeight', 'bold');
end