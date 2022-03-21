function class1903
    A = [0, 1; 1, 0];
    B = [0; 1];

    h = log(2);

    Ad = [cosh(h), sinh(h); sinh(h), cosh(h)];
    Bd = [sinh(h); cosh(h)];

    theta = -place(Ad, Bd, [-0.5, 0.5]);

    nx = size(Ad, 1); % количество строчек
    % задаем время отрисовки графиков
    TIME = 10.0;
    % задаем параметры метода решения системы дифференциальных уравнений
    options = odeset('RelTol', 1e-5, 'AbsTol', 1e-5 * ones(1, nx));
    % задаем начальные условия 
    x0 = zeros(1, nx);   

    x0(1) = 0.02;
    x0(2) = 0.03;

    discr = 0 : h : TIME;
    tlst = [];
    xlst = [];
    tulst = [];
    ulst = [];
    x0lst = [];

    for i = 1 : size(discr, 2) - 1
        ticks = discr(i) : 0.001 : discr(i + 1);
        [TL, YL] = ode45(@(t, X)(A * X), ticks, x0, options);
        tlst = [tlst; TL];
        xlst = [xlst; YL];
        if (i == 1)
            ulst = [ulst, 0];
        else
            ulst = [ulst, theta * x0'];
        end
        tulst = [tulst, h * (i - 1)];
        x0 = YL(end, :);
        x0lst = [x0lst, YL(end, :)'];
        x0 = (B * theta + eye(size(A, 1), size(A, 2))) * x0';
        x0lst = [x0lst, x0];
        x0 = x0';
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
        plot(xlst(:, 1), xlst(:, 2), 'r', x0lst(1, :), x0lst(2, :), 'ob', 'LineWidth', 1.0)
        grid on;
        xlabel('x(t)', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('y(t)', 'FontSize', 12, 'FontWeight', 'bold');
    subplot(2, 2, 4)
        plot(tulst, ulst(1, :), 'or', 'LineWidth', 1.0)
        grid on;
        xlabel('t', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('u(t)', 'FontSize', 12, 'FontWeight', 'bold');

end