function types_of_control
    h = 0.5;
    A = [1, h; 0, 1];
    B = [0.5 * h^2; h];

    C = [0.5 * h^2, 1.5 * h^2; h, h];
    C = inv(C);
    chiA = A^2 - 0.25 * A^0;
    theta = -[0, 1] * C * chiA;
    theta = -place(A, B, [-0.1, -0.2]);

    disp(eig(A + B * theta));

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
    xlst = [;];

%     for i = 1 : size(discr, 2) - 1
%         ticks = discr(i) : 0.001 : discr(i + 1);
%         [TL, YL] = ode45(@(t, X)(A * X + B * theta * x0), ticks, x0', options);
%         tlst = [tlst; TL];
%         xlst = [xlst; YL];
%         x0 = YL(end, :)';
%     end

    tlst = [];
    xlst = [;];
    eps = 0.001;
    for i = 1 : size(discr, 2) - 1
        tend = discr(i) + eps * h;
        ticks = discr(i) : 0.00001 : tend;
        [TL, YL] = ode45(@(t, X)(A * X + B * theta * x0), ticks, x0', options);
        tlst = [tlst; TL];
        xlst = [xlst; YL];
        x0 = YL(end, :)';

        ticks = tend : 0.001 : discr(i + 1);
        [TL, YL] = ode45(@(t, X)(A * X), ticks, x0', options);
        tlst = [tlst; TL];
        xlst = [xlst; YL];
        x0 = YL(end, :)';
    end
% 
%     tlst = [];
%     xlst = [;];
%     delta = h * 0.25;
%     for i = 1 : size(discr, 2) - 1
%         ticks = discr(i) : 0.001 : (discr(i) + delta);
%         [TL, YL] = ode45(@(t, X)(A * X + B * theta * x0), ticks, x0', options);
%         tlst = [tlst; TL];
%         xlst = [xlst; YL];
%         x0 = YL(end, :)';
% 
%         ticks = (discr(i) + delta) : 0.001 : discr(i + 1);
%         [TL, YL] = ode45(@(t, X)(A * X), ticks, x0', options);
%         tlst = [tlst; TL];
%         xlst = [xlst; YL];
%         x0 = YL(end, :)';
%     end

    fhandle = figure;
    subplot(2, 1, 1)
        plot(tlst, xlst(:, 1), 'b', 'LineWidth', 2.0)
        grid on;
        xlabel('t', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('x(t)', 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('x_1^0 = %0.3f', x0(1)));
    subplot(2, 1, 2)
        plot(tlst, xlst(:, 2), 'b', 'LineWidth', 2.0)
        grid on;
        xlabel('t', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('\xi(t)', 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('x_2^0 = %0.3f', x0(2)));

end