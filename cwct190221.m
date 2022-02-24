function cwct190221

    A = [1, 0; 0, 0];
    B = [1; 1];
    C = [1, 1];
    
    CM = [C', A * C'];
    chiA = A^2 + 6 * A + 409 * A^0;
    
    theta = -[0, 1] * inv(CM) * (chiA);
    L = -theta';
    
    AC = [A, B * theta; L * C, A - L * C + B * theta];

    nx = size(AC, 1); % количество строчек
    nu = size(B, 2); % количество столбцов

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
    [TL, YL] = ode45(@(t, X)(AC * X), ticks, X0, options);


    % --------------------------------------------------------------------
	% СТРОИМ ГРАФИКИ РЕШЕНИЙ
    % --------------------------------------------------------------------  

    fhandle = figure;
    subplot(2, 1, 1)
        plot(TL, YL(:, 1), 'b', TL, YL(:, 3), 'r', 'LineWidth', 2.0)
        grid on;
        xlabel('t', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('x(t)', 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('x_1^0 = %0.3f', X0(1)));
    subplot(2, 1, 2)
        plot(TL, YL(:, 2), 'b', TL, YL(:, 4), 'r', 'LineWidth', 2.0)
        grid on;
        xlabel('t', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('\xi(t)', 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('x_2^0 = %0.3f', X0(2)));               

end
