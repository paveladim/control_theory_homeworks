function hwct260221
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

    % вектор C
    C = [0, 0, 1, 0; 0, 0, 0, 1];

    % матрица наблюдаемости
    OBS = [C; C * A; C * A * A; C * A * A * A];
    
    if (rank(OBS) == size(A, 1))
        disp("system is observed");
    else
        disp("system is unobserved");
    end

    % вектор C
    C = [1, 0, 0, 0; 0, 1, 0, 0];

    % матрица наблюдаемости
    OBS = [C; C * A; C * A * A; C * A * A * A];
    
    if (rank(OBS) == size(A, 1))
        disp("system is observed");
    else
        disp("system is unobserved");
    end

    newA = A';
    newB = C';
    newB = newB(:, 1);
    
    CTR = [newB, newA * newB, newA * newA * newB, newA * newA * newA * newB];
    
    chiA = newA^4 + 18.0189 * newA^3 + 104.6974139 * newA^2 + 193.5929306 * newA;
    theta = -[0, 0, 0, 1] * inv(CTR) * chiA;
    L = -theta';
    L = [L, zeros(size(L, 1), 1)];

    AA = [A, zeros(4, 4) + B * theta; L * C, A - L * C + B * theta];
    
    nx = size(AA, 1); % количество строчек
    nu = size(B, 2); % количество столбцов

    % задаем время отрисовки графиков
    TIME = 10.0;
    % задаем точки, в которых нужно вычислять решение системы
    ticks = 0 : 0.001 : TIME;
    % задаем параметры метода решения системы дифференциальных уравнений
    options = odeset('RelTol', 1e-5, 'AbsTol', 1e-5 * ones(1, nx));
    % задаем начальные условия 
    X0 = zeros(1, nx);   

    X0(1) = 0.1;
    X0(2) = -0.1;

    % находим решение замкнутой линейной системы
    [TL, YL] = ode45(@(t, X)(AA * X), ticks, X0, options);

    % --------------------------------------------------------------------
	% СТРОИМ ГРАФИКИ РЕШЕНИЙ
    % --------------------------------------------------------------------  

    label1 = 'x(t)';
    label2 = '\dot{x}(t)';
    label3 = '\xi(t)';
    label4 = '\dot{\xi}(t)';

    fhandle = figure;
    subplot(2, 2, 1)
        plot(TL, YL(:, 1), 'b', TL, YL(:, 5), 'r', 'LineWidth', 2.0)
        grid on;
        xlabel('t', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel(label1, 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('x_1^0 = %0.3f', X0(1)));
    subplot(2, 2, 2)
        plot(TL, YL(:, 2), 'b', TL, YL(:, 6), 'r', 'LineWidth', 2.0)
        grid on;
        xlabel('t', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel(label2, 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('x_2^0 = %0.3f', X0(2)));        
    subplot(2, 2, 3)
        plot(TL, YL(:, 3), 'b', TL, YL(:, 7), 'r', 'LineWidth', 2.0)
        grid on;
        xlabel('t', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel(label3, 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('x_2^0 = %0.3f', X0(3)));  
    subplot(2, 2, 4)
        plot(TL, YL(:, 4), 'b', TL, YL(:, 8), 'r', 'LineWidth', 2.0)
        grid on;
        xlabel('t', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel(label4, 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('x_2^0 = %0.3f', X0(4)));  

end