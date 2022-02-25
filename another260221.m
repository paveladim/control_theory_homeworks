function another260221
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
    
    % вектор выхода
    C = [0, 0, 1, 0; 0, 0, 0, 1];

    % перейдём к двойственной системе, чтобы построить наблюдатель
    A = A';
    C = C';
    B = B';

    % матрица управляемости
    controlM = [C, A * C, A^2 * C, A^3 * C];
    
    if (rank(controlM) == size(A, 1))
        disp("system is contollable");
    else
        disp("system is uncontrollable")
    end

    vec = [1; 0; 0; 0];
    controlM1 = [controlM(:, 1), controlM(:, 3), controlM(:, 5), vec];
    controlM2 = [controlM(:, 2), controlM(:, 4), controlM(:, 6), vec];
    
    % передвинем собственные числа итеративно, если это возможно
    [V, D] = eig(A');
    vec = V(:, 2);
    
    S = [vec'; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
    newA = S * A * inv(S);
    newC = S * C;
    newA = newA + newC * [0, 0, 0, 0; 15, 0, 0, 0];

    theta = [0, 0, 0, 0; 15, 0, 0, 0] * S;
    
    % сам наблюдатель
    L = -theta';
    theta = [19.6330, -64.6908, 21.4010, -9.2987];

    A = A';
    B = B';
    C = C';
    AC = [A, zeros(4, 4) + B * theta; L * C, A - L * C + B * theta];

    nx = size(AC, 1); % количество строчек

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
        plot(TL, YL(:, 3), 'b', TL, YL(:, 7), 'r', 'LineWidth', 2.0)
        grid on;
        xlabel('t', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('x(t)', 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('x_1^0 = %0.3f', X0(1)));
    subplot(2, 1, 2)
        plot(TL, YL(:, 4), 'b', TL, YL(:, 8), 'r', 'LineWidth', 2.0)
        grid on;
        xlabel('t', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('\xi(t)', 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('x_2^0 = %0.3f', X0(2)));               


end