function hwct190221
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
    % собственные числа матрицы
    eigen_values_1 = eig(A);
    [V, D] = eig(A);
    disp(V);
    vec = V(:, 4);
    disp(A * vec);
    % характеристический многочлен
    p = poly(A);

    % матрица управляемости
    C = [B, A * B, A^2 * B, A^3 * B];

    if (rank(C) == size(A, 1)) 
        disp("system is controllable")
    else
        disp("system is uncontrollable")
    end

    % переместим корни 0 -> -4; 6.6011 -> -6
    ChiA = A^4 + 21.4178 * A^3 + 167.5053743 * A^2 + 567.3009432 * A + 703.8569837 * A^0;
    theta = -[0, 0, 0, 1] * inv(C) * ChiA;
    disp('regulator');
    disp(theta);

    Ac = A + B * theta;
    eigen_values_ac = eig(Ac);
    disp(eigen_values_ac);

    Ac = [0, 0, 1, 0; 0, 0, 0, 1; -131.4332, 14.4426, -26.3799, 25.3925;
          -40.1635, 3.3355, -5.9868, 5.8642];

    % --------------------------------------------------------------------
    % ПРИМЕНЯЕМ ПОСТРОЕННЫЕ РЕГУЛЯТОРЫ К ЛИНЕЙНОЙ СИСТЕМЕ
    % --------------------------------------------------------------------    
    
    nx = size(A, 1); % количество строчек
    nu = size(B, 2); % количество столбцов

    % задаем время отрисовки графиков
    TIME = 3.5;
    % задаем точки, в которых нужно вычислять решение системы
    ticks = 0 : 0.001 : TIME;
    % задаем параметры метода решения системы дифференциальных уравнений
    options = odeset('RelTol', 1e-5, 'AbsTol', 1e-5 * ones(1, nx));
    % задаем начальные условия 
    X0 = zeros(1, nx);   
    X0(2) = -0.4;

    % находим решение замкнутой линейной системы
    [TL, YL] = ode45(@(t, X)((A + B * theta) * X), ticks, X0, options); 

    for k = 1 : length(YL)
        UL(:, k) = theta * YL(k, 1 : nx)';
    end


    % --------------------------------------------------------------------
	% СТРОИМ ГРАФИКИ РЕШЕНИЙ
    % --------------------------------------------------------------------  

    fhandle = figure;
    subplot(3, 1, 1)
        plot(TL, YL(:, 1), 'b', 'LineWidth', 2.0)
        grid on;
        xlabel('t', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('x(t)', 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('x_1^0 = %0.3f', X0(1)));
    subplot(3, 1, 2)
        plot(TL, YL(:, 2), 'b', 'LineWidth', 2.0)
        grid on;
        xlabel('t', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('\phi(t)', 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('x_2^0 = %0.3f', X0(2)));            
    subplot(3, 1, 3)
        plot(TL, UL(1, :), 'b', 'LineWidth', 2.0) 
        grid on;
        xlabel('t', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('u(t)', 'FontSize', 12, 'FontWeight', 'bold');   

end