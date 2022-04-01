function class2603
    A  = [0, -1; 0, -1];
    Bu = [0; 1];
    Bw = [0; 1];
    
    h = 0.1;
    sigmaw = 1;

    Ad = expm(A * h);
    f  = @(s)(expm(-A * s) * Bu);
    Bd = integral(f, 0, h, "ArrayValued", true);

    g = @(s)(expm(A * s) * Bw * Bw' * expm(A' * s));
    W = sigmaw * integral(g, 0, h, "ArrayValued", true);

    Q = eye(2, 2);
    R = 1;

    TIME = 0 : h : 50;
    len = length(TIME);

    matrixes   = [];
    regulators = [];
    XN = zeros(2, 2);
    U = -inv(R + Bd' * XN * Bd) * Bd' * XN * Ad;
    matrixes   = [matrixes;  XN];
    regulators = [regulators; U];

    for i = 1 : len - 1
        XN1 = Q + Ad' * XN * Ad - Ad' * XN * Bd * inv(R + Bd' * XN * Bd) * Bd' * XN * Ad;
        XN = XN1;
        U = -inv(R + Bd' * XN * Bd) * Bd' * XN * Ad;
        matrixes   = [matrixes;  XN];
        regulators = [regulators; U];
    end

    regulators = flipud(regulators);

    % сгенерируем случайные величины
    vk = [];
    for i = 1 : len - 1
        s = rng;
        % RR = chol(sigmaw)
        matr = randn(100, 100);
        tmp = matr(23, 71);
        %tmp = repmat(0, 50000, 1) + tmp * sqrt(2) * RR;
        vk = [vk, tmp];
    end

    nx = size(Ad, 1); % количество строчек
    % задаем параметры метода решения системы дифференциальных уравнений
    options = odeset('RelTol', 1e-5, 'AbsTol', 1e-5 * ones(1, nx));
    % задаем начальные условия 
    x0 = zeros(1, nx);   

    x0(1) = 0.0;
    x0(2) = -0.3;

    tlst = [];
    xlst = [];
    ulst = [];

    for i = 1 : len - 1
        ticks = TIME(i) : 0.001 : TIME(i + 1);
        theta = regulators(i, :);
        [TL, YL] = ode45(@(t, X)(A * X + Bu * theta * x0' + Bw * vk(i)), ticks, x0, options);
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
        ylabel('\phi(t)', 'FontSize', 12, 'FontWeight', 'bold');
    subplot(3, 1, 2)
        plot(tlst, xlst(:, 2), 'b', 'LineWidth', 2.0)
        grid on;
        xlabel('t', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('\psi(t)', 'FontSize', 12, 'FontWeight', 'bold');
    subplot(3, 1, 3)
        plot(tlst, ulst(1, :), 'b', 'LineWidth', 2.0)
        grid on;
        xlabel('t', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('u(t)', 'FontSize', 12, 'FontWeight', 'bold');

end