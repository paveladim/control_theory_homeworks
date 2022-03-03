function cw260221
    h = 0.5;
    A = [1, h; 0, 1];
    B = [0.5 * h^2; h];

    C = [0.5 * h^2, 1.5 * h^2; h, h];
    C = inv(C);
    chiA = A^2 - 0.25 * A^0;
    theta = -[0, 1] * C * chiA;
    theta1 = -place(A, B, [-0.5, 0.5]);

    %theta = [-0.75 / h^2, -0.5 / h + 0.75 * 0.5 / h];

    disp(eig(A + B * theta));
    AC = A + B * theta;

    nx = size(AC, 1); % количество строчек

    x = [0.1; 0.2];
    vec(:, 1) = x;
    ticks(:, 1) = 0;
 
    i = 2;
    while i < 50
        x = AC * x;
        vec(:, i) = x;
        ticks(:, i) = ticks(:, i - 1) + h;
        i = i + 1;
    end
    
    fhandle = figure;
    subplot(2, 1, 1)
        plot(ticks, vec(1, :), 'b', 'LineWidth', 2.0)
        grid on;
        xlabel('t', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('x(t)', 'FontSize', 12, 'FontWeight', 'bold');
    subplot(2, 1, 2)
        plot(ticks, vec(2, :), 'b', 'LineWidth', 2.0)
        grid on;
        xlabel('t', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('\xi(t)', 'FontSize', 12, 'FontWeight', 'bold');     

end