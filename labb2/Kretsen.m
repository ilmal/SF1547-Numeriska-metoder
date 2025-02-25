format long;
% -----------------------------
% Uppgift 3a

% Differentialekvation: 
% Lq''(t) + Rq'(t) + (1/C)q(t) = 0 
% q(0) = 1 
% dq/dt (0) = 0

% Strömmen ges av: 
% i = dq/dt

% Andraderivatan skrivs om till: 
% Lq'' + Rq' + (1/C)q = 0 
% Lq'' = -(1/C)q(t)-Rq'
% q'' = -((1/C)q - Rq')/L

% Omskrivning:
% q'' = f(t, q, q') 
% -> u(t) = q(t), u(0) = 1
%    v(t) = q'(t), v(0) = 0

% -> u' = q' = I, u'(0) = 0
%    v' = q'' = f(t, q, q') = -((1/C)u - Rv)/L

% y = [u; v];
% F(t, U) = [u'; v'];
% F(t,U) = [   v; 
%     -((1/C)*u + R*v)/L]; 
% U0 = [1; 0];

% -----------------------------
% Uppgift 3b

function F = returnVectorF(t, y, R, L, C)
    % substituera så u = q, v = q'
    u = y(1);
    v = y(2);

    % Beräknad från 3a
    F = [v;
        -((1/C)*u + R*v)/L
    ];
end

% -----------------------------
% Uppgift 3c

% startvärden
y0 = [1; 0];
tspan = [0 20];

% y(:,1) → Innehåller värden för laddningen q(t)
% y(:,2) → Innehåller värden för laddningen i(t)

% ode45(funktion, integrationsgräns, begynnelsevillkor)
% diffekv av formen y' = f(t, y) -> kalla på returnvector med givna RLC-värden

% Dämpad svängning 
[t, y] = ode45(@(t, y) returnVectorF(t, y, 1, 2, 0.5), tspan, y0);

% Plotta resultatet
figure;
plot(t, y(:,1), 'b-', 'LineWidth', 1.5); hold on;
plot(t, y(:,2), 'r--', 'LineWidth', 1.5);
xlabel('Tid t');
ylabel('Q och I');
legend('Laddning Q', 'Ström I');
title('Lösning av RLC-kretsen med odämpad svängning');
grid on;

% Odämpad svängning
[t, y] = ode45(@(t, y) returnVectorF(t, y, 0, 2, 0.5), tspan, y0);

% plotta resultat
figure;
plot(t, y(:,1), 'b-', 'LineWidth', 1.5); hold on;
plot(t, y(:,2), 'r--', 'LineWidth', 1.5);
xlabel('Tid t');
ylabel('Q och I');
legend('Laddning Q', 'Ström I');
title('Lösning av RLC-kretsen med odämpad svängning');
grid on;

% -----------------------------
% Uppgift 3d
% Euler fram för system

% Givna värden på R, L, C
L = 2; C = 0.5; R = 1;

% Vektor F med derivatorna u' och v'
F = @(t, U) [U(2); -((1/C)*U(1) + R*U(2))/L];
U0 = [1; 0];

N = [40 80 160 320]; % antal tidssteg
t0 = 0; Tslut = 40;

% disp(F);

figure; 
for ii = 1:length(N)
    N0 = N(ii);
    h = (Tslut-t0)/N0; % steglängd
    t = t0;
    U = U0;

    % vektorer för att spara alla beräknade punkter (n0+1 rader, 1 kolumn)
    tvec = zeros(1, N0+1);
    Uvec = zeros(2, N0+1);

    Uvec(:,1) = U; % sätt 1st kolumn i matris till U0 värden
    tvec(1) = t0;

    for jj = 1:N0
        U = U + h*F(t0, U);
        t = t+h;

        Uvec(:, jj+1) = U;
        tvec(jj+1) = t;
    end

    subplot(2, 2, ii);
    tspan2 = [0 40];

    % Plotta ode45 lösn
    [t, y] = ode45(@(t, y) returnVectorF(t, y, 1, 2, 0.5), tspan2, y0);
    plot(t, y(:,1), 'b-', 'LineWidth', 1.5); hold on;
    plot(t, y(:,2), 'r-', 'LineWidth', 1.5); hold on;
    
    % plotta euler lösn
    plot(tvec, Uvec(1,:), '-o', tvec, Uvec(2,:), '-or')

    % formattering + text o sånt
    title(sprintf('Euler Method vs. ODE45 (N = %d)', N0));
    xlabel('Tid t');
    ylabel('Laddning Q och Ström I');
    legend({'ODE45 Q(t)', 'ODE45 I(t)', 'Euler Q(t)', 'Euler I(t)'}, 'Location', 'Best');
    grid on;

end