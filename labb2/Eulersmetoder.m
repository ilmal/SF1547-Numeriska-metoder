% Differentialekvation: 
%                       dy/dt = sin(3t) - 2y 
%                       y(0) = 1.2, t ∈ [0,8]
% --------------------------------------------------
% Uppgift 2a
% Analytisk lösning:
% y(t) = (93/65)e^(-2t)-(3/13)cos(3t)+(2/13)sin(3t)

y_analytic = @(t) (93/65)*exp(-2*t)-(3/13)*cos(3*t)+(2/13)*sin(3*t);
t_analytic = linspace(0, 8);

figure;
plot(t_analytic, y_analytic(t_analytic), 'Color', 'Red'); hold on;

% --------------------------------------------------
% Uppgift 2b

T = 8;
n = [50 100 200 400]; % antal steg h = 8/n

% Euler fram: y(i+1) = y(i) + hf(t(i), y(i))
f = @(t, y) sin(3*t)-2*y;

% Genererar 4 st distinkta färger för olika steglängd i n
colors = lines(length(n));  

% --------------------------------------------------
% för Uppgift 2c
y_errors = [];
h_errors = [];
% --------------------------------------------------

% Iterera över varje steglängd i n
for jj = 1:length(n)
    % Resetta och initiera startvärden 
    t0 = 0; 
    y0 = 1.2;
    
    % Välj steglängd från n
    n0 = n(jj);
    h = (T-t0)/n0;

    h_errors = [h_errors, h]; % uppg 2c
    
    % Initiera vektorer som sparar punkter för t och y
    tvec = zeros(n0+1, 1);
    yvec = zeros(n0+1, 1);

    % Sätt in startpunkt
    tvec(1) = t0; 
    yvec(1) = y0;

    % printa vilken iteration det är
    fprintf('\n| Börjar Euler-framåt med steglängd %d |\n', n0);
    
    % inre for-loop for Euler-fram
    for ii = 1:n0
        % Uppdatera y och t med Euler-fram
        y = y0 + h*f(t0, y0); 
        t = t0 + h;
        
        % Lägg till nya värden i vektorn for t och y punkter
        tvec(ii+1) = t;
        yvec(ii+1) = y;
        
        % Uppdatera y0 och t0
        y0 = y;
        t0 = t;

        % Printa ut resultaten
        disp([ii y y0 t])
    end
    y_errors = [y_errors, abs(y_analytic(8) - yvec(8))]; % uppg 2c

    % Plotta grafen med tvec och yvec och en godtycklig färg
    plot(tvec, yvec, 'o-', 'Color', colors(jj,:), 'LineWidth', 1.5, 'MarkerSize', 4);
end

legend('Analytisk lösning', 'Euler-framåt (50)', 'Euler-framåt (100)', ...
       'Euler-framåt (200)', 'Euler-framåt (400)', 'Location', 'best');
xlabel('t');
ylabel('y(t)');
title('Euler-framåt VS Analytisk lösning');
grid on;

% --------------------------------------------------
% Uppgift 2c
% Konvergens: Noggrannhetsordning 1

figure;
loglog(h_errors, y_errors, 'o-', 'LineWidth', 1.5)
xlabel('Steglängd h');
ylabel('Fel |y(8) - yh(8)|');
title('Konvergens av Euler-framåt');
grid on;

% --------------------------------------------------
% Uppgift 2d
% Euler bakåt 

% Plotta analytisk lösning
figure;
plot(t_analytic, y_analytic(t_analytic), 'Color', 'Red'); hold on;

Tb = 8;
nb = [50 100 200 400]; % antal steg h = 8/n

% Euler fram: y(i+1) = y(i) + hf(t(i), y(i))
f = @(t, y) sin(3*t)-2*y;

% Genererar 4 st distinkta färger för olika steglängd i n
colorsb = lines(length(n));  

% Räkna fel
y_errorsb = [];
h_errorsb = [];


% Iterera över varje steglängd i n
for jj = 1:length(nb)
    % Resetta och initiera startvärden 
    t0 = 0; 
    y0 = 1.2;
    
    % Välj steglängd från n
    n0 = nb(jj);
    h = (Tb-t0)/n0;

    h_errorsb = [h_errorsb, h]; % uppg 2c
    
    % Initiera vektorer som sparar punkter för t och y
    tvec = zeros(n0+1, 1);
    yvec = zeros(n0+1, 1);

    % Sätt in startpunkt
    tvec(1) = t0; 
    yvec(1) = y0;

    % printa vilken iteration det är
    fprintf('\n| Börjar Euler-framåt med steglängd %d |\n', n0);
    
    % inre for-loop for Euler-fram
    for ii = 1:n0
        % Uppdatera y och t med Euler-fram
        y = y0 + h*f(t, y); 
        t = t0 + h;
        
        % Lägg till nya värden i vektorn for t och y punkter
        tvec(ii+1) = t;
        yvec(ii+1) = y;
        
        % Uppdatera y0 och t0
        y0 = y;
        t0 = t;

        % Printa ut resultaten
        disp([ii y y0 t])
    end
    y_errorsb = [y_errorsb, abs(y_analytic(8) - yvec(8))]; 

    % Plotta grafen med tvec och yvec och en godtycklig färg
    plot(tvec, yvec, 'o-', 'Color', colorsb(jj,:), 'LineWidth', 1.5, 'MarkerSize', 4);
end

legend('Analytisk lösning', 'Euler-bakåt (50)', 'Euler-bakåt (100)', ...
       'Euler-bakåt (200)', 'Euler-bakåt (400)', 'Location', 'best');
xlabel('t');
ylabel('y(t)');
title('Euler-bakåt VS Analytisk lösning');
grid on;

