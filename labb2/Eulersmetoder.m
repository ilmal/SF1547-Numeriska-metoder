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
plot(t_analytic, y_analytic(t_analytic), 'Color', 'Black', 'LineWidth', 2.0); hold on;

% --------------------------------------------------
% Uppgift 2b

T = 8; % slutvärde i intervallet [0, 8] = [t0, T]
n = [50 100 200 400]; % några olika antal steg -> h = T-t0/n

% Euler fram: y(i+1) = y(i) + hf(t(i), y(i))
f = @(t, y) sin(3*t)-2*y;

% Genererar 4 st distinkta färger för olika steglängd i n
colors = lines(length(n));  

% --------------------------------------------------
% ...för Uppgift 2c
y_errors = [];
h_errors = [];
% --------------------------------------------------

% Iterera över varje olika steglängd n (50, 100, 200, 400)
for jj = 1:length(n)
    % Resetta och initiera startvärden 
    t0 = 0; % Startvärde i intervall [0, 8]
    y0 = 1.2; % givet startvärde av uppg y(0) = 1.2
    
    % Beräkna steglängd för euler fram
    n0 = n(jj); % Välj steglängd från n med index jj
    h = (T-t0)/n0; % beräkna h

    h_errors = [h_errors, h]; % uppg 2c
    
    % Initiera vektorer som sparar punkter för t och y
    tvec = zeros(n0+1, 1); % Kolumnvektor för t med n+1 nollor
    yvec = zeros(n0+1, 1); % Kolumnvektor för y med n+1 nollor

    % Sätt in startpunktens värden
    tvec(1) = t0; 
    yvec(1) = y0;

    % printa vilken iteration det är
    % fprintf('\n| Börjar Euler-framåt med steglängd %d |\n', n0);
    
    % Inre for-loop for Euler-fram
    % Iterera från 1 till värdet på n0 (50 100 200 400 beroende på omgång)
    for ii = 1:n0
        % Uppdatera y och t med Euler-fram
        y = y0 + h*f(t0, y0); % Uppdatera y
        t = t0 + h; % Uppdatera t
        
        % Lägg till nya värden i vektorn for t och y punkter
        tvec(ii+1) = t; % Sätt nytt t-värde i tvec
        yvec(ii+1) = y; % Sätt nytt y-värde i yvec
        
        % Uppdatera y0 och t0
        y0 = y; % Sätt gamla y till nya y
        t0 = t; % Sätt gamla t till nya t

        % Printa ut resultaten
        % disp([ii y y0 t])
    end
    y_errors = [y_errors, abs(y_analytic(T) - yvec(end))]; % uppg 2c

    % Plotta grafen med tvec och yvec och en godtycklig färg (med given n0 för iterationen)
    plot(tvec, yvec, 'o-', 'Color', colors(jj,:), 'LineWidth', 1.5, 'MarkerSize', 4);
end

% Formatgrejer
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

% loglog sätter naturliga logaritmen automatiskt på båda axlar
% eh = Ch^p -> ln eh = p * ln h + ln C
loglog(h_errors, y_errors, 'o-', 'LineWidth', 1.5)
xlabel('Steglängd h');
ylabel('Fel |y(8) - yh(8)|');
title('Konvergens av Euler-framåt');
grid on;

% Polyfit ger koefficienter på polynom av grad 1 i detta fall
% ln eh = p * ln h + ln C -> p(1) ger p, om p ≈ 1 så noggrannhet 1
p = polyfit(log(h_errors), log(y_errors), 1);
disp(['Noggrannhetsordning Euler-framåt, p ≈ ', num2str(p(1))]);

% --------------------------------------------------
% Uppgift 2d
% Euler bakåt 

% Plotta analytisk lösning
figure;
plot(t_analytic, y_analytic(t_analytic), 'Color', 'Black'); hold on;

% Startvärden
Tb = 8; % Slutvärde i intervallet [0, 8] = [t0, T]
nb = [50 100 200 400]; % antal steg h = 8/n

% Genererar 4 st distinkta färger för olika steglängd i n
colorsb = lines(length(nb));  

% Räkna fel
y_errorsb = [];
h_errorsb = [];


% Iterera över varje steglängd i n
for jj = 1:length(nb)
    % Resetta och initiera startvärden 
    t0 = 0; % Startvärde på intervall
    y0 = 1.2; % Givet startvärde på y(0)

    % Välj steglängd från n
    n0 = nb(jj);
    h = (Tb-t0)/n0;

    h_errorsb = [h_errorsb, h]; % Lägg till nuvarande steglängd h i errorvektor

    % Initiera vektorer som sparar punkter för t och y
    tvec = zeros(n0+1, 1);
    yvec = zeros(n0+1, 1);

    % Sätt in startpunkt
    tvec(1) = t0; 
    yvec(1) = y0;

    % printa vilken iteration det är
    % fprintf('\n| Börjar Euler-bakåt med steglängd %d |\n', n0);

    % inre for-loop for Euler-fram
    for ii = 1:n0
        % Nya t = t(i+1)
        t = t0 + h;
        
        % BERÄKNING AV EULER BAK FORMEL
        % y(i+1) = y(i) + h*sin(3t(i+1)) - h*2*y(i+1)
        % y(i+1) + 2h*y(i+1) = y(i) + h*sin(3t(i+1))
        % y(i+1)(1+2h) = y(i) + h*sin(3t(i+1))
        % y(i+1) = y(i) + h*sin(3t(i+1)) / (1+2h)

        % Uppdatera y med Euler-bakåt
        y = (y0 + h*sin(3*t)) / (1+2*h); 

        % Lägg till nya värden i vektorn for t och y punkter
        tvec(ii+1) = t;
        yvec(ii+1) = y;

        % Uppdatera y0 och t0
        y0 = y;
        t0 = t;

        % Printa ut resultaten
        % disp([ii y y0 t])
    end
    % Beräkna felet
    y_errorsb = [y_errorsb, abs(y_analytic(Tb) - yvec(end))]; 

    % Plotta grafen med tvec och yvec och en godtycklig färg
    plot(tvec, yvec, 'o-', 'Color', colorsb(jj,:), 'LineWidth', 1.5, 'MarkerSize', 4);
end

legend('Analytisk lösning', 'Euler-bakåt (50)', 'Euler-bakåt (100)', ...
       'Euler-bakåt (200)', 'Euler-bakåt (400)', 'Location', 'best');
xlabel('t');
ylabel('y(t)');
title('Euler-bakåt VS Analytisk lösning');
grid on;

% Plotta felen
figure;
loglog(h_errorsb, y_errorsb, 'o-', 'LineWidth', 1.5)
xlabel('Steglängd h');
ylabel('Fel |y(8) - yh(8)|');
title('Konvergens av Euler-bakåt');
grid on;

pb = polyfit(log(h_errorsb), log(y_errorsb), 1);
disp(['Noggrannhetsordning Euler-bakåt, p ≈ ', num2str(p(1))]);

% --------------------------------------------------
% Uppgift 2e
% Euler framåt med [0, 80] och n = [50 100 400 800]

T_e = 80;
n_e = [50 100 400 800]; % antal steg h = 8/n

% Genererar 4 st distinkta färger för olika steglängd i n
colors = lines(length(n_e));  

% Iterera över varje steglängd i n
figure;
for jj = 1:length(n_e)
    % Resetta och initiera startvärden 
    t0 = 0; 
    y0 = 1.2;
    
    % Välj steglängd från n
    n0 = n_e(jj);
    h = (T_e-t0)/n0;
    
    % Initiera vektorer som sparar punkter för t och y
    tvec = zeros(n0+1, 1);
    yvec = zeros(n0+1, 1);

    % Sätt in startpunkt
    tvec(1) = t0; 
    yvec(1) = y0;

    % printa vilken iteration det är
    % fprintf('\n| Börjar Euler-framåt med steglängd %d |\n', n0);
    
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
        % disp([ii y y0 t])
    end

    % Skapa subplot för aktuell iteration
    subplot(2, 2, jj);
    
    % Plotta analytisk lösning
    t_analytic_vals = linspace(0, T_e, 1000);
    plot(t_analytic_vals, y_analytic(t_analytic_vals), 'Color', 'Black', 'LineWidth', 2); hold on;

    % Plotta numerisk lösning
    plot(tvec, yvec, 'o-', 'Color', colors(jj,:), 'LineWidth', 1.5, 'MarkerSize', 4);
    
    % Lägg till etiketter och titel
    xlabel('t');
    ylabel('y(t)');
    title(sprintf('Euler framåt, n = %d', n0));
    grid on;
    legend('Analytisk lösning', sprintf('Euler-framåt (n=%d)', n0), 'Location', 'best');
    hold off;
    
end

% --------------------------------------------------
% Uppgift 2f
% Euler bakåt med [0, 80] och n = [50 100 400 800]

T_f = 80;
n_f = [50 100 400 800]; % antal steg h = 8/n

% Genererar 4 st distinkta färger för olika steglängd i n
colors = lines(length(n_f));  

% Iterera över varje steglängd i n
figure;
for jj = 1:length(n_f)
    % Resetta och initiera startvärden 
    t0 = 0; 
    y0 = 1.2;
    
    % Välj steglängd från n
    n0 = n_f(jj);
    h = (T_f-t0)/n0;
    
    % Initiera vektorer som sparar punkter för t och y
    tvec = zeros(n0+1, 1);
    yvec = zeros(n0+1, 1);

    % Sätt in startpunkt
    tvec(1) = t0; 
    yvec(1) = y0;

    % printa vilken iteration det är
    % fprintf('\n| Börjar Euler-framåt med steglängd %d |\n', n0);
    
    % inre for-loop for Euler-fram
    for ii = 1:n0
         % Nya t
        t = t0 + h;

        % y(i+1) = y(i) + h*sin(3t(i+1)) - h*2*y(i+1)
        % y(i+1) + 2h*y(i+1) = y(i) + h*sin(3t(i+1))
        % y(i+1)(1+2h) = y(i) + h*sin(3t(i+1))
        % y(i+1) = y(i) + h*sin(3t(i+1)) / (1+2h)

        % Uppdatera y och  med Euler-bakåt
        y = (y0 + h*sin(3*t)) / (1+2*h); 
        
        % Lägg till nya värden i vektorn for t och y punkter
        tvec(ii+1) = t;
        yvec(ii+1) = y;
        
        % Uppdatera y0 och t0
        y0 = y;
        t0 = t;

        % Printa ut resultaten
        % disp([ii y y0 t])
    end

    % Skapa subplot för aktuell iteration
    subplot(2, 2, jj);
    
    % Plotta analytisk lösning
    t_analytic_vals = linspace(0, T_f, 1000);
    plot(t_analytic_vals, y_analytic(t_analytic_vals), 'Color', 'Black', 'LineWidth', 2); hold on;

    % Plotta numerisk lösning
    plot(tvec, yvec, 'o-', 'Color', colors(jj,:), 'LineWidth', 1.5, 'MarkerSize', 4);
    
    % Lägg till etiketter och titel
    xlabel('t');
    ylabel('y(t)');
    title(sprintf('Euler bakåt, n = %d', n0));
    grid on;
    legend('Analytisk lösning', sprintf('Euler-bakåt (n=%d)', n0), 'Location', 'best');
end

% SLUTSATS: Euler-bakåt fungerar bättre för större steglängd h och är mer
% stabil