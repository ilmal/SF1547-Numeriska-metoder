% -----------------
% Uppgift 2a
% Newtons metod: x(n+1) = x(n) - (f(xn)/f'(xn)), n = 1, 2, ...

% så att alla print har tillräckligt många decimaler
format long

% värde på H 
H = 0.5;
% ------
% för Uppgift 2c 
% initialisera tomma arrays för att plotta konvergens
errors_newton = [];
errors_secant = [];
% ------

% anonym modell där H = y = 0.5
% 0.5 = 8e^(-t/2)*cos(3t) -> f(x) = 0 ger 0 = 8e^(-t/2)*cos(3t) - 0.5
y = @(t) 8.*(exp(1)).^(-t/2).*cos(3.*t) - H;
y_prim = @(t) -24.*(exp(1)).^(-t/2).*sin(3.*t) - 4.*exp(1).^(-t/2).*cos(3.*t);

% plotta funktion för att gissa startvärde
% t = linspace(1, 100);
% y_plot = y(t);
% figure;
% plot(t, y_plot);


% startvärde baserad på plottad graf och tolerans från uppgiften
start_x = 4.5;
tolerance = 10^(-8);

% skillnaden mellan nya och gamla x, antal iterationer, och max iterationer
difference_x = 1; 
iterations = 0; 
max_iterations = 1000;

disp('Newtons metod:')

% iterering med newtons metod
while difference_x > tolerance && iterations < max_iterations
    iterations = iterations+1; % öka iterationscounter
    new_x = start_x - y(start_x)/y_prim(start_x); % beräkna nya x med newtonmetoden
    difference_x = abs(new_x - start_x); % beräkna skillnad mellan nya och gamla x
    errors_newton = [errors_newton, difference_x]; % uppg 2c
    start_x = new_x; % sätt nya x som start x
    disp([iterations new_x difference_x]) % display
end

% -----------------
% Uppgift 2b
% Samma H men med sekantmetoden
start1_x = 4.4; % x0 - startgissning 1
start2_x = 4.6; % x1 - startgissning 2
difference_sekant_x = 1; % skillnaden
iterations_sekant = 0; % antalet iterationer

disp('Sekantmetoden:')

% iterera med sekantmetoden
while difference_sekant_x > tolerance && iterations_sekant < max_iterations
    iterations_sekant = iterations_sekant+1; % öka itereringscounter
    new_x_sekant = start2_x - y(start2_x)*((start2_x - start1_x)/(y(start2_x) - y(start1_x))); % beräkna nytt x med sekantmetoden
    difference_sekant_x = abs(new_x_sekant - start2_x); % beräkna felet
    errors_secant = [errors_secant, difference_sekant_x]; % uppg 2c
    start1_x = start2_x; % uppdatera start1 först
    start2_x = new_x_sekant; % uppdatera start2 sen
    disp([iterations_sekant new_x_sekant start2_x difference_sekant_x]) % display
end

% Newtons metod löser problemet på 3 iterationer medan sekantmetoden tar 5 iterationer
% Anledningen till detta är att Netwtons metod har en 
% konvergenshastighet som är kvadratisk till skillnad från den superlinjära hastigheten hos sekantmetoden. 
% Detta leder till att Newtons metod kommer vara snabbare i detta fall då funktionen är snäll (deriverbar)



% -----------------
% Uppgift 2c
% Plotta vardera metods konvergenshastighet
figure; % ny figur
semilogy(1:length(errors_newton), errors_newton, 'b-o', 'DisplayName', 'Newton'); % newton
hold on;
semilogy(1:length(errors_secant), errors_secant, 'r-o', 'DisplayName', 'Secant'); % sekant
xlabel('Iteration');
ylabel('|T(n+1) - T(n)|');

% -----------------
% Uppgift 2d
H2 = 2.8464405473; % nytt H-värde
y_2d = @(t) 8.*(exp(1)).^(-t/2).*cos(3*t) - H2; % samma modell med nytt H
y_2d_prim = @(t) -24*(exp(1)).^(-t/2)*sin(3*t) - 4*exp(1).^(-t/2)*cos(3*t); % derivatan

% startgissning baserad på plottad graf och tolerans från uppgift
start2d_x = 2.1;
tolerance_2d = 10^(-8);
errors_2d = []; % för att plotta konvergens

% antal iterationer + max antal iterationer
iterations_2d = 0;
max_iterations_2d = 1000;
difference_2d = 1;

disp("Newton med nytt H-värde: ");

% iterera med newtonmetoden
while iterations_2d < max_iterations_2d && difference_2d > tolerance_2d
    iterations_2d = iterations_2d+1; % öka iterationscounter
    new_x_2d = start2d_x - (y_2d(start2d_x)/y_2d_prim(start2d_x)); % beräkna nytt x med newtonmetoden
    difference_2d = abs(new_x_2d - start2d_x); % beräkna differens
    errors_2d = [errors_2d difference_2d]; % för att plotta konvergens
    start2d_x = new_x_2d; % sätt nya x som start
    disp([iterations_2d new_x_2d start2d_x difference_2d]); % display
end

% plotta figur med nytt H gentemot gammalt H
figure;
semilogy(1:length(errors_2d), errors_2d, 'b-o', 'DisplayName', 'Newton');
hold on;
semilogy(1:length(errors_newton), errors_newton, 'r-o', 'DisplayName', 'Newton')
xlabel('Iteration');
ylabel('|T(n+1) - T(n)|');

% Anledningen till att newtons metod går långsammare är eftersom det finns
% en rot med högre multiplicitet än 1 nära nollpunkten vi söker. I detta 
% fall kommer metoden att konvergera linjärt istället för kvadratiskt. 

% -----------------