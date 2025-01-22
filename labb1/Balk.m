% Uppgift 2a
% Newtons metod: x(n+1) = x(n) - (f(xn)/f'(xn)), n = 1, 2, ...

format long

% startvärden
H = 0.5;

% anonym modell: H = y = 0.5
% 0.5 = 8e^(-t/2)*cos(3t) -> f(x) = 0 ger 0 = 8e^(-t/2)*cos(3t) - 0.5
y = @(t) 8*(exp(1)).^(-t/2)*cos(3*t) - H;
y_prim = @(t) -24*(exp(1)).^(-t/2)*sin(3*t) - 4*exp(1).^(-t/2)*cos(3*t);

% startvärde baserad på plot och tolerans från uppg.
start_x = 4.5;
tolerance = 10^(-8);

% skillnaden mellan nya och gamla x, antal iterationer, och max iterationer
difference_x = 1; 
iterations = 0; 
max_iterations = 1000;

disp('Newtons metod:')

% iterering med newtons metod
while difference_x > tolerance && iterations < max_iterations
    iterations = iterations+1;
    new_x = start_x - y(start_x)/y_prim(start_x);
    difference_x = abs(new_x - start_x);
    start_x = new_x;
    disp([iterations new_x difference_x])
end

% Uppgift 2b
% Samma H men med sekantmetoden
start1_x = 4.4; % x0
start2_x = 4.6; % x1
difference_sekant_x = 1;
iterations_sekant = 0;

disp('Sekantmetoden:')

while difference_sekant_x > tolerance && iterations_sekant < max_iterations
    iterations_sekant = iterations_sekant+1;
    new_x_sekant = start2_x - y(start2_x)*((start2_x - start1_x)/(y(start2_x) - y(start1_x)));
    difference_sekant_x = abs(new_x_sekant - start2_x);
    start1_x = start2_x; % uppdatera start1 först
    start2_x = new_x_sekant;
    disp([iterations_sekant new_x_sekant start2_x difference_sekant_x])
end

% Uppgift 2c
% Plotta vardera metods konvergenshastighet

