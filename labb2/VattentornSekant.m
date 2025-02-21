format long;

% Uppgift 1b
% -----------------------------
% Hitta ett beta så att volymen blir 1500m^3

% F(beta) = pi * integral(0, 20, y(x; beta)^2)dx - 1500 = 0

y = @(x, beta) (exp(beta*x)+8)./(1+(x/5).^3); % y
y_2 = @(x, beta) ((exp(beta*x)+8)./(1+(x/5).^3)).^2; % y^2

a = 0; b = 20; % integrationsgränser
f = @(beta) pi*integral(@(x) y_2(x, beta), a, b) - 1500;

% Plotta för att hitta startgissningar
% figure;
% x = linspace(0.1, 0.3, 100); % 0.1 ≤ beta ≤ 0.3
% fx = arrayfun(f, x);
% plot(x, fx, "LineWidth", 1.5); hold on;
% yline(0, '--r', 'y=0', 'LineWidth', 1.5);

x0 = 0.262; x1 = 0.263; % startgissningar baserat på plot
tol = 10^(-12); % tolerans

iter = 0; maxiter = 100; % iterationer

errors = zeros(1, maxiter); % errors 

fprintf("\nIter |      x1       |       x0       |       p  ")
fprintf("\n----------------------------------------------------")
fprintf("\n%d    |  %.10f |  %.10f  |  ----", iter, x1, x0);

% Sekantmetoden
while abs(x1-x0) > tol && iter < maxiter 
    iter = iter+1; % uppdatera iter
    x = x1 - f(x1) * (x1-x0) / (f(x1)-f(x0)); % sekantmetoden

    errors(iter) = abs(x-x1); % spara errors
    
    % beräkna konvergensen om villkor stämmer
    if iter > 2
        p = log(errors(iter)-errors(iter-1)) / log(errors(iter-1)-errors(iter-2));
        fprintf("\n%d    |  %.10f |  %.10f  |  %.4f", iter, x1, x0, p);
    % printa ---- om inget p-värde kan beräknas    
    else
        fprintf("\n%d    |  %.10f |  %.10f  |  ---- ", iter, x1, x0);
    end
    
    % uppdatera x0 och x1
    x0 = x1;
    x1 = x;
    
end

% Sätt x1 som resultat och printa
res = x1;
fprintf("\n \nResultat = %.10f\n", res);