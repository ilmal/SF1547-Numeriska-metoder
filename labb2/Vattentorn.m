% Uppgift 1a
% -----------------------------
% Fördefinierade värden & funktion
% 0.1 ≤ beta ≤ 0.3
% 0 ≤ x ≤ 20
% V = π integral(0, 20, y(x;beta)^2)dx

beta = 0.2; % givet värde på beta
y = @(x) (exp(beta*x)+8)/(1+(x/5)^3); % y
y_2 = @(x) ((exp(beta*x)+8)./(1+(x/5).^3)).^2; % y^2

% Integrationsgränser och olika värden på n
startVal = 0; endVal = 20;
n = [10 20 40 80 100 200 500 1000 1280];

% Exakt värde för felanalys
I = integral(y_2, startVal, endVal);

% 1. Trapetsmetoden
%    Approximera med räta linjer och steglängden h

fprintf('\n-----------------------------------------')
fprintf('\n Approximation using Trapezoidal Method')
fprintf('\n-----------------------------------------')
fprintf('\n|     n      |       Result (V)         |')
fprintf('\n|------------|--------------------------|')

% För felanalys
errors_trapeze = []; 
errors_trapeze_h = []; 

for ii = 1:length(n)
    % Beräkna steglängd baserat på värde från n
    n0 = n(ii);
    h = (endVal-startVal) / n0;
    errors_trapeze_h = [errors_trapeze_h, h]; 
    sum = 0; % sätt startsumma till 0

    % Beräkna summan av alla termer gånger h
    for jj = 1:n0
        x = startVal + h*(jj-1); % Uppdatera till nya punkten
        term = y_2(x); % hitta termen
        
        % Om det är första eller sista term, dela med två
        if jj == 1||jj == n0+1
            term = term/2;
        end
        
        % Uppdatera summan
        sum = sum+term;
        
    end
    
    % Hitta res genom att multiplicera med h och pi
    res = sum * h * pi;
    
    % Lägg till felet för n0 i en array
    error = abs(pi*I - res);
    errors_trapeze = [errors_trapeze, error];

    fprintf('\n| %6d     |    %20.6f  |', n0, res)

end

fprintf('\n-----------------------------------------\n')

% Felanalys med figur
% eh ≤ Ch^2 -> ln eh ≤ ln C + 2 * ln h
figure;
loglog(errors_trapeze_h, errors_trapeze, '-o', 'LineWidth', 1.5);
xlabel('Steglängd h');
ylabel('Felet |I - Th|');
title('Konvergens av Trapetsmetoden');
grid on;

p1 = polyfit(log(errors_trapeze_h), log(errors_trapeze), 1);
disp(['Noggrannhetsordning Trapetsmetoden, p ≈ ', num2str(p1(1))]);

% -----------------------------

% 2. Simpsons metod
%    Approximera med andragradspolynom och steglängden h

fprintf('\n-----------------------------------------')
fprintf('\n Approximation using Simpsons Method')
fprintf('\n-----------------------------------------')
fprintf('\n|     n      |       Result (V)         |')
fprintf('\n|------------|--------------------------|')

errors_simpson = []; % För felanalys
errors_simpson_h = [];

for ii = 1:length(n) 
    n0 = n(ii);
    h = (endVal-startVal)/n0;
    errors_simpson_h = [errors_simpson_h, h];
    sum = 0;

    for jj = 1:n0
        x = startVal + h*(jj-1);
        term = y_2(x);
        
        % Om det är en jämn term OCH inte första termen
        if mod(jj, 2) == 0 && jj ~= 1
            term = 4*term;
        end
        
        % Om det är udda term OCH inte sista termen
        if mod(jj, 2) == 1 && jj~=n0
            term = 2*term;
        end

        sum = sum+term;

    end
    
    res = (h/3) * sum * pi;
    
    error = abs(I - (h/3)*sum);
    errors_simpson = [errors_simpson, error];

    fprintf('\n| %6d     |    %20.6f  |', n0, res)
end

fprintf('\n-----------------------------------------\n')

% Felanalys med figur
% eh ≤ Ch^2 -> ln eh ≤ ln C + 2 * ln h
figure;
loglog(errors_simpson_h, errors_simpson, '-o', 'LineWidth', 1.5);
xlabel('Steglängd h');
ylabel('Felet |I - Th|');
title('Konvergens av Simpsons metod');
grid on;

p1 = polyfit(log(errors_trapeze_h), log(errors_trapeze), 1);
disp(['Noggrannhetsordning Simpsons metod, p ≈ ', num2str(p1(1))]);