format long;

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

for ii = 1:length(n)
    % Beräkna steglängd baserat på värde från n
    n0 = n(ii);                     % välj antal punkter
    h = (endVal - startVal) / n0;   % beräkna steglängd
    sum = 0;                        % sätt startsumma till 0

    % Beräkna summan av alla termer gånger h
    for jj = 1:n0
        x = startVal + h*(jj-1); % uppdatera till nya punkten
        term = y_2(x);           % hitta termen
        
        % Om det är första eller sista term, dela med två
        if jj == 1||jj == n0
            term = term/2;
        end
        
        % Uppdatera summan
        sum = sum+term;
        
    end
    
    % Hitta res genom att multiplicera med h och pi
    res = sum * h * pi;
    fprintf('\n| %6d     |    %20.6f  |', n0, res) % printa resultat

end

fprintf('\n-----------------------------------------\n')

% -----------------------------

% 2. Simpsons metod
%    Approximera med andragradspolynom och steglängden h

fprintf('\n-----------------------------------------')
fprintf('\n Approximation using Simpsons Method')
fprintf('\n-----------------------------------------')
fprintf('\n|     n      |       Result (V)         |')
fprintf('\n|------------|--------------------------|')

for ii = 1:length(n) 
    n0 = n(ii);                 % välj antal punkter
    h = (endVal-startVal)/n0;   % beräkna steglängd
    sum = 0;                    % sätt startsumma till 0
    
    % loopa till n0+1 eftersom n0 är jämnt och vi vill ha jämnt antal
    % delintervall -> n0 jämnt gör udda antal delintervall
    for jj = 1:n0+1
        x = startVal + h*(jj-1);
        term = y_2(x);
        
        if jj == 1 || jj == n0+1  % Första och sista term
            term = term;          % Ingen multiplikation (ska vara 1)
        elseif mod(jj, 2) == 0  % Jämn term
            term = 4 * term;
        elseif mod(jj, 2) == 1  % Udda term
            term = 2 * term;
        end

        sum = sum+term;

    end
    
    res = (h/3) * sum * pi;
    
    fprintf('\n| %6d     |    %20.6f  |', n0, res)
end

fprintf('\n-----------------------------------------\n')
