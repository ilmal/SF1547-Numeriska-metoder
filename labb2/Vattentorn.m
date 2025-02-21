format long;

% Uppgift 1a
% -----------------------------

beta = 0.2; % givet värde på beta
y = @(x) (exp(beta*x)+8)./(1+(x/5).^3); % y
y_2 = @(x) ((exp(beta*x)+8)./(1+(x/5).^3)).^2; % y^2

a = 0; b = 20; % Integrationsgränser
N = [100 500 800 1280 2000]; % Några olika antal steg

I = integral(y_2, a, b); % Faktiskt integralvärde att jämföra med

% TRAPETSMETODEN

errors_T = zeros(size(N)); % spara fel för varje N-värde

fprintf('\n-------------------------------------');
fprintf('\nTrapetsmetoden\n');

% Trapetsmetoden för varje värde N i for-loop
for ii = 1:length(N)
    h = (b-a)/N(ii);

    x = a:h:b; % punkterna som ska räknas ut
    yx = y_2(x);  % räkna ut funktionen i x-värdena

    Th = h*(sum(yx) - 0.5*(yx(1)+yx(end))); % trapetsregeln
    res = pi*Th; % multiplicera resultat med pi
    
    errors_T(ii) = abs(res-pi*I); % Spara felet på index ii

    fprintf('N = %d,  res = %.6f, eh = %.6e\n', N(ii), res, errors_T(ii));
end


fprintf('\n');

% Beräkna noggrannhetsordningen p
for ii = 1:length(N)-1

    % två intilliggande h-värden
    h1 = (b-a)/N(ii);
    h2 = (b-a)/N(ii+1);

    % två intilliggande fel
    e1 = errors_T(ii);
    e2 = errors_T(ii+1);
    
    % eh ≈ Ch^p -> log eh = log C + p * log h
    % två fel: e1/e2 = (h1/h2)^p -> log (e1/e2) = log (h1/h2) * p
    % alltså:
    p = log(e1/e2) / log(h1/h2); % beräkna noggrannhetsordning p

    fprintf('p mellan N=%d och N=%d: %.4f\n', N(ii), N(ii+1), p);
end
fprintf('-------------------------------------\n');


% SIMPSONS METOD
fprintf('Simpsons metod\n');

errors_S = zeros(size(N)); % spara felet för varje N-värde 

% Simpsons metod för varje värde N i for-loop
for ii = 1:length(N)
    h = (b-a)/N(ii); % beräkna steglängd för nuvarande N

    x = a:h:b; % punkterna som ska räknas ut
    yx = y_2(x);  % räkna ut funktionen i x-värdena

    Sh = h/3 * (yx(1) + 4*sum(yx(2:2:end-1)) + 2*sum(yx(3:2:end-2)) + yx(end));
    res = pi*Sh; % multiplicera resultat med pi
    
    errors_S(ii) = abs(res-pi*I); % Spara felet på index ii

    fprintf('N = %d, res = %.6f, eh = %.6e\n', N(ii), res, errors_S(ii));
end

fprintf('\n');

% Beräkna noggrannhetsordningen p
for ii = 1:length(N)-1

    % två intilliggande h-värden
    h1 = (b-a)/N(ii);
    h2 = (b-a)/N(ii+1);

    % två intilliggande fel
    e1 = errors_S(ii);
    e2 = errors_S(ii+1);
    
    % eh ≈ Ch^p -> log eh = log C + p * log h
    % två fel: e1/e2 = (h1/h2)^p -> log (e1/e2) = log (h1/h2) * p
    % alltså:
    p = log(e1/e2) / log(h1/h2); % beräkna noggrannhetsordning p

    fprintf('p mellan N=%d och N=%d: %.4f\n', N(ii), N(ii+1), p);
end

fprintf('-------------------------------------\n');