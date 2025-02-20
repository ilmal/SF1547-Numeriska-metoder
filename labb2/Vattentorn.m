format long;

% Uppgift 1a
% -----------------------------

% Trapetsmetoden

beta = 0.2; % givet värde på beta
y = @(x) (exp(beta*x)+8)./(1+(x/5).^3); % y
y_2 = @(x) ((exp(beta*x)+8)./(1+(x/5).^3)).^2; % y^2

a = 0; b = 20; % Integrationsgränser
N = [100 500 800 1280 2000]; % Antalet steg

errors = zeros(size(N));

I = integral(y_2, a, b); % Faktiskt integralvärde att jämföra med

for ii = 1:length(N)
    h = (b-a)/N(ii);

    x = a:h:b; % punkterna som ska räknas ut
    yx = y_2(x);  % räkna ut funktionen i x-värdena

    Th = h*(sum(yx) - 0.5*(yx(1)+yx(end))); % trapetsregeln
    res = pi*Th; % multiplicera resultat med pi
    
    errors(ii) = abs(res-pi*I);

    fprintf('N = %d, eh = %.6e\n', N(ii), errors(ii));
end

% Beräkna noggrannhetsordningen p
fprintf('\n');
for ii = 1:length(N)-1
    % två intilliggande h-värden
    h1 = (b-a)/N(ii);
    h2 = (b-a)/N(ii+1);

    % två intilliggande fel
    e1 = errors(ii);
    e2 = errors(ii+1);
    
    % eh ≈ Ch^p -> log eh = log C + p * log h
    % två fel: e1/e2 = (h1/h2)^p -> log (e1/e2) = log (h1/h2) * p
    % alltså:
    p = log(e1/e2) / log(h1/h2); % beräkna noggrannhetsordning p

    fprintf('p mellan N=%d och N=%d: %.4f\n', N(ii), N(ii+1), p);
end


% Simpsons metod

