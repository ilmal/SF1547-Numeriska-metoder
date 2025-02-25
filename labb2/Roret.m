format long; 

% --------------------------------------------------
% Uppgift 4c

% givna värden
a = 1; k = 1; Ti = 450; Te = 20; r0 = 1; rslut = 2; 
N = 25; % ska vara N+1 delintervall

tol = 10^(-2.5); % tolerans för att garantera första decimalen korrekt
diff = 0.1;    % godtyckligt startvärde på diff
prev_T_end = 0; % godtyckligt startvärde på tidigare T(end)


fprintf('----------------------------------------------------------\n')
fprintf('Uppgift 4c: Tend där N fördubblas varje iteration\n \n')

% while-loop tills given tolerans är uppnådd
while diff > tol
    h = (rslut-r0)/(N+1); % steglängd
    
    % Ax=b behandlar endast inre punkter och ej randpunkter
    A = zeros(N, N); % skapa default 0-matris med dimensioner NxN
    b = zeros(N, 1); % skapa vektorn b i ekvationen Ax = b
    
    ri = (r0+h:h:rslut-h)'; % alla inre punkter r med steglängd h
    
    % sätt kolumn 1-2 i rad 1 (rad 1 behandlas separat)
    A(1, 1:2) = [-2*ri(1)/h^2, (ri(1)/h^2)+(1/(2*h))]; % beräknade värden
    b(1) = -Ti*((ri(1)/h^2)-(1/(2*h))); % beräknade värden

    % for-loop genom inre rader med ihopsamlade termer
    for ii=2:N-1
        % sätt rad ii och kolumn ii-1
        A(ii,ii-1) = (ri(ii)/h^2)-(1/(2*h));
        
        % sätt rad ii och kolumn ii
        A(ii, ii) = -2*ri(ii)/h^2;
    
        % sätt rad ii och kolumn ii+1
        A(ii, ii+1) = (ri(ii)/h^2)+(1/(2*h));
    end
    
    % sista raden N behandlas separat, sätt sista två kolumnerna N-1 och N
    A(N, N-1:N) = [(ri(end)/h^2)-(1/(2*h)), (-2*ri(end)/h^2)+(k/(k+a*h))*((ri(end)/h^2)+(1/(2*h)))];
    b(N) = -a*h*Te/(k+a*h)*((ri(end)/h^2)+(1/(2*h)));
    
    A = sparse(A); % lagra som gles matris
    
    T = A\b; % lös systemet och spara i T

    Tend = (a*h*Te+k*T(end))/(k+a*h); % beräkna T(N+1), randpunkten
    
    fprintf('Iteration with N = %d: T(N+1) = %.10f\n', N, Tend);

    % beräkna diff absolutbelopp på sista punkten Tend
    diff = abs(Tend-prev_T_end);
    prev_T_end = Tend; % uppdatera senaste Tend

    fprintf('Difference at the moment: %.10f\n \n', diff);
        
    % dubbla N
    N = 2 * N;
end

fprintf('Final T(end) = %.10f\n', Tend);

% --------------------------------------------------
% Uppgift 4d

k = 1; N_a = 20; % värde på k samt antal delintervall som N_a
a0 = 0; aslut = 10; % alpha intervall
N = 10000; % värde på N för beräkning

h_a = (aslut-a0)/N_a; % steglängd för alpha

ai = a0:h_a:aslut; % skapa alla intervaller för alpha
Tres = zeros(1, length(ai)); % initiera vektor för y-axel

fprintf('----------------------------------------------------------\n')
fprintf('Uppgift 4d: Tend baserat på alpha = [0, 10]\n \n')

% for-loop genom alla värden på alpha med steglängd h_a
for jj = 1:length(ai)
    h = (rslut-r0)/(N+1); % steglängd för N
    
    % Ax=b behandlar endast inre punkter och ej randpunkter
    A = zeros(N, N); % skapa default 0-matris med dimensioner NxN
    b = zeros(N, 1); % skapa vektorn b i ekvationen Ax = b
    
    ri = (r0+h:h:rslut-h)'; % alla inre punkter r med steglängd h
    
    % sätt kolumn 1-2 i rad 1 (rad 1 behandlas separat)
    A(1, 1:2) = [-2*ri(1)/h^2, (ri(1)/h^2)+(1/(2*h))]; % beräknade värden
    b(1) = -Ti*((ri(1)/h^2)-(1/(2*h))); % beräknade värden

    % for-loop genom inre rader med ihopsamlade termer
    for ii=2:N-1
        % sätt rad ii och kolumn ii-1
        A(ii,ii-1) = (ri(ii)/h^2)-(1/(2*h));
        
        % sätt rad ii och kolumn ii
        A(ii, ii) = -2*ri(ii)/h^2;
    
        % sätt rad ii och kolumn ii+1
        A(ii, ii+1) = (ri(ii)/h^2)+(1/(2*h));
    end
    
    % sista raden N behandlas separat, sätt sista två kolumnerna N-1 och N
    A(N, N-1:N) = [(ri(end)/h^2)-(1/(2*h)), (-2*ri(end)/h^2)+(k/(k+ai(jj)*h))*((ri(end)/h^2)+(1/(2*h)))];
    b(N) = -ai(jj)*h*Te/(k+ai(jj)*h)*((ri(end)/h^2)+(1/(2*h)));
    
    A = sparse(A); % lagra som gles matris
    
    T = A\b; % lös systemet och spara i T

    Tend = (a*h*Te+k*T(end))/(k+ai(jj)*h); % beräkna T(N+1), randpunkten
    Tres(jj) = Tend; % lägg till Tend till Tres för att plotta
    
    % printa resultatet på Tend för nuvarande alpha-värde
    fprintf('Iteration med a = %.1f: T(N+1) = %.6f\n', ai(jj), Tend);
end

% plotta
figure;
plot(ai, Tres, '-o', LineWidth=1.5);
legend('Temperaturen som funktion av alpha');
xlabel('Alpha'); ylabel('Temperatur');