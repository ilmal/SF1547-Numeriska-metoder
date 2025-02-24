format long; 

% Nytt försök
% --------------------------------------------------
% Uppgift 4c

% givna värden
a = 1; k = 1; Ti = 450; Te = 20; r0 = 1; rslut = 2; 
N = 25; % N+1 delintervall

tol = 10^(-2);
diff = 0.1;
prev_T_end = 0;

while diff > tol
    h = (rslut-r0)/(N+1); % steglängd

    A = zeros(N, N); % skapa default 0-matris med rätt dimensioner
    b = zeros(N, 1); % skapa vektorn b i ekvationen Ax = b
    
    ri = (r0:h:rslut)'; % alla punkter r med steglängd h
    
    % sätt kolumn 1-2 i rad 1 (rad 1 behandlas separat)
    A(1, 1:2) = [-2*ri(1)/h^2, (ri(1)/h^2)+(1/(2*h))];
    b(1) = -Ti*((ri(1)/h^2)-(1/(2*h)));

    % for-loop genom inre rader
    for ii=2:N-1
        % sätt rad ii och kolumn ii-1
        A(ii,ii-1) = (ri(ii)/h^2)-(1/(2*h));
        
        % sätt rad ii och kolumn ii
        A(ii, ii) = -2*ri(ii)/h^2;
    
        % sätt rad ii och kolumn ii+1
        A(ii, ii+1) = (ri(ii)/h^2)+(1/(2*h));
    end

    % sista raden N, sätt sista två kolumnerna N-1 och N
    A(N, N-1:N) = [(ri(end)/h^2)-(1/(2*h)), (-2*ri(end)/h^2)+(1/(k+a*h))*((ri(end)/h^2)+(1/(2*h)))];
    b(N) = -a*h*Te/(k+a*h)*((ri(end)/h^2)+(1/(2*h)));
    
    A = sparse(A); % lagra som gles matris
    
    T = A\b;
    
    fprintf('Iteration with N = %d: T(end) = %.10f\n', N, T(end));

    diff = abs(T(end)-prev_T_end);
    prev_T_end = T(end);
    
    fprintf('Difference at the moment: %.10f\n', diff);
        
    % dubbla N
    N = 2 * N;
end

fprintf('\nFinal T(end) = %.10f\n', T(end));




