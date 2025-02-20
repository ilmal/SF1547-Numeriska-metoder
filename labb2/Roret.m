a = 1; k = 1; Ti = 450; Te = 20; r0 = 1; rslut = 2; 
N = 25;
h = (rslut-r0)/(N+1);

Tresult = rand(a, N, k, Ti, Te, r0, rslut);
Tend = (a*h*Te*Tresult(end))/(k+a*h);
disp(['Temperatur vid r = 2: ', num2str(Tend)]);

function Tnew = rand(a, N, k, Ti, Te, r0, rslut)
    h = (rslut-r0)/(N+1); % steglängd
    ri = (r0+h:h:rslut-h)'; % inre punkter 

    A = zeros(N, N); % fylld med nollor
    A(1,1:2) = [(-2*ri(1))/h^2, (ri(1)/h^2)+(1/(2*h))]; % första raden på A

    % inre punkter
    for ii = 2:N-1
        A(ii, ii-1) = (ri(ii)/h^2)-(1/(2*h)); % för T(i-1)
        A(ii, ii) = (-2*ri(ii))/h^2; % T(i)
        A(ii, ii+1) = (ri(ii)/h^2)+(1/(2*h)); % T(i+1)
    end

    % Sista raden utanför for-loop, rad N
    A(N, N-1:N) = [(ri(end)/h^2)-(1/(2*h)), (-2*ri(N))/h^2 + (1/(k+a*h))*(ri(N)/h^2+(1/(2*h)))];

    A = sparse(A);
    b = zeros(N, 1);
    b(1) = -Ti * (ri(1)/h^2 - 1/(2*h));
    b(end) = ((-a*h*Te)/(k+a*h))*(ri(end)/h^2+(1/(2*h)));

    T = A\b;
    Tnew = T;

end


