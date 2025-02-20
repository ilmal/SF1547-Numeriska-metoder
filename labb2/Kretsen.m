format long;

% Uppgift 3a
% -----------------------------
% Spänningen i kretsen ges av: 
% Lq''(t) + Rq'(t) + (1/C)q(t) = 0 
% q(0) = 1 
% dq/dt (0) = 0

% Strömmen ges av: 
% i = dq/dt

% Andraderivatan skrivs om till: q'' = -((1/C)q - Rq')/L

% Omskrivning:
% q'' = f(t, q, q') 
% -> Q(t) = q(t), I(t) = q'(t)
% -> Q' = q' = I, I' = q'' = f(t, q, q') = -((1/C)Q - RI)/L

U = [Q; I];
F = [   I; 
    -((1/C)*Q + R*I)/L]; % F(t, U)

U0 = [1; 0]; % U(t0) = [c1, c2]

% -----------------------------
% Uppgift 3b
% Använd ode45

function F = returnVectorF(t, y, R, L, C)
    F = [y(2);
        -((1/C)*y(1) + R*y(2))/L
    ];
    
end




