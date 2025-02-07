% Differentialekvation: 
%                       dy/dt = sin(3t) - 2y 
%                       y(0) = 1.2, t ∈ [0,8]
% --------------------------------------------------
% Uppgift 2a
% Analytisk lösning:
% y(t) = (93/65)e^(-2t)-(3/13)cos(3t)+(2/13)sin(3t)

% --------------------------------------------------
% Uppgift 2b

t = 0; T = 8;
n = 50; % antal steg h = 8/n
h = (T-t)/n;
y0 = 1.2;

tvec = zeros(n+1, n); yvec = zeros(n+1, n);
tvec(1) = t; yvec(1) = y0;

% Euler fram: y(i+1) = y(i) + hf(t(i), y(i))
f = @(t) (93/65)*exp(-2*t)-(3/13)*cos(3*t)+(2/13)*sin(3*t);

for ii = 1:n
    y = y + h*f(y);
    t = t + h;
    tvec(ii+1) = t;
    yvec(ii+1) = y;
    y0 = y;
    disp([ii y y0 t])
end


