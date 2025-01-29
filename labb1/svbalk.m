

e = exp(1);
% formel
y = @(t) 8.*e.^(-t./2).*cos(3.*t);
y_prim = @(t) -4.*e.^(-t./2).*cos(3.*t) - 24.*e.^(-t./2)*sin(3.*t);

H = 0.5;

x = 4.5;
tol = 0.00000001;

diffx = 1; iter = 0; maxiter = 100;

disp(["iteration" "time" "diff" "value"])

while diffx > tol && iter < maxiter
    iter  = iter + 1;
    xnew = x - ( y(x) - H ) ./ y_prim(x);
    diffx = abs(xnew-x);
    x = xnew;
    disp([iter xnew diffx y(x)])
end


