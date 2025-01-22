load STHLMTEMP.mat

% UPPGIFT 1a
% värde av k
k = 2*pi/365;

% anonym funktion T(t)
T = @(c, t) c(1) + c(2)*sin(k*t) + c(3)*cos(k*t) + c(4)*sin(2*k*t) + c(5)*cos(2*k*t);

% t = tiden från dag 1 till dag 98251
% Tdm = b -> temp i grader
t = (1:98251)';
A = [ones(size(t)) sin(k*t) cos(k*t) sin(2*k*t) cos(2*k*t)];

% \ löser med A' automatiskt
c = A\Tdm;

% printa ut resultat
disp("Koefficientvärden c: ")
fprintf("c0: %d\n", c(1));
fprintf("c1: %d\n", c(2));
fprintf("c2: %d\n", c(3));
fprintf("c3: %d\n", c(4));
fprintf("c4: %d\n", c(5));

% plotta datan i figur 1
figure(1);
x = linspace(1, 98251, 98251);
p = T(c, x);
plot(x, p), hold on
plot(x, Tdm, 'o')
xlabel('t, tiden i dygn')
ylabel('T, temperaturen i grader')

% UPPGIFT 1b
% residualvektor
r = Tdm - (p)';

% plotta residual mot tid i figur 2
figure(2);
plot(x, r)
xlabel('t, tiden i dygn')
ylabel('y, residual')

% minstakvadratsumma
s = sum(r.^2);
fprintf("\nMinstakvadratsumman: %d\n\n", s);

% UPPGIFT 1c
new_T = @(a, t) a(1) + a(2)*t + a(3)*t.^2 + a(4)*sin(k*t) + a(5)*cos(k*t) + a(6)*sin(2*k*t) + a(7)*cos(2*k*t);
new_A = [ones(size(t)) t t.^2 sin(k*t) cos(k*t) sin(2*k*t) cos(2*k*t)];

a = new_A\Tdm;

disp("Koefficientvärden a: ")
fprintf("a0: %d\n", a(1));
fprintf("a1: %d\n", a(2));
fprintf("a2: %d\n", a(3));
fprintf("a3: %d\n", a(4));
fprintf("a4: %d\n", a(5));
fprintf("a5: %d\n", a(6));
fprintf("a6: %d\n", a(7));

% plotta datan i figur 3
figure(3);
new_p = new_T(a, x);
plot(x, new_p), hold on
plot(x, Tdm, 'o')
xlabel('t, tiden i dygn')
ylabel('T, temperaturen i grader')

% UPPGIFT 1d
new_r = Tdm - (new_p)';
new_s = sum(new_r.^2);
fprintf("\nMinstakvadratsumman: %d\n", new_s);

% UPPGIFT 1e
figure(4);
plot(x, new_p)
