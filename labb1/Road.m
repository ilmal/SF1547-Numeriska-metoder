format long

% -----------------
% Uppgift 3a

% Startgissningar i matrisformat baserat på plottad graf
x0 = [204 1001; 458 2457; 712 1748]; 

% Vektorer med värden för (xa, ya) samt la (längden)
xa = [175, 410, 675]; ya = [950, 2400, 1730]; 
la = [60, 75, 42];

% Vektorer med värden för (xb, yb) samt lb (längden)
xb = [160, 381, 656]; yb = [1008, 2500, 1760]; 
lb = [45, 88, 57];
 
% Plotta cirklarna (nivåkurva 0) för att hitta högra skärningspunkt (gissningar)
figure;
hold on;
colors = ['r', 'g', 'b']; % Färger för olika cirkel-'par'

% Loop genom alla par av cirklar 
for i = 1:length(xa)
    % Ekvationerna för cirklarna
    p_a = @(xp, yp) (xp - xa(i)).^2 + (yp - ya(i)).^2 - la(i).^2;
    p_b = @(xp, yp) (xp - xb(i)).^2 + (yp - yb(i)).^2 - lb(i).^2;

    % Plotta cirklarna med olika färger beroende på set
    fimplicit(p_a, [xa(i)-la(i)-10, xa(i)+la(i)+10, ya(i)-la(i)-10, ya(i)+la(i)+10], colors(i), 'LineWidth', 1.5);
    fimplicit(p_b, [xb(i)-lb(i)-10, xb(i)+lb(i)+10, yb(i)-lb(i)-10, yb(i)+lb(i)+10], colors(i), 'LineWidth', 1.5);

    % Markera centrum i cirklarna
    plot(xa(i), ya(i), 'o', 'MarkerSize', 8, 'MarkerFaceColor', colors(i));
    plot(xb(i), yb(i), 'o', 'MarkerSize', 8, 'MarkerFaceColor', colors(i));
    
    % Text för vilken punkt det är (P)
    text(xa(i) + 10, ya(i), ['Point ' num2str(i)], 'Color', colors(i), 'FontSize', 12);
    text(xb(i) + 10, yb(i), ['Point ' num2str(i)], 'Color', colors(i), 'FontSize', 12);
end

% punkter P 1-3 för for-loop (slippa kod-duplikation)
for i = 1:length(xa)
    % Printa vilken punkt som ska hittas i denna iteration
    fprintf('\n=== Solving for Point P%d ===\n', i);

    % Funktionerna för cirklarna i vektorform F
    F = @(xp, yp) [
        (xp - xa(i)).^2 + (yp - ya(i)).^2 - la(i).^2;
        (xp - xb(i)).^2 + (yp - yb(i)).^2 - lb(i).^2
    ];

    % Jakobianen till F (alltså derivatorna till F's ekvationer med avseende på x eller y)
    J = @(xp, yp) [
        2 * (xp - xa(i)), 2 * (yp - ya(i)); 
        2 * (xp - xb(i)), 2 * (yp - yb(i))
    ];
    
    % Startgissning är radvektorn i x0 transponerad till kolonnvektor
    x_guess = x0(i, :)'; % ':' tar alla kolumner på rad nr 'i'
    
    % Rimlig tolerans och iterationsräknare
    tol = 10^(-6); iter = 0; max_iter = 1000;
    diff = 1; % Startdifferens som är större än toleransen

    % Print header för resultaten
    fprintf('Iter |      x       |      y       |   diff   \n');
    fprintf('---------------------------------------------\n');
    
    % printa fel
    errors = [];

    % Newtons metod applicerad på system
    while iter < max_iter && diff > tol
        iter = iter+1; % inkrementera räknaren för iterationer
        h = - J(x_guess(1), x_guess(2)) \ F(x_guess(1), x_guess(2)); % Newtons metod för system
        new_x = x_guess + h; % Hitta nya x
        diff = norm(h); % differensen är normen av h
        x_guess = new_x; % sätt nya x som gissning
        errors = [errors, diff]; % spara felet

        % Print resultat för varje iteration
        fprintf('%4d | %10.6f | %10.6f | %8.2e \n', iter, new_x(1), new_x(2), diff);
    end
    
    % plotta felet norm(h) -> kvadratisk konvergenshastighet
    figure; 
    semilogy(1:length(errors), errors, 'b-o', 'DisplayName', 'Newton'); % newton
    xlabel('Iteration');
    ylabel('norm(h)');
end

% Konvergenshastigheten för felen är kvadratisk eftersom kvoten är ca 64
% Enligt teori: |x(n+1) - x(n)| = |h(n+1)| ≈ K•|h(n)|^2
% |h(n+1)|/|h(n)|^2 ≈ 64
% då kvoten stämmer väl överens med teorin så konvergerar newton som det
% ska

% -----------------
% Uppgift 3b

% Anonym funktion av fjärde grad
p4 = @(c, x) c(1) + c(2)*x + c(3)*x.^2 + c(4)*x.^3 + c(5)*x.^4;

% Koordinaterna för punkterna P1, P2, P3
x = [204.6 458.1 712.1];
y = [1002.2 2457.6 1749.7]'; 

% Matrisen A
A = [ones(size(x))' x' x'.^2 x'.^3 x'.^4];

% Lös c
c = A\y;

% Printa koefficienternas värde
fprintf('\n');
disp("Koefficientvärden c: ");
fprintf("c0: %d\n", c(1));
fprintf("c1: %d\n", c(2));
fprintf("c2: %d\n", c(3));
fprintf("c3: %d\n", c(4));
fprintf("c4: %d\n", c(5));

% Plotta interpolationen och punkterna i en graf
x_plot = linspace(1, 1020, 1020);
func_plot = p4(c, x_plot);
figure;
plot(x_plot, func_plot), hold on
plot(x, y, 'o');
