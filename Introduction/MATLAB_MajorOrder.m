n = 10000;
matrice = rand(n, n);

%row
tic;
for i = 1:n
    for j = 1:n
        valeur = matrice(i, j);
    end
end
fprintf('Ligne : %.3f secondes\n', toc);

% col
tic;
for i = 1:n
    for j = 1:n
        valeur = matrice(j, i);
    end
end
fprintf('Colonne : %.3f secondes\n', toc);
