clear all
a = rand();
b = rand();

%% sum
t=0;
n=1000000;
for i = 1:n
    tic();
    sum(a, b);
    % prod(a,b);
    % logic(a,b);
    t = t + toc();
end
fprintf("Temps moyen (s) : %.12f", t/n);





%% functions
function res = sum(a,b)
    res = a+b;
end

function res = prod(a,b)
    res = a * b;
end

function res = logic(a,b)
    if a>0.5 && b<0.5
        res = a;
    elseif a<0.25 || b>0.5
        res = b;
    else
        res = 0;
    end
end

%%