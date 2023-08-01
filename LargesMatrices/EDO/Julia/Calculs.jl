function Sol(x,y)::Float64
    t2 = x ^ 2
    t6 = cos(pi * y)
    return sin(pi * x + t6 * (t2 - x) * beta)
end

function champMagNormaliseRenverse(u::SVector{2,Float64}, p, t)
    x = u[1];
    y = u[2];
    
    t1 = beta * pi;
    t3 = -1 + x;
    t4 = pi * y;
    t5 = sin(t4);
    t7 = beta ^ 2;
    t8 = x ^ 2;
    t10 = pi * (t8 - x);
    t11 = 2 * x;
    t15 = cos(t4);
    t16 = t15 ^ 2;
    t23 = pi ^ 2;
    t24 = t3 ^ 2;
    t30 = sqrt(-t16 * (t10 - t11 + 0.1e1) * (t10 + t11 - 0.1e1) * t7 + 0.4e1 * t15 * (x - 0.1e1 / 0.2e1) * t1 + (t7 * t24 * t8 + 0.1e1) * t23);
    MatBx = 0.1e1 / t30 * t5 * t3 * x * t1;
    
    t1 = 2 * x;
    t5 = cos(pi * y);
    t8 = beta ^ 2;
    t9 = x ^ 2;
    t11 = pi * (t9 - x);
    t15 = t5 ^ 2;
    t23 = pi ^ 2;
    t25 = (-1 + x) ^ 2;
    t31 = sqrt(-t15 * (t11 - t1 + 0.1e1) * (t11 + t1 - 0.1e1) * t8 + 0.4e1 * t5 * (x - 0.1e1 / 0.2e1) * pi * beta + (t8 * t25 * t9 + 0.1e1) * t23);
    MatBy = 0.1e1 / t31 * (t5 * (t1 - 1) * beta + pi);
    
    return @SVector [-MatBx; -MatBy]
end