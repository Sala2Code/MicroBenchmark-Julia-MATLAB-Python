function SolBetareturn = SolBeta(x, y, beta, epsilon, nPeriod)
  t1 = pi * x;
  t2 = x ^ 2;
  t5 = pi * y;
  t7 = cos(0.3e1 * t5);
  t11 = sin((t1 + t7 * (t2 - x) * beta) * nPeriod);
  t13 = cos(0.2e1 * t5);
  t15 = sin(t1);
  SolBetareturn = epsilon * t13 * t15 + t11;
