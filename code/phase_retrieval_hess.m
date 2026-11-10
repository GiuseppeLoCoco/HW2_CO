function hessfxu = phase_retrieval_hess(x, u, A, y)
m = size(A, 1);
Ax = A * x;
Av = A * u;
coeff = 3 * (Ax).^2 - y.^2;
hessfxu = (1 / m) * A' * (coeff .* Av);
end