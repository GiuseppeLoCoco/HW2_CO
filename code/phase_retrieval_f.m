function fx = phase_retrieval_f(x, A, y)
m = size(A, 1);
Ax = A * x;
fx = (1 / (4 * m)) * sum(((Ax).^2 - y.^2).^2);
end