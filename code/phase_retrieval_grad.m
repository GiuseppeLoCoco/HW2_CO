function gradfx = phase_retrieval_grad(x, A, y)
m = size(A, 1);
Ax = A * x;
residual = (Ax).^2 - y.^2;
gradfx = (1 / m) * A' * (residual .* Ax);
end