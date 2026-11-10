function gradfx = rosenbrockgrad(x)

y = x(2:end) - x(1:end-1).^2;

gradfx = -2 * (1 - x);

gradfx(1:end-1) = gradfx(1:end-1) - 400 * x(1:end-1) .* y;
gradfx(2:end) = gradfx(2:end) + 200 * y;

gradfx(end) = 200 * (x(end) - x(end-1)^2);

end