function hessfxu = rosenbrockhess(x, u)

a = 1200 * x.^2;
a(1) = a(1) + 2;
a(1:end-1) = a(1:end-1) - 400 * x(2:end);
a(2:end-1) = a(2:end-1) + 202;
a(end) = 200;

b = -400 * x(1:end-1);

hessfxu = a .* u;
hessfxu(1:end-1) = hessfxu(1:end-1) + b .* u(2:end);
hessfxu(2:end) = hessfxu(2:end) + b .* u(1:end-1);

end