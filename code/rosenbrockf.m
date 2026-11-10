function fx = rosenbrockf(x)

fx = sum(100*(x(2:end) - x(1:end-1).^2).^2 + (1 - x(1:end-1)).^2, "all");

end