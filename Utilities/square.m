function y = square(x)

inp = sin(x) >= 0;
y(~inp) = -1;
y(inp) = 1;

end