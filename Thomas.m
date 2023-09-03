function h = Thomas(a,b,c,r)

% Implements Thomas algorithm

N = length(b);
w = zeros(N, 1);
g = zeros(N, 1);
w(1) = c(1)/b(1);
g(1) = r(1)/b(1);
if isrow(c)
    c = c';
end
if isrow(a)
    a = a';
end
c = [c; 0]; a = [0; a];
for i=2:N
    w(i) = c(i)/(b(i)-a(i)*w(i-1));
    g(i) = (r(i)-a(i)*g(i-1))/(b(i)-a(i)*w(i-1));
end
h = zeros(N, 1);
h(N) = g(N);
for i=N-1:-1:1
    h(i) = -w(i)*h(i+1)+g(i);
end
end