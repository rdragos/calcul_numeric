m = 10;
A = zeros(m,m);
b = zeros(m,m);
bb = ones(1,m);
p = 100;
eps = 1e-10;
for i = 1:1:m    

    A(i,i) = 2;

end

for i = 1:1:m-1
   
    A(i,i+1) = 1;
    A(i+1,i) = 1;
    
end

n = 0;
no = 0;
t = 0.0;
h = 0.0;
s = 0.0;
z = zeros(1,8);
Q = zeros(1,8);
q = 0;

for i = 1:1:m
    Q(i) = 0;
    
    for j = 1:1:m
        Q(i) = Q(i) + abs(A(i,j));
    end
    
    if(q<Q(i))
        q = Q(i);
    end
    
end

t = 2/q;
h = t/p;

for k = 1:1:p

    s = k * h;
    for i = 1:1:m
        for j = 1:1:m
            if( i == j)
                b(i,j) = 1 - s*A(i,j);
            else
                b(i,j) = -s*A(i,j);
            end
        end
        bb(i) = s *A(i,m);
    end
    
    n = 0;
    for i = 1:1:m
        x(i) = 0.0;
    end

    do
        for i = 1:1:m
            y(i) = bb(i);
            for j = 1:1:m
                y(i) = y(i)+b(i,j)*x(j);
            end
        end
        er = 0;
        for i = 1:1:m
            er = er + A(i,i)*(x(i)-y(i))*(x(i)-y(i));
        end
        er = sqrt(er);
        n = n+1;
        for i = 1:1:m
            x(i) = y(i);
        end
        if(k == 1)
            no = n;
            for i = 1:1:m
                z(i) = x(i);
            end
        end
    until(er > eps)

    if(n<no)
        no = n;
        for i = 1:1:m
            z(i) = x(i);
        end
    end
end

x
