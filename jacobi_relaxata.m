% m = 10;
% A = zeros(m,m);
% b = zeros(m,m);
% bb = zeros(1,m);
% p = 10;
% er = 0.0;
% eps = 0.0;
% eps = 1e-10;
% 
% for i = 1:1:m-1    
% 
%     A(i,i) = 2;
% 
% end
% 
% for i = 1:1:m-2
%    
%     A(i,i+1) = 1;
%     A(i+1,i) = 1;
%     
% end
% 
% for i = 1:1:m
% 
%     A(i,m) = 1;
% 
% end
% 
% n = 0;
% no = 0;
% t = 0.0;
% h = 0.0;
% s = 0.0;
% z = zeros(1,m);
% Q = zeros(1,m);
% x = zeros(1,m);
% y = zeros(1,m);
% q = 0.0;
% for i = 1:1:m
%     Q(i) = 0;
%     
%     for j = 1:1:m
%         Q(i) = Q(i) + abs(A(i,j));
%     end
%     
%     if(q < Q(i))
%         q = Q(i);
%     end
%     
% end
% 
% t = 2 / q;
% h = t / p;
% 
% for k = 2:1:p
% 
%     s = k * h;
%     for i = 1:1:m
%         for j = 1:1:m
%             if( i == j)
%                 b(i,j) = 1 - s * A(i,i);
%             else
%                 b(i,j) = -s * A(i,j);
%             end
%         end
%         bb(i) = s *A(i,m);
%     end
%     
%     n = 0;
%     for i = 1:1:m
%         x(i) = 0.0;
%     end
% 
%     for i = 1:1:m
%         y(i) = bb(i);
%         for j = 1:1:m
%             y(i) = y(i)+b(i,j)*x(j);
%         end
%     end
%     er = 0.0;
%     for i = 1:1:m
%         er = er + A(i,i)*(x(i)-y(i))*(x(i)-y(i));
%     end
%     er = sqrt(er);
%     n = n+1;
%     for i = 1:1:m
%         x(i) = y(i);
%     end
%     if(k == 1)
%         no = n;
%         for i = 1:1:m
%             z(i) = x(i);
%         end
%     end
%     while(er > eps)
%         for i = 1:1:m
%             y(i) = bb(i);
%             for j = 1:1:m
%                 y(i) = y(i)+b(i,j)*x(j);
%             end
%         end
%         er = 0.0;
%         for i = 1:1:m
%             er = er + A(i,i)*(x(i)-y(i))*(x(i)-y(i));
%         end
%         er = sqrt(er);
%         n = n+1;
%         for i = 1:1:m
%             x(i) = y(i);
%         end
%         if(k == 1)
%             no = n;
%             for i = 1:1:m
%                 z(i) = x(i);
%             end
%         end 
%     end
% 
%     if(n<no)
%         no = n;
%         for i = 1:1:m
%             z(i) = x(i);
%         end
%     end
% end
% 
% x
% z





%   0.8x - 0.4y = 41
%  -0.4x + 0.8y - 0.4z = 25
% -0.4y + 0.8z = 105

n = 10;
A = zeros(n);
b = ones(n,1);

for i=1:1:n
    A(i,i) = 2;
end

for i=1:1:n-1
    A(i,i+1) = 1;
    A(i+1,i) = 1;
end

% error tolerance

tol = pow2(10,-10);

%initial guess:
x0 = zeros(n,1);


%  Jacobi method
%---------------

xnew=x0;
error=1;

pas_j = 0;
while error>tol
    xold=xnew;
    pas_j = pas_j + 1;

    for i=1:length(xnew)
        off_diag = [1:i-1 i+1:length(xnew)];
        xnew(i) = 1/A(i,i)*( b(i)-sum(A(i,off_diag)*xold(off_diag)) );
    end
    error=norm(xnew-xold)/norm(xnew);
end
x_jacobian=xnew

pas_j


%Gauss?Seidel:
%---------------

lambda=1.2;
n=length(x0);
x=x0;
error=1;
iter = 0;
maxiter = 1000;
while (error>tol && iter<maxiter)
    xold=x;
    for i=1:n
        I = [1:i-1 i+1:n];
        x(i) = (1-lambda)*x(i)+lambda/A(i,i)*( b(i)-A(i,I)*x(I) );
    end
    error = norm(x-xold)/norm(x);
    iter = iter+1;
end

x_siedal =x



% gradient
maxiter = 1000;

 [gradient_m,flag,relres,iter,resvec] = pcg(A,b,tol,maxiter);
 
 gradient_m
 %iter
 
 rotatii_m = eig(A);
 rotatii_m
