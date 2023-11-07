function [Acomplete,difference] = AltMin(AA1,MASK,r,maxIter,tol,show)


% Iteration log will show if show = true.
m=225;
n=225;
V=rand(r,n);
U=rand(m,r);
A=U*V;

[m,n]=size(AA1);
NumOfSample = sum(MASK(:));
ROW = Row_Nonzero(A); % Nonzeros entry index for each row
COLUMN = Column_Nonzero(A); % Nonzeros entry index for each column

[U,~,~] = svd(A);
U = U(:,1:r);
difference = [];

for t = 1:maxIter
    if show
        display(t);
    end 
    V = zeros(r,n);
    for j = 1:n
        pos = cut_zero(COLUMN(:,j));
        V(:,j) = U(pos,:) \ AA1(pos,j);
        if r == 1 && isinf(V(j))
            V(j) = 0;
        end
    end
   temp = norm((U*V-AA1),'fro')/norm(AA1,'fro');
    difference = [difference temp];
    if temp < tol
        break;
    end
    
    U = zeros(m,r);
    for i = 1:m
        pos = cut_zero(ROW(i,:));
        u = V(:,pos)' \ AA1(i,pos)';
        U(i,:) = u';
        if r == 1 && isinf(U(i))
            U(i) = 0;
        end
    end
   temp = norm((U*V-AA1),'fro')/norm(AA1,'fro');
    difference = [difference temp];
    if temp < tol
        break;
    end
    
end

Acomplete = U*V;