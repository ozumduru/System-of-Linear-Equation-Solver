function[x]=Jacobi(xi,A,b,accuracy,MaxIter)

error = 1;
j=1;

[B_JM,C] = B_matrix_and_C(A,b);
if abs(det(B_JM)) < 1
    
while error > accuracy

xi(:,j+1)=B_JM*xi(:,j) + C;

difference = xi(:,j+1) - xi(:,j);
error = sqrt(difference.^2);

j=j+1;
if j > MaxIter
    break
end

end
x=xi(:,j);
else
'B_JM not converge'
end
end