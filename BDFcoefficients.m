function BDFcoefficients(k)

% The BDF coefficients if we have b_0=1 and calculating with operator nabla

A=ones(k+1,k+1);
A=tril(A);

A=A([2:end],:);
for i=2:k+1-1
    for j=2:i
        A(i,j)=(nchoosek(i-1,j-2)+nchoosek(i-1,j-1));
    end
end
D=spdiags(1./[1:1:k+1]',0,k+1,k+1-1);
BDFvector=sum(D*A,1);
BDFvector(2:2:end)=-BDFvector(2:2:end);
rats(BDFvector)
