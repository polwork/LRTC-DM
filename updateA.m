function [A,iter,Res]=updateA(Y_p,A_true,X_true,lambda,rho1,miter)
%%
Z=zeros(size(A_true));M=zeros(size(A_true));
res=1;
Res=[];
Rel=[];
tol=10^(-3);
iter=0;
mu1=20;
A_p=zeros(size(A_true));
A_p(1,1)=1;
while res>tol && iter<miter
    temp=X_true'*X_true+(rho1+mu1)*eye(size(X_true'*X_true));
    A= pinv(temp)*(X_true'*Y_p+rho1*A_true+mu1*Z-M);
        Z= solve_l2l1(A+M/mu1, lambda/mu1);
    res=norm(A-A_p,'fro')/norm(A_p,'fro');
    Res=[Res,res];
    A_p=A;
M=M+mu1*(A-Z);
    %%% updating mu1  
% mu1=min(mu1*1.5,1e6);

    iter=iter+1;
end




