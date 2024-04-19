function [X,iter,Res]=updateX(Y_p,A_p,X_true,lambda,rho1,miter)
%%
Z=zeros(size(X_true));M=zeros(size(X_true));
res=1;
Res=[];
Rel=[];
tol=10^(-3);
iter=0;
mu1=20;
X_p=zeros(size(X_true));
X_p(1,1)=1;
while res>tol && iter<miter

    temp=A_p*A_p'+(rho1+mu1)*(eye(size(A_p*A_p')));
    X=(Y_p*A_p'+rho1*X_true+mu1*Z-M)*pinv(temp);
Z= solve_l1l2(X+M/mu1, lambda/mu1);
    res=norm(X-X_p,'fro')/norm(X_p,'fro');
    Res=[Res,res];
    X_p=X;
M=M+mu1*(X-Z);
    %%% updating mu1  
% mu1=min(mu1*1.5,1e6);

    iter=iter+1;
end




