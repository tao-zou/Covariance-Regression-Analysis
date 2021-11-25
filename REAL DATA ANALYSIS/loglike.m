function v=loglike(beta,W,Y)
[p,n]=size(Y);

for i=1:n
    WI(:,:,:)=W(:,:,i,:);
    
    Sigma=Ximat(beta,WI);
    [Gamma,Lambda]=eig(Sigma);
    diagLambda=diag(Lambda);
    Omega=Gamma*diag(1./diagLambda)*Gamma';
    
    quar(i)=sum(log(diagLambda))+Y(:,i)'*Omega*Y(:,i);
    
end


v=sum(quar);