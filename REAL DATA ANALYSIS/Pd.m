function AA=Pd(d,W,Y,Sigma)
[p,n]=size(Y);
sizeW=size(W);
K=sizeW(4);


for i=1:n
    Sigmad=Sigma(:,:,i)^(d/2);
    aa=zeros(K,K);
    for k1=1:K
        for k2=1:k1
            aa(k1,k2)=sum(diag(Sigmad*W(:,:,i,k1)*Sigmad).*diag(Sigmad*W(:,:,i,k2)*Sigmad));
        end
    end
    aa=aa+tril(aa,-1).';
    aaa(:,:,i)=aa;
end

AA=mean(aaa,3)/p;



