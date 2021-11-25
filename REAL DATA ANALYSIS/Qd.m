function AA=Qd(d,W,Y,Sigma)
[p,n]=size(Y);
sizeW=size(W);
K=sizeW(4);


for i=1:n
    Sigmad=Sigma(:,:,i)^d;
    aa=zeros(K,K);
    for k1=1:K
        for k2=1:k1
            aa(k1,k2)=trace(Sigmad*W(:,:,i,k1)*Sigmad*W(:,:,i,k2));
        end
    end
    aa=aa+tril(aa,-1).';
    aaa(:,:,i)=aa;
end

AA=mean(aaa,3)/p;



