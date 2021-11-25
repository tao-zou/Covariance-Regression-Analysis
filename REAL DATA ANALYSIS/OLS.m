function [AA,hatbeta]=OLS(W,Y)
[p,n]=size(Y);
sizeW=size(W);
K=sizeW(4);


for i=1:n
    aa=zeros(K,K);
    for k1=1:K
        for k2=1:k1
            aa(k1,k2)=trace(W(:,:,i,k1)*W(:,:,i,k2));
        end
    end
    aa=aa+tril(aa,-1).';
    aaa(:,:,i)=aa;
end

AA=sum(aaa,3);





BB=zeros(K,1);
for k=1:K
    for i=1:n
        bb(i)=Y(:,i)'*W(:,:,i,k)*Y(:,i);
    end
    BB(k)=sum(bb);
end

hatbeta=AA\BB;
