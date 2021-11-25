function matrix=Ximat(beta,W)
for i=1:length(beta)
    TW(:,:,i)=beta(i)*W(:,:,i);
end
matrix=sum(TW,3);