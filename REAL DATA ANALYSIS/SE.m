function result=SE(type,W,Y,hatSigma)
[p,n]=size(Y);


for i=1:n
    hatSigmaM1o2=hatSigma(:,:,i)^(-1/2);
%    if issymmetric(hatSigmaM1o2)==0
%        hatSigmaM1o2=(hatSigmaM1o2+hatSigmaM1o2')/2;
%    end
    tildeY(:,i)=hatSigmaM1o2*Y(:,i);  
end

mu4=mean(mean(tildeY.^4));

if strcmpi(type , 'OLS')
    P1=Pd(1,W,Y,hatSigma);
    Q1=Qd(1,W,Y,hatSigma);
    Q0=Qd(0,W,Y,hatSigma);
    Q0inv=inv(Q0);
    result= sqrt(diag(2*Q0inv*Q1*Q0inv+(mu4-3)*Q0inv*P1*Q0inv)/n/p);
    
elseif strcmpi(type , 'GLS')
    Pm1=Pd(-1,W,Y,hatSigma);
    Qm1=Qd(-1,W,Y,hatSigma);
    Qm1inv=inv(Qm1);
%    if issymmetric(Qm1inv)==0
%        Qm1inv=(Qm1inv+Qm1inv')/2;
%    end
    result= sqrt(diag(2*Qm1inv+(mu4-3)*Qm1inv*Pm1*Qm1inv)/n/p);
    
elseif strcmpi(type , 'MLE') 
    Pm1=Pd(-1,W,Y,hatSigma);
    Qm1=Qd(-1,W,Y,hatSigma);
    Qm1inv=inv(Qm1);
%    if issymmetric(Qm1inv)==0
%        Qm1inv=(Qm1inv+Qm1inv')/2;
%    end
    result= sqrt(diag(2*Qm1inv+(mu4-3)*Qm1inv*Pm1*Qm1inv)/n/p);    

else
    error('The estimator does not exit!');
end

