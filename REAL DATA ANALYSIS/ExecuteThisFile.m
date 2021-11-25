%% 20211125
clear all;

%% Read DATA
Y=csvread('Ystock.csv',1,0);
Y=Y-mean(Y);

W_IND=csvread('W_IND.csv',1,0);
W_LOCATION=csvread('W_LOCATION.csv',1,0);
W_SIZE=csvread('W_SIZE.csv',1,0);
W_PRICE=csvread('W_PRICE.csv',1,0);




%% Parameters
p=length(Y);
n=1;


%OLS+ & GLS+
epsilon=1e-5;
mu=0.05;
xi=1e-8;
Lambda0=zeros(p,p,n);
iterMAX=1000;



%MLE
optim_options = optimset('TolX',1e-2,'TolFun',1e-3,'FunValCheck','on');



%% PreCompute
K=4;

%%


W1=speye(p);

for i=1:n
    W(:,:,i,1)=full(W1);
    W(:,:,i,2)=full(W_IND);
    W(:,:,i,3)=full(W_SIZE);
    W(:,:,i,4)=full(W_PRICE);
end




%for loop=1:1
    t=cputime;
 
    
    t1=cputime;
    
    %OLS
    [AA,hatbetaOLS]=OLS(W,Y);
    
    %OLS+
    [hatbetaOLSp,hatSigmaOLSp,iterOLSp]=OLSp(W,Y,AA,hatbetaOLS,Lambda0,mu,epsilon,xi,iterMAX);
    
    
    SEOLSp=SE('OLS',W,Y,hatSigmaOLSp);
    
    %GLS
    for i=1:n
        hatSigmaM1o2=hatSigmaOLSp(:,:,i)^(-1/2);
        if issymmetric(hatSigmaM1o2)==0
            hatSigmaM1o2=(hatSigmaM1o2+hatSigmaM1o2')/2;
        end
        tildeY(:,i)=hatSigmaM1o2*Y(:,i);
        for k=1:K
            tildeW(:,:,i,k)=hatSigmaM1o2*W(:,:,i,k)*hatSigmaM1o2;
            if issymmetric(tildeW(:,:,i,k))==0
                tildeW(:,:,i,k)=(tildeW(:,:,i,k)+tildeW(:,:,i,k)')/2;
            end
        end       
    end
    


    [AA,hatbetaGLS]=OLS(tildeW,tildeY);

    %GLS+
    [hatbetaGLSp,hatSigmaOLSpt,iterGLSp]=OLSp(tildeW,tildeY,AA,hatbetaGLS,Lambda0,mu,epsilon,xi,iterMAX);

    for i=1:n
        WI(:,:,:)=W(:,:,i,:);
        hatSigmaGLSp(:,:,i)=Ximat(hatbetaGLSp,WI);
    end    
    
    SEGLSp=SE('GLS',W,Y,hatSigmaGLSp);
    
    
    [iterOLSp,iterGLSp]
    
    [hatbetaOLS,hatbetaOLSp,hatbetaGLS,hatbetaGLSp]
    
    [SEOLSp,SEGLSp]
    
    [iterOLSp,iterGLSp]
    [hatbetaGLSp*100,SEGLSp*100,hatbetaGLSp./SEGLSp,(1-cdf('Normal',abs(hatbetaGLSp./SEGLSp),0,1))]
    
    %MLE
    %Without Derivatives
    %objfun = @(param)(loglike(param,W,Y));
    %[hatbetaMLE,fval,exitflag] = fminsearch(objfun,[1;0;0;0;0;0],optim_options);
    
    %for i=1:n
    %    WI(:,:,:)=W(:,:,i,:);
    %    hatSigmaMLE(:,:,i)=Ximat(hatbetaMLE,WI);
    %end    
    
    %SEMLE=SE('MLE',W,Y,hatSigmaMLE);
    
    %t2=cputime-t1;
    
    %[iterOLSp,iterGLSp,exitflag,t2]
    
    %[hatbetaOLS,hatbetaOLSp,hatbetaGLS,hatbetaGLSp,hatbetaMLE]
    
    %[SEOLSp,SEGLSp,SEMLE]

    
   


%end
