function [Hstar,error_HPstar]=DLEFT_MKC(KH,numclass,lambda1,lambda2,option)
numker = size(KH,3);
num = size(KH,2); %sample number
sX = [num,numclass,numker];
qnorm = 2;

%% initial the coefficients
gamma = ones(numker+1,1)/(numker+1);
avgSigma = ones(numker,1)/numker;

%% initial H0,HP,TP, H*,HP*
Kmatrix = sumKbeta(KH,avgSigma.^2);
H0 = mykernelkmeans(Kmatrix,numclass);
Hstar = zeros(num,numclass);
HP = zeros(num,numclass,numker+1);
HPstar = zeros(num,numclass,numker+1);
TP = zeros(numclass,numclass,numker+1);
A1 = zeros(sX);Y1 = zeros(sX);

for p=1:numker
   HP(:,:,p) = mykernelkmeans(KH(:,:,p),numclass);
   TP(:,:,p) = eye(numclass); 
   A1(:,:,p) = HP(:,:,p);
   Y1(:,:,p) = zeros(num,numclass);
end
p = numker + 1;
HP(:,:,p) = H0;
TP(:,:,p) = eye(numclass); 

Hmatrix = zeros(num,numclass);
for p =1:numker+1
    Hmatrix = Hmatrix + gamma(p).^2*HP(:,:,p)*TP(:,:,p);
end
[U,~,V] = svd(Hmatrix,'econ');
Hstar = U*V';


    
%% algorithm

    
Isconverg = 0;
epson = 1e-4;
iter = 1;
mu = 10e-1; 
max_mu = 10e10; pho_mu = 2;

while(Isconverg == 0)
%     fprintf('----processing iter %d--------\n', iter+1);
    HPstarOld=HPstar;
    gammaOld=gamma;
        %% step1 (initial) update Hp*      
    for p =1:numker   
        HPmatrix = gamma(p).^2*Hstar*TP(:,:,p)'+lambda1*gamma(p).^2*HP(:,:,p)+mu*(A1(:,:,p)-1/mu*Y1(:,:,p));
        [U,~,V] = svd(HPmatrix,'econ');
        HPstar(:,:,p) = U*V'; 
    end
    p = numker + 1;
    HPstar(:,:,p) = HP(:,:,p);
    
%% step2 update gamma and H*     by Reduced Gradient Descent Algorithm 

    Hmatrix = zeros(num,numclass);
H2matrix = zeros(numclass,numclass);
for p =1:numker+1
    Hmatrix = Hmatrix + gamma(p).^2*HPstar(:,:,p)*TP(:,:,p);
    H2matrix = H2matrix + gamma(p).^2*HPstar(:,:,p)'*HP(:,:,p);
end
[U,~,V] = svd(Hmatrix,'econ');
Hstar = U*V';
     costOld = trace(Hstar'*Hmatrix+lambda1*H2matrix);
     loop = 1;
nloop = 1;
obj = [];
obj(nloop) = costOld;
while loop
    nloop = nloop + 1;
    gamma_rgd_old = gamma;
     [grad] = DLEFT_MKC_grad(gamma,Hstar,HPstar,TP,HP,lambda1);    
     [gamma,Hstar,obj(nloop)] = DLEFT_MKC_update(Hstar,HPstar,HP,TP,gamma,grad,obj(nloop-1),numclass,lambda1,option);
  
      if max(abs(gamma-gamma_rgd_old))<option.seuildiffsigma 
        loop = 0;
    end
    if (nloop>3 &&(obj(nloop-1)-obj(nloop))/obj(nloop)<1e-4 && (obj(nloop-2)-obj(nloop-1))/obj(nloop-1)<1e-4 )
        loop = 0;
    end
end
   
%% step3 update TP       
for p =1:numker+1   
        [U,~,V] = svd(HPstar(:,:,p)'*Hstar,'econ');
        TP(:,:,p) = U*V'; 
end
%% step4  update A

    a = wshrinkObj(HPstar(:,:,1:numker)+Y1*1/mu,lambda2/mu,sX,0,3);
    A1 = reshape(a,sX);
    
    
     %% Update Lagrange multiplier
    Y1 = Y1 + mu*(HPstar(:,:,1:numker)-A1);
    mu = mu*pho_mu;
    
    %% coverge condition    
    HPs1 = HPstar(:,:,1:numker);
    diffA1 = max(abs(A1(:)-HPs1(:)));
    diffHPstar = max(abs(HPstarOld(:)-HPstar(:)));
    diffGamma = max(abs(gammaOld(:)-gamma(:)));
    maxDiff = max([diffA1,diffHPstar,diffGamma]);
    error_HPstar(iter) = diffHPstar;
   
     if maxDiff < epson
        Isconverg  = 1;
    end   
    if (iter>30)
        Isconverg  = 1;
    end
    
   
   iter = iter + 1;
end


end

