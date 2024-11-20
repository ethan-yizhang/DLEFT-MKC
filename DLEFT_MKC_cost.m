function [cost,Hstar] = DLEFT_MKC_cost(HPstar,HP,TP,StepSigma,DirSigma,gamma,numclass,lambda1,option)

global nbcall
nbcall=nbcall+1;

gamma = gamma+ StepSigma * DirSigma;
gamma(gamma<option.numericalprecision)=0;
gamma=gamma/sum(gamma);

num = size(HPstar,1);
numker = size(HPstar,3)-1;
Hmatrix = zeros(num,numclass);
H2matrix = zeros(numclass,numclass);
for p =1:numker+1
    Hmatrix = Hmatrix + gamma(p).^2*HPstar(:,:,p)*TP(:,:,p);
    H2matrix = H2matrix + gamma(p).^2*HPstar(:,:,p)'*HP(:,:,p);
end
[U,~,V] = svd(Hmatrix,'econ');
Hstar = U*V';
cost = trace(Hstar'*Hmatrix+lambda1*H2matrix);
