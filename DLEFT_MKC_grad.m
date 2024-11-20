function [grad] = DLEFT_MKC_grad(gamma,Hstar,HPstar,TP,HP,lambda1)

d=length(gamma);
grad=zeros(d,1);
for p=1:d
     grad(p) = 2*gamma(p)*trace(Hstar'*HPstar(:,:,p)*TP(:,:,p)+lambda1*HPstar(:,:,p)'*HP(:,:,p));  
end
grad = grad / (d-1);