function [U,V,obj]= mykernelsvd(P,cluster_count)


opt.disp = 0;
[U,~,V] = svds(P,cluster_count);
obj = trace(U' * P * V);
% H_normalized = H;