function KH = imputeKH_ISMKKM_DRGM_normalized(KH,WP,mis_set,obs_set)

numker = size(KH,3);
num = size(KH,1);
for p = 1 : numker
%     mis_setp = mis_set{p};
%     obs_setp = obs_set{p};
%     Wp = zeros(length(obs_set),length(mis_set));
%     Wp = WP{p};
    % WP{p} =  [eye(length(mis_set));zeros(length(obs_set)-length(mis_set),length(mis_set))];
    KH(obs_set{p},obs_set{p},p) = KH(obs_set{p},obs_set{p},p);
    KH(obs_set{p},mis_set{p},p) = KH(obs_set{p},obs_set{p},p)*WP{p};
    KH(mis_set{p},obs_set{p},p) = KH(obs_set{p},mis_set{p},p)';
    KH(mis_set{p},mis_set{p},p) = KH(mis_set{p},obs_set{p},p)*WP{p};
    
    %% diag to be 1
    tempx = sqrt(diag(KH(mis_set{p},mis_set{p},p) ));% m * 1
    KH(:,mis_set{p},p) = KH(:,mis_set{p},p) ./ repmat(tempx',num,1);
    KH(mis_set{p},:,p) = KH(mis_set{p},:,p) ./ repmat(tempx,1,num);
%     KH(mis_set{p},mis_set{p},p) =  KH(mis_set{p},mis_set{p},p) ./repmat(tempx,1,length(mis_set{p}));
    
    KH(:,:,p) = (KH(:,:,p) + KH(:,:,p)')/2;
%     KH(mis_set{p},mis_set{p},p) = eye(length(mis_set{p}));
end