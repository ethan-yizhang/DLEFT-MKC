clear
% clc
warning off;


    dataName = '*';
	load(['F:\work2015\datasets\',dataName,'_Kmatrix'],'KH','Y');
    numclass = length(unique(Y));
    Y(Y<1)=numclass;
    numker = size(KH,3);
    num = size(KH,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    KH = kcenter(KH);
    KH = knorm(KH);
    options.seuildiffsigma=1/numker * 1e-4;        % stopping criterion for weight variation
    %------------------------------------------------------
    % Setting some numerical parameters
    %------------------------------------------------------
    options.goldensearch_deltmax=1e-3; % initial precision of golden section search
    options.numericalprecision=1e-16;   % numerical precision weights below this value
    % are set to zero
    %------------------------------------------------------
    % some algorithms paramaters
    %------------------------------------------------------
    options.firstbasevariable='first'; % tie breaking method for choosing the base
    % variable in the reduced gradient method
    options.nbitermax=500;             % maximal number of iteration
    options.seuil=0;                   % forcing to zero weights lower than this
    options.seuilitermax=10;           % value, for iterations lower than this one
    
    options.miniter=0;                 % minimal number of iterations
    options.threshold = 1e-4;
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    qnorm = 2;
    options.display = 0;
    
    
    
	%% Optimal Parameters:
%	            lambda:		rho:
%	bbcsport:	100 		32
%	CCV:		100			32
%	flower102:	100			32
%	liver:		1000		128
%	plant:		1			8
%	proteinFold:0.1			4
%	psortNeg:	0.001		8
%	Reuters:	1000		128
%	scene15:	0.01		16
%	willow:		1			4
	
    %% DLEFT-MKC
    lambda1set = 10.^[-3:5];
    lambda2set = 2.^[-1:7];
    res_allmean = zeros(4,length(lambda2set),length(lambda1set));
    res_allstd = zeros(4,length(lambda2set),length(lambda1set));
    tic;
    for il1 = 1:length(lambda1set)
        for il2 = 1:length(lambda2set)
            [Hstar,~] = DLEFT_MKC(KH,numclass,lambda1set(il1),lambda2set(il2),options);
            [res_mean1,res_std1] = myNMIACCV2(Hstar,Y,numclass);
            res_allmean(:,il2,il1)=res_mean1;
            res_allstd(:,il2,il1)=res_std1;
        end
    end
    [~,max_idx]=max(res_allmean(1,:,:),[],'all','linear');
    res_mean(:,1) = res_allmean(:,max_idx);
    res_std(:,1) = res_allstd(:,max_idx);
    timecost(1) = toc/length(lambda2set)/length(lambda1set);
    
    precision = reshape(res_allmean(1, :, :), size(res_allmean,2), size(res_allmean,3)); % 将第一行转换为一个行向量
    [~,max_idx2]=max(precision,[],'all','linear');% 找出最大聚类精度对应的位置
    [max_il2 ,max_il1]=ind2sub(size(precision),max_idx2);% 推导出超参数 p 和 q 的大小
    fprintf('ACC: %.4g, Paramter:  %.5g , %.5g &',res_mean(1,1),lambda1set(max_il1),lambda2set(max_il2));
    
   