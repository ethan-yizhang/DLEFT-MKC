function [Sigma,Hstar,CostNew] = DLEFT_MKC_update(Hstar,HPstar,HP,TP,Sigma,GradNew,CostOld,numclass,lambda1,option)

%------------------------------------------------------------------------------%
% Initialize
%------------------------------------------------------------------------------%
gold = (sqrt(5)+1)/2 ;
SigmaInit = Sigma ;
SigmaNew  = SigmaInit; 

%---------------------------------------------------------------
% Compute reduced Gradient and descent direction
%%--------------------------------------------------------------
switch option.firstbasevariable
    case 'first'
        [val,coord] = max(SigmaNew) ;
        %[val,coord] = max(trSTp) ;
    case 'random'
        [val,coord] = max(SigmaNew) ;
        coord=find(SigmaNew==val);
        indperm=randperm(length(coord));
        coord=coord(indperm(1));
    case 'fullrandom'
        indzero=find(SigmaNew~=0);
        if ~isempty(indzero)
            [mini,coord]=min(GradNew(indzero));
            coord=indzero(coord);
        else
            [val,coord] = max(SigmaNew) ;
        end
end

GradNew = GradNew - GradNew(coord);
desc = - GradNew.* ( (SigmaNew>0) | (GradNew<0) );
desc(coord) = - sum(desc);  % NB:  GradNew(coord) = 0

%----------------------------------------------------
% Compute optimal stepsize
%-----------------------------------------------------
stepmin  = 0;
costmin  = CostOld;
costmax  = 0;

%-----------------------------------------------------
% maximum stepsize
%-----------------------------------------------------
ind = find(desc<0);
stepmax = min(-(SigmaNew(ind))./desc(ind));
deltmax = stepmax;
if isempty(stepmax) || stepmax==0
    Sigma = SigmaNew;
    CostNew = CostOld;
    return
end

[costmax,~] = DLEFT_MKC_cost(HPstar,HP,TP,stepmax,desc,SigmaInit,numclass,lambda1,option);
%-----------------------------------------------------
%  Linesearch
%-----------------------------------------------------
Step = [stepmin stepmax];
Cost = [costmin costmax];
coord = 0;
% optimization of stepsize by golden search
while (stepmax-stepmin)>option.goldensearch_deltmax*(abs(deltmax)) && stepmax > eps
    stepmedr = stepmin+(stepmax-stepmin)/gold;
    stepmedl = stepmin+(stepmedr-stepmin)/gold;
    [costmedr,~] = DLEFT_MKC_cost(HPstar,HP,TP,stepmedr,desc,SigmaInit,numclass,lambda1,option);
    [costmedl,~] = DLEFT_MKC_cost(HPstar,HP,TP,stepmedl,desc,SigmaInit,numclass,lambda1,option);
    Step = [stepmin stepmedl stepmedr stepmax];
    Cost = [costmin costmedl costmedr costmax];
    [~,coord] = min(Cost);
    switch coord
        case 1
            stepmax = stepmedl;
            costmax = costmedl;
        case 2
            stepmax = stepmedr;
            costmax = costmedr;
        case 3
            stepmin = stepmedl;
            costmin = costmedl;
        case 4
            stepmin = stepmedr;
            costmin = costmedr;
    end
end
%---------------------------------
% Final Updates
%---------------------------------
[~,coord] = min(Cost);
CostNew = Cost(coord);
step = Step(coord);
if CostNew < CostOld
    [CostNew,Hstar] = DLEFT_MKC_cost(HPstar,HP,TP,step,desc,SigmaNew,numclass,lambda1,option);
    Sigma = SigmaNew + step * desc;
Sigma(Sigma<option.numericalprecision)=0;
Sigma=Sigma/sum(Sigma);
else
    Sigma = SigmaInit;
    CostNew = CostOld;
end
