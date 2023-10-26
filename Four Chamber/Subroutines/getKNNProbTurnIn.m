function [turnBias,XX,YY,turnBias_baseline,turnBias_during_baseline] =  ...
    getKNNProbTurnIn(f_orco, s, keyNdx, xGrid, yGrid, zGrid, ratio, Thresh, K, state, genotype, plotFig)%f_orco, condNdx
% condNdx3 = zeros(size(f_orco.spk));
% condNdx3(condNdx<6) = 1;
% condNdx3(condNdx==0) = 3;% not an condition case
% contIdx = 2;% index of the control group (before light on period)
% fe = f_orco.getFirstEntry('H',1.3);
% for i = 1:f_orco.nFly
%     condNdx3(i,1:fe(i)-1) = contIdx;% before case
% end
% key = cellstr(num2str(unique(condNdx3)));thresh = 0.1;
% [s, keyNdx] = f_orco.getKinematicsCond(condNdx3,key,thresh);

% calculate time from first entry for the beginning of each state
fe = f_orco.getFirstEntry('H',1.3);
for fly = 1:size(s.ndx,1)
    for j = 1:size(s.ndx,2)
        s.timeFromFe{fly,j} = s.ndx{fly,j}-fe(fly);
    end
end

if plotFig
    sbpltY = ceil(sqrt(numel(state)));
    sbpltX = round(sqrt(numel(state)));
    figure;set(gcf,'Position',[6 39 736.*sbpltX./sbpltY 736])
end
for currState = 1:numel(state)

    l = cellfun(@(c)strcmp(c,keyNdx),{state(currState)},'UniformOutput',false);
    l2 = cellfun(@(c)strcmp(c,keyNdx),{'1'},'UniformOutput',false);
    l3 = cellfun(@(c)strcmp(c,keyNdx),{'2'},'UniformOutput',false);
    typeIdx = any(l{1},2) & any(l2{1},2);
    baselineIdx = any(l{1},2) & any(l3{1},2);
    
    if strcmpi(state{currState},'stops')
        goodNdx_before = find(~cellfun(@isempty,s.optimal(:,baselineIdx)));
        goodNdx_during = find(~cellfun(@isempty,s.optimal(:,typeIdx)));
        
        turnBias_before = sign(cell2mat(s.optimal(goodNdx_before,baselineIdx)));
        turnBias_during = sign(cell2mat(s.optimal(goodNdx_during,typeIdx)));
        turnBias_before(turnBias_before==-1) = 0;
        turnBias_during(turnBias_during==-1) = 0;
%         turnBias_before = -sign(cell2mat(s.optimalStop(goodNdx_before,baselineIdx)));
%         turnBias_during = -sign(cell2mat(s.optimalStop(goodNdx_during,typeIdx)));
%         turnBias_before(turnBias_before==-1) = 0;
%         turnBias_during(turnBias_during==-1) = 0;
    else
        goodNdx_before = find(~cellfun(@isempty,s.optimal(:,baselineIdx)));
        goodNdx_during = find(~cellfun(@isempty,s.optimal(:,typeIdx)));
        turnBias_before = (cell2mat(s.optimal(goodNdx_before,baselineIdx)));
        turnBias_during = (cell2mat(s.optimal(goodNdx_during,typeIdx)));
    end

    allSpk = cell2mat(s.spkStart(goodNdx_during,typeIdx));
    alldSpk = cell2mat(s.dSpkStart(goodNdx_during,typeIdx));
    allTime = cell2mat(s.timeFromFe(goodNdx_during,typeIdx));
    
    turnOpt_during_baseline = turnBias_during(abs(allSpk-4.701)<0.001 & abs(alldSpk)<0.001);
    turnBias_during_baseline{currState} = sum(turnOpt_during_baseline)./numel(turnOpt_during_baseline);
    turnBias_during_spkRate = turnBias_during(abs(allSpk-4.701)>=0.001 | abs(alldSpk)>=0.001);
    alldSpk_spkRate = alldSpk(abs(allSpk-4.701)>=0.001 | abs(alldSpk)>=0.001);
    allSpk_spkRate = allSpk(abs(allSpk-4.701)>=0.001 | abs(alldSpk)>=0.001);
    allTime_time = allTime(abs(allSpk-4.701)>=0.001 | abs(alldSpk)>=0.001);
    
%     gprMd2 = fitrgp([alldSpk_spkRate, allSpk_spkRate],double(turnBias_during_spkRate),'Basis','linear',...
%       'FitMethod','exact','PredictMethod','exact');
%     [XX,YY] = meshgrid(xGrid,yGrid);
%     [ypred2,ystd2] = predict(gprMd2,[XX(:),YY(:)]);
%     figure;surf(XX,YY,reshape(ypred2,size(XX)));hold on;
%     surf(XX,YY,reshape(ypred2+ystd2,size(XX)));
%     surf(XX,YY,reshape(ypred2-ystd2,size(XX)));
%     zlim([0 1])
    
    [turnBias(:,currState),XX{currState},YY{currState}] = getTurnBias(turnBias_during_spkRate,alldSpk_spkRate, allSpk_spkRate,allTime_time,xGrid,yGrid,zGrid,ratio,Thresh,K);
    turnBias_baseline{currState} = sum(turnBias_before)./numel(turnBias_before);
    
    if plotFig
        [X,Y] = meshgrid(xGrid,yGrid);
        subplot(sbpltY,sbpltX,currState);
        surf(X,Y,turnBias{1,currState}(X,Y));hold on;
        points=[[min(xGrid) min(xGrid) max(xGrid) max(xGrid)]' ...
            [min(yGrid) max(yGrid) max(yGrid) min(yGrid)]' ...
            turnBias_baseline{currState}.*[1 1 1 1]'];
        points2=[[min(xGrid) min(xGrid) max(xGrid) max(xGrid)]' ...
            [min(yGrid) max(yGrid) max(yGrid) min(yGrid)]' ...
            turnBias_during_baseline{currState}.*[1 1 1 1]'];
        h = fill3(points(:,1),points(:,2),points(:,3),'k');
        h2 = fill3(points2(:,1),points2(:,2),points2(:,3),'r');
        grid on
        alpha(h, 0.1);alpha(h2, 0.1);zlim([0 1]);caxis([0 1])
        title([genotype ' ' state{currState} ' turn optimality (t=0)'])
        xlabel('df');ylabel('f');zlabel('Prob')
    end
end

end

function [Vq,XX,YY] = getTurnBias(turnBias_during,alldSpk, allSpk,allTime,xGrid,yGrid,zGrid,ratio,Thresh,K)
%---------
% baselineNdx = (allSpk<0.01) & (abs(alldSpk) < 1);
% turnBias_during = turnBias_during(~baselineNdx);
% X2 = [alldSpk(~baselineNdx), allSpk(~baselineNdx), allTime(~baselineNdx)]./ratio;
%---------

X2 = [alldSpk, allSpk, allTime]./ratio;
[XX,YY,ZZ] = ndgrid(xGrid,yGrid,zGrid);
Y2 = [XX(:),YY(:),ZZ(:)]./ratio;
[idx2, D2] = knnsearch(X2,Y2,'K',K);

% X2 = [alldSpk, allSpk]./ratio(1:2);
% [XX,YY] = meshgrid(xGrid,yGrid);
% Y2 = [XX(:),YY(:)]./ratio(1:2);
% [idx2, D2] = knnsearch(X2,Y2,'K',K);

turnBiasDat = double(turnBias_during(idx2));
mask2 = D2<Thresh;
turnBiasDat(~mask2) = nan;
turnBiasDat(sum(mask2,2)<15,:) = nan;
mask2(sum(mask2,2)<15,:) = 0;

turnBias = sum(turnBiasDat==1,2)./sum(mask2,2);
vals = reshape(turnBias,size(XX));

emptyZGrid = squeeze(sum(~isnan(vals),[1,2]))'<50;
notEmptyZGrid = find(~emptyZGrid);
emptyZGrid = find(emptyZGrid);
if ~isempty(emptyZGrid)
    [closestPtBef,closestPtAft,idxB] = findBeforeAfter(emptyZGrid,notEmptyZGrid,'both');
    dToNonEmptyGrid = abs([closestPtBef-emptyZGrid; closestPtAft-emptyZGrid]);
    [~,ndx] = min(dToNonEmptyGrid);
    closestZGrid = [];
    for k = 1:size(idxB,1)
        closestZGrid(k) = idxB(k,ndx(k));
    end
    closestZGrid = notEmptyZGrid(closestZGrid);
    
    vals(:,:,emptyZGrid) = vals(:,:,closestZGrid);
end

for slice = 1:numel(zGrid)
    [Vq{slice,1},~,~,~] = interpretLandscape(vals(:,:,slice)',XX(:,:,slice)',YY(:,:,slice)',2,false);    
end 
end


