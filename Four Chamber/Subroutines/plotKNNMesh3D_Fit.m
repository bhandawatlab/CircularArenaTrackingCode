function [fitFun,pHat,XX,YY,pHatBase,CorrMat1,CorrMat2,qBounds1,qBounds2,...
    baselineDat,allDat,allSpk,alldSpk,allTime,allTime2,allFly,allFlyBaseline,turnOpt_during] = plotKNNMesh3D_Fit(f_orco, s, keyNdx, xGrid, yGrid, ...
    zGrid, ratio, Thresh, K, state, kin, genotype, plotFig)%f_orco, condNdx

fe = f_orco.getFirstEntry('H',1.3);

% calculate time from first entry for the beginning of each state
for fly = 1:size(s.ndx,1)
    for j = 1:size(s.ndx,2)
        s.timeFromFe{fly,j} = s.ndx{fly,j}-fe(fly);
    end
end

l = cellfun(@(c)strcmp(c,keyNdx),{state},'UniformOutput',false);
l2 = cellfun(@(c)strcmp(c,keyNdx),{'1'},'UniformOutput',false);
l3 = cellfun(@(c)strcmp(c,keyNdx),{'2'},'UniformOutput',false);
typeIdx = any(l{1},2) & any(l2{1},2);
baselineIdx = any(l{1},2) & any(l3{1},2);

if strcmpi(kin, 'spd')
    baselineDat = (cell2mat(s.avgSpd(:,baselineIdx)));
    allDat = abs(cell2mat(s.avgSpd(:,typeIdx)));
%     baselineDat = (cell2mat(s.diffSpd(:,baselineIdx)));
%     allDat = (cell2mat(s.diffSpd(:,typeIdx)));
    zlab = 'Speed (mm/s)';
    x = [0:1:30];
    %fitFun = 'norm';
    fitFun = 'logN';
elseif strcmpi(kin, 'avgCurv')
    baselineDat = abs(cell2mat(s.avgCurv(:,baselineIdx))*180./pi).*f_orco.fs;%deg/s
    allDat = abs(cell2mat(s.avgCurv(:,typeIdx))*180./pi).*f_orco.fs;%deg/s
%     baselineDat = abs(cell2mat(s.maxCurv(:,baselineIdx))*180./pi);%deg/s
%     allDat = abs(cell2mat(s.maxCurv(:,typeIdx))*180./pi);%deg/s
%     baselineDat = abs(cell2mat(s.stdSpd(:,baselineIdx)));%deg/s
%     allDat = abs(cell2mat(s.stdSpd(:,typeIdx)));%deg/s
    zlab = 'Avg Curv';
    x = [0:0.5:10];
    fitFun = 'logN';
elseif strcmpi(kin, 'totCurv')
    baselineDat = abs(cell2mat(s.totCurv(:,baselineIdx))*180./pi);
    allDat = abs(cell2mat(s.totCurv(:,typeIdx))*180./pi);
    baselineDat = mod(baselineDat,360);
    allDat = mod(allDat,360);
    badNdx = allDat<0.01;allDat(badNdx) = 0.01;
    zlab = 'Tot Curv';
    x = [0:10:360];
    fitFun = 'logN';
elseif strcmpi(kin, 'dur')
%     baselineDat = abs(cell2mat(s.dur(:,baselineIdx)))./f_orco.fs;% in seconds
%     allDat = abs(cell2mat(s.dur(:,typeIdx)))./f_orco.fs;% in seconds
%     zlab = 'Duration relative to before';
%     x = 0:0.2:30;
%     fitFun = 'exp';
    baselineDat = abs(cell2mat(s.dur(:,baselineIdx)))./f_orco.fs;% in seconds
    allDat = abs(cell2mat(s.dur(:,typeIdx)))./f_orco.fs;% in seconds
    zlab = 'Duration relative to before';
    x = 0:0.2:30;
    fitFun = 'logN';
end
%---------------------------
goodNdx_before = find(~cellfun(@isempty,s.optimal(:,baselineIdx)));
goodNdx_during = find(~cellfun(@isempty,s.optimal(:,typeIdx)));
turnOpt_before = (cell2mat(s.optimal(goodNdx_before,baselineIdx)));
turnOpt_during = (cell2mat(s.optimal(goodNdx_during,typeIdx)));

turnBias_before = (cell2mat(s.turnBias(goodNdx_before,baselineIdx)));
turnBias_during = (cell2mat(s.turnBias(goodNdx_during,typeIdx)));

f_After = (cell2mat(s.spkEnd(goodNdx_during,typeIdx)));
f_Before = (cell2mat(s.spkStart(goodNdx_during,typeIdx)));
df_After = (cell2mat(s.dSpkEnd(goodNdx_during,typeIdx)));
df_Before = (cell2mat(s.dSpkStart(goodNdx_during,typeIdx)));

allSpk = cell2mat(s.spkStart(:,typeIdx));
alldSpk = cell2mat(s.dSpkStart(:,typeIdx));
idx = allSpk<0.01 & abs(alldSpk)<0.001;

tmpTurnBias = sum(turnBias_during&idx)./sum(idx);
tmpTurnOpt = sum(turnOpt_during&idx)./sum(idx);
tmpTurnHighF = sum(f_After>f_Before & idx)./sum(idx);
tmpTurnNCF = sum(f_After==f_Before & idx)./sum(idx);
tmpTurnLowF = sum(f_After<f_Before & idx)./sum(idx);
tmpTurnHighdF = sum(df_After>0 & idx)./sum(idx);
tmpTurnNCdF = sum(df_After==0 & idx)./sum(idx);
tmpTurnLowdF = sum(df_After<0 & idx)./sum(idx);

if (strcmpi(kin,'totCurv') || strcmpi(kin,'avgCurv')) && plotFig
    if strcmpi(kin,'totCurv')
        maxCurv = 360;ylab = 'curv (deg)';
    elseif strcmpi(kin,'avgCurv')
        maxCurv = 90;ylab = 'curv (deg/s)';
    end
    
    theta_during = (cell2mat(s.theta(goodNdx_during,typeIdx)));
    figure;set(gcf,'position',[849 49 824 918]);subplot(3,1,1);
    scatter(theta_during(idx),allDat(idx));hold on;
    scatter(theta_during(turnBias_during & idx),allDat(turnBias_during & idx),'r*');
    scatter(theta_during(~turnBias_during & idx),allDat(~turnBias_during & idx),'g*');
    plot([180 180],[0 maxCurv]);
    text(30,maxCurv/360*250,'Optimal');
    text(210,maxCurv/360*250,'Non-optimal')
    xlim([0 360]);ylim([0 maxCurv])
    xlabel('theta');ylabel(ylab);
    title(['turnBias: ' num2str(tmpTurnBias),' optimality: ' num2str(tmpTurnOpt)])
    legend({'all inhibition turns','turn inward turns','turn outward turns'})

    subplot(3,1,2);scatter(theta_during(idx),allDat(idx));hold on;
    scatter(theta_during(f_After>f_Before & idx),allDat(f_After>f_Before & idx),'r*');
    scatter(theta_during(f_After==f_Before & idx),allDat(f_After==f_Before & idx),'g*');
    scatter(theta_during(f_After<f_Before & idx),allDat(f_After<f_Before & idx),'c*');
    plot([180 180],[0 maxCurv]);
    text(30,maxCurv/360*250,'Optimal');
    text(210,maxCurv/360*250,'Non-optimal')
    xlim([0 360]);ylim([0 maxCurv])
    xlabel('theta');ylabel(ylab);
    title(['higher spk rate: ' num2str(tmpTurnHighF),' NC spk rate: ' ...
        num2str(tmpTurnNCF),' lower spk rate: ' num2str(tmpTurnLowF)])
    legend({'all inhibition turns','ends in higher firing rate',...
        'ends in same firing rate','ends in lower firing rate'})

    subplot(3,1,3);scatter(theta_during(idx),allDat(idx));hold on;
    scatter(theta_during(df_After>0 & idx),allDat(df_After>0 & idx),'r*');
    scatter(theta_during(df_After==0 & idx),allDat(df_After==0 & idx),'g*');
    scatter(theta_during(df_After<0 & idx),allDat(df_After<0 & idx),'c*');
    plot([180 180],[0 maxCurv]);
    text(30,maxCurv/360*250,'Optimal');
    text(210,maxCurv/360*250,'Non-optimal')
    xlim([0 360]);ylim([0 maxCurv])
    xlabel('theta');ylabel(ylab);
    title(['+ delta spk rate: ' num2str(tmpTurnHighdF),' NC delta spk rate: ' ...
        num2str(tmpTurnNCdF),' - delta spk rate: ' num2str(tmpTurnLowdF)])
    legend({'all inhibition turns','ends in increasing firing rate',...
        'ends in no change in firing rate','ends in decreasing firing rate'})
    suptitle([genotype ', ' state ', ' kin ', K=' num2str(K) ', Thresh=' num2str(Thresh)])
end
%---------------------------

%---------------------------
if strcmpi(state,'sharp turns')
    spd1 = (abs(cell2mat(s.avgSpd(:,baselineIdx))));
    %spd1 = ((cell2mat(s.diffSpd(:,baselineIdx))));
    curv1 = (abs(cell2mat(s.totCurv(:,baselineIdx))*180./pi));
    curv1 = mod(curv1,360);
    dur1 = (abs(cell2mat(s.dur(:,baselineIdx)))./f_orco.fs);
    
    spd2 = (abs(cell2mat(s.avgSpd(:,typeIdx))));
    curv2 = (abs(cell2mat(s.totCurv(:,typeIdx))*180./pi));
    curv2 = mod(curv2,360);
    dur2 = (abs(cell2mat(s.dur(:,typeIdx)))./f_orco.fs);
    
elseif strcmpi(state,'curved walks')
    spd1 = (abs(cell2mat(s.avgSpd(:,baselineIdx))));
    %spd1 = ((cell2mat(s.diffSpd(:,baselineIdx))));
    curv1 = (abs(cell2mat(s.avgCurv(:,baselineIdx))*180./pi));
    dur1 = (abs(cell2mat(s.dur(:,baselineIdx)))./f_orco.fs);
    
    spd2 = (abs(cell2mat(s.avgSpd(:,typeIdx))));
    curv2 = (abs(cell2mat(s.avgCurv(:,typeIdx))*180./pi));
    dur2 = (abs(cell2mat(s.dur(:,typeIdx)))./f_orco.fs);
else
    curv1 = (abs(cell2mat(s.totCurv(:,baselineIdx))*180./pi));
    curv1 = mod(curv1,360);
    dur1 = (abs(cell2mat(s.dur(:,baselineIdx)))./f_orco.fs);
    
    curv2 = (abs(cell2mat(s.totCurv(:,typeIdx))*180./pi));
    curv2 = mod(curv2,360);
    dur2 = (abs(cell2mat(s.dur(:,typeIdx)))./f_orco.fs);
end
q1 = 0.001;q2 = 0.999;
%q1 = 0.01;q2 = 0.99;
if ~strcmpi(state,'stops')
    
    CorrMat1 = corrcoef([spd1,curv1,dur1]);
    qBounds1 = [quantile(spd1,q1)  quantile(spd1,q2);
        quantile(curv1,q1)  quantile(curv1,q2);
        quantile(dur1,q1)  quantile(dur1,q2)];
    
    CorrMat2 = corrcoef([spd2,curv2,dur2]);
    qBounds2 = [quantile(spd2,q1)  quantile(spd2,q2);
        quantile(curv2,q1)  quantile(curv2,q2);
        quantile(dur2,q1)  quantile(dur2,q2)];
else
    CorrMat1 = corrcoef([curv1,dur1]);
    qBounds1 = [quantile(curv1,q1)  quantile(curv1,q2);
        quantile(dur1,q1)  quantile(dur1,q2)];
    
    CorrMat2 = corrcoef([curv2,dur2]);
    qBounds2 = [quantile(curv2,q1)  quantile(curv2,q2);
        quantile(dur2,q1)  quantile(dur2,q2)];
end
%---------------------------

% baselineAvg = mean(baselineDat);
% stdBase = std(baselineDat);
%allDat = allDat-baselineAvg;

locDat = cell2mat(s.locStart(:,typeIdx));

allSpk = cell2mat(s.spkStart(:,typeIdx));
alldSpk = cell2mat(s.dSpkStart(:,typeIdx));
allTime = max(0,cell2mat(s.timeFromFe(:,typeIdx)));
allTime2 = max(0,cell2mat(s.ndx(:,typeIdx)));
allFly = repelem([1:f_orco.nFly],cellfun(@numel, s.ndx(:,typeIdx), 'UniformOutput',true))';
allFlyBaseline = repelem([1:f_orco.nFly],cellfun(@numel, s.ndx(:,baselineIdx), 'UniformOutput',true))';

% baselineNdx = (allSpk<0.01) & (abs(alldSpk) < 1);
% allDat2 = allDat(~baselineNdx);

% [F,x] = ecdf(baselineDat);figure;plot(x,F,'k','Linewidth',1);hold on;
% [F,x] = ecdf(allDat(allSpk<4.7));plot(x,F,'r','Linewidth',1);
% [F,x] = ecdf(allDat(alldSpk<-10));plot(x,F,'g','Linewidth',1);
% [F,x] = ecdf(allDat(alldSpk>10));plot(x,F,'b','Linewidth',1);
% hold off;
% legend({'baseline','inhibition','df<-10','df>10'})
% xlabel('speed (mm/s)');ylabel('cdf');title('std speed in curved walk state')

% 3D
%--------------------------------------------------------------------------
%X = [alldSpk(~baselineNdx), allSpk(~baselineNdx), allTime(~baselineNdx)]./ratio;
X = [alldSpk, allSpk, allTime]./ratio;
[XX,YY,ZZ] = ndgrid(xGrid,yGrid,zGrid);
Y = [XX(:),YY(:),ZZ(:)]./ratio;
[idx, D] = knnsearch(X,Y,'K',K);
mask = D<Thresh;

maskedDat = allDat(idx);
maskedDat(~mask) = nan;

try
    if strcmpi(fitFun, 'logN')
        pHat = nan(size(maskedDat,1),2);
        p = nan(size(maskedDat,1),2);
        for i = 1:size(maskedDat,1)%6
            if sum(mask(i,:))>15
                dat = maskedDat(i,:);
                dat(isnan(dat))= [];
                pHat(i,:) = lognfit(dat+eps);
            end
        end
        pHatBase = lognfit(baselineDat+eps);
        
    elseif strcmpi(fitFun, 'exp')
        pHat = nan(size(maskedDat,1),1);
        for i = 1:size(maskedDat,1)
            if sum(mask(i,:))>15
                dat = maskedDat(i,:);
                dat(isnan(dat))= [];
                pHat(i,:) = expfit(dat);
            end
        end
        pHatBase = expfit(baselineDat);
    elseif strcmpi(fitFun, 'norm')
        pHat = nan(size(maskedDat,1),2);
        p = nan(size(maskedDat,1),2);
        for i = 1:size(maskedDat,1)%6
            if sum(mask(i,:))>15
                dat = maskedDat(i,:);
                dat(isnan(dat))= [];
                [pHat(i,1),pHat(i,2)] = normfit(dat);
            end
        end
        pHatBase = normfit(baselineDat);
    end
catch
    disp('error with fitFun name, use logN, exp, or norm');
end
end