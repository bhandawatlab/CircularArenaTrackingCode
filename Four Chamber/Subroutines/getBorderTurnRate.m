function [f_orco] = getBorderTurnRate(f_orco,state,plotFig)
fprintf('Generating negative dF transition prob\n')

zGrid = f_orco.model.params{1, 1}.KNN.Time;
dz = f_orco.model.params{1, 1}.KNN.ratio(3);

% calculate smoothed change in firing rate
dF = f_orco.calcDeltaFR;
dF(:,end-1:end) = repmat(dF(:,end-2),1,2);
hist = ceil(0.2.*f_orco.fs)-1;% 200 ms sampling rate
dFSmooth = dF(:,1:end-hist);
FSmooth = f_orco.spk(:,1:end-hist);
for i = 2:hist+1
    dFSmooth = dF(:,i:end-(hist+1)+i)+dFSmooth;
    FSmooth = f_orco.spk(:,i:end-(hist+1)+i)+FSmooth;
end
dFSmooth = [zeros(size(dFSmooth,1),hist),dFSmooth]./(hist+1);% average
FSmooth = [zeros(size(FSmooth,1),hist),FSmooth]./(hist+1);% average

%%
ThreshOn = 15;
ThreshOff = -15;
OnState = dFSmooth>ThreshOn;
OffState = dFSmooth<ThreshOff;

t_off = nan(size(OffState));
t_on = nan(size(OffState));
for i = 2:f_orco.nPt
    t_off(:,i) = t_off(:,i-1)+1;
    t_off(OffState(:,i),i) = 0;
    t_on(:,i) = t_on(:,i-1)+1;
    t_on(OnState(:,i),i) = 0;
end
t_off(isnan(t_off)) = inf;
t_on(isnan(t_on)) = inf;

k = 1;nPt = 600;
t_trackOn = nan(1,nPt);
t_trackOff = nan(1,nPt);
for fly = 1:f_orco.nFly
    [startNdx,endNdx,type] = startEndSeq(f_orco.states.ndx(fly,:)==state);
    startNdx(type==0 | endNdx<5400) = [];
    endNdx(type==0 | endNdx<5400) = [];
    dur = endNdx-startNdx+1;
    for track = 1:numel(startNdx)
        t_trackOn(k,:) = [t_on(fly,startNdx(track):endNdx(track)) nan(1,nPt-dur(track))];
        t_trackOff(k,:) = [t_off(fly,startNdx(track):endNdx(track)) nan(1,nPt-dur(track))];
        k = k+1;
    end
end
badTracks = (max(t_trackOn,[],2)>1500) | (max(t_trackOff,[],2)>1500);
t_trackOn(badTracks,:) = [];
t_trackOff(badTracks,:) = [];

y = [];tOn = [];tOff = [];
for i = 1:size(t_trackOn,1)
    t_trackOnTmp = t_trackOn(i,~isnan(t_trackOn(i,:)));
    t_trackOffTmp = t_trackOff(i,~isnan(t_trackOn(i,:)));
    y = [y;zeros(numel(t_trackOnTmp)-1,1);1];
    tOn = [tOn;t_trackOnTmp'];
    tOff = [tOff;t_trackOffTmp'];
end
tOn = tOn./f_orco.fs;%keep it in per second
tOff = tOff./f_orco.fs;%keep it in per second

[NTrans,xBin,yBin] = histcounts2(tOn(y==1),tOff(y==1),[0:5:1505]./f_orco.fs,[0:5:1505]./f_orco.fs);
[NAll,xBin,yBin] = histcounts2(tOn,tOff,[0:5:1505]./f_orco.fs,[0:5:1505]./f_orco.fs);
figure;set(gcf,'Position',[2 42 838 924])
subplot(2,2,1);scatter(tOn,tOff);xlabel('Time (s) since high df');ylabel('Time (s) since low df')
subplot(2,2,2);imagesc(xBin(1:end-1),yBin(1:end-1),log10(NAll))
set(gca,'YDir','normal');colorbar;colormap(hot);
xlabel('Time (s) since high df');ylabel('Time (s) since low df');title('density of all CW tracks (log scale)')
subplot(2,2,3);scatter(tOn(y==1),tOff(y==1));xlabel('Time (s) since high df');ylabel('Time (s) since low df')
subplot(2,2,4);imagesc(xBin(1:end-1),yBin(1:end-1),log10(NTrans))
set(gca,'YDir','normal');colorbar;colormap(hot);
xlabel('Time (s) since high df');ylabel('Time (s) since low df');title('density of transitions (log scale)')

%%
fe = f_orco.getFirstEntry('H',1.3);
%Thresh = 15;currentState = 'on';
Thresh = -15;currentState = 'off';
if strcmpi(currentState,'on')
    OffState = dFSmooth>Thresh;
else
    OffState = dFSmooth<Thresh;
end
rndBefState = rand(size(dFSmooth))>0.85;
rndBefState(:,5400:end) = false;

t = nan(size(OffState));t_bef = nan(size(rndBefState));
tSinceFE = nan(size(OffState));
for i = 2:f_orco.nPt
    t(:,i) = t(:,i-1)+1;
    t(OffState(:,i),i) = 0;
    tSinceFE(:,i) = i-fe;
    
    t_bef(:,i) = t_bef(:,i-1)+1;
    t_bef(rndBefState(:,i),i) = 0;
end
t(isnan(t)) = inf;
t_bef(isnan(t)) = inf;

k = 1;nPt = 600;
tBef_track = nan(1,nPt);
for fly = 1:f_orco.nFly
    [startNdx,endNdx,type] = startEndSeq(f_orco.states.ndx(fly,:)==state);
    startNdx(type==0 | endNdx>5400) = [];
    endNdx(type==0 | endNdx>5400) = [];
    dur = endNdx-startNdx+1;
    for track = 1:numel(startNdx)
        tBef_track(k,:) = [t_bef(fly,startNdx(track):min(endNdx(track),startNdx(track)+nPt-1))...
            nan(1,nPt-dur(track))];
        k = k+1;
    end
end
badTracks = min(tBef_track,[],2)>600;
tBef_track(badTracks,:) = [];


k = 1;nPt = 600;
t_track = nan(1,nPt);
t_trackAll = nan(1,nPt);
t_fe = nan(1,nPt);
for fly = 1:f_orco.nFly
    [startNdx,endNdx,type] = startEndSeq(f_orco.states.ndx(fly,:)==state);
    startNdx(type==0 | endNdx<5400) = [];
    endNdx(type==0 | endNdx<5400) = [];
    dur = endNdx-startNdx+1;
    for track = 1:numel(startNdx)
        t_track(k,:) = [t(fly,startNdx(track):endNdx(track)) nan(1,nPt-dur(track))];
        t_trackAll(k,:) = [t(fly,startNdx(track):min(startNdx(track)+nPt-1,f_orco.nPt)),...
             nan(1,(startNdx(track)+nPt-1)-f_orco.nPt)];
         
         t_fe(k,:) = [tSinceFE(fly,startNdx(track):endNdx(track)) nan(1,nPt-dur(track))];
        k = k+1;
    end
end
badTracks = min(t_track,[],2)>600;

t_track(badTracks,:) = [];
t_trackAll(badTracks,:) = [];
t_fe(badTracks,:) = [];

%%
[y,t,tSinceFE] = getTransitionPts(f_orco,t_track,t_fe);
[yBef,tBef,~] = getTransitionPts(f_orco,tBef_track,[]);

x = [1:1:600]./f_orco.fs;
xx = linspace(0,20,100);
NTrans = histcounts(tBef(yBef==1),[0:1:601]./f_orco.fs);
NAll = histcounts(tBef,[0:1:601]./f_orco.fs);
peakBef = NTrans(1)./NAll(1);
[f0Bef,f1Bef,f2Bef,ss_NTransBef] = fitTransitionProb(f_orco,NTrans,NAll,x,xx);
f_orco.model.borderChoice.peakBef(1,state) = peakBef;
f_orco.model.borderChoice.fBef{state} = @(x) f2Bef(x)./f1Bef(x);

for slice = 1:numel(zGrid)
    %display(['slice:'  num2str(slice)])
    tSinceFETmp = tSinceFE((tSinceFE>zGrid(slice)-dz) & (tSinceFE<zGrid(slice)+dz));
    yTmp = y((tSinceFE>zGrid(slice)-dz) & (tSinceFE<zGrid(slice)+dz));
    tTmp = t((tSinceFE>zGrid(slice)-dz) & (tSinceFE<zGrid(slice)+dz));

    NTrans = histcounts(tTmp(yTmp==1),[0:1:601]./f_orco.fs);
    NAll = histcounts(tTmp,[0:1:601]./f_orco.fs);
    
	[f0,f1,f2,ss_NTrans] = fitTransitionProb(f_orco,NTrans,NAll,x,xx);

    if plotFig
        figure;set(gcf,'Position',[2 42 838 924])
        subplot(3,1,1);plot(x,NTrans(2:end)-ss_NTrans,'o',xx,f2(xx),'r-','Linewidth',2);
        subplot(3,1,2);plot(x,NAll(2:end),'o',xx,f1(xx),'r-','Linewidth',2);
        subplot(3,1,3);plot([0 xx],[NTrans(1)./NAll(1); f2(xx)./(f1(xx))],'r-','Linewidth',2);hold on;
        plot([0 xx],[peakBef; f2Bef(xx)./(f1Bef(xx))],'g-','Linewidth',2);
        ylim([0 0.1]);
        suptitle([num2str(Thresh) ', ' currentState])
    end
    f_orco.model.borderChoice.peak(slice,state) = NTrans(1)./NAll(1);
    f_orco.model.borderChoice.f{slice,state} = @(x) f2(x)./f1(x);
    f_orco.model.borderChoice.state{state} = f_orco.states.key{state};
    f_orco.model.borderChoice.Thresh(state) = Thresh;
    f_orco.model.borderChoice.type{state} = currentState;
    
end
f_orco.model.borderChoice.tt = zGrid./f_orco.fs;
f_orco.model.borderChoice.dt = dz./f_orco.fs;

% for i = 1:2
%     figure(i)
%     print('-opengl','-dpsc2','BorderChoice040621.ps','-loose','-append');
% end
% ps2pdf('psfile', 'BorderChoice040621.ps', 'pdffile', 'BorderChoice040621.pdf', 'gspapersize', 'letter',...
% 'gscommand','C:\Program Files\gs\gs9.50\bin\gswin64.exe',...
% 'gsfontpath','C:\Program Files\gs\gs9.50\lib',...
% 'gslibpath','C:\Program Files\gs\gs9.50\lib');
end

function [y,t,tSinceFE] = getTransitionPts(f_orco,t_track,t_fe)

y = [];t = [];tSinceFE = [];
for i = 1:size(t_track,1)
    t_trackTmp = t_track(i,~isnan(t_track(i,:)));
    if ~isempty(t_trackTmp)
        y = [y;zeros(numel(t_trackTmp)-1,1);1];
        t = [t;t_trackTmp'];
        if ~isempty(t_fe)
            t_trackTmpFE = t_fe(i,~isnan(t_fe(i,:)));
            tSinceFE = [tSinceFE;t_trackTmpFE'];
        end
    end
end
t = t./f_orco.fs;%keep it in per second

end

function [f0,f1,f2,ss_NTrans] = fitTransitionProb(f_orco,NTrans,NAll,x,xx)
ss_NTrans = sum(NTrans(end-f_orco.fs+1:end))./f_orco.fs;
ss_NAll = sum(NAll(end-f_orco.fs+1:end))./f_orco.fs;

g1 = fittype('b*exp(-c*x)');
g2 = fittype('a+b*exp(-c*x)');
% f0 = fit(x',NTrans(2:end)',g,'StartPoint',[0;30;2]);
% f1 = fit(x',NAll(2:end)',g,'StartPoint',[20;350;2]);
%     f0 = fit(x',NTrans(2:end)',g2,'StartPoint',[0;30;2]);
%     f1 = fit(x',NAll(2:end)',g2,'StartPoint',[0;350;2]);
f0 = fit(x',NTrans(2:end)',g2,'StartPoint',[ss_NTrans+eps;NTrans(2);1]);
f1 = fit(x',NAll(2:end)',g2,'StartPoint',[ss_NAll+eps;NAll(2);1]);
f2 = fit(xx',(f0(xx)-f0.a),g1,'StartPoint',[f0.b;f0.c]);

ss_NTrans = f0.a;
end



