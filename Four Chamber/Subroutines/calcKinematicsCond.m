function [s,keyNdx] = calcKinematicsCond(self,condNdx,key,thresh)

nKey = numel(key);nStates = numel(self.states.key);
bins = 1:nKey+1;
spd = self.calcSpd;
[~, thetaAll] = self.calcCurv;
thetaAll = thetaAll+pi;
[ang,dAng] = self.calcAng;

hist = ceil(0.2.*self.fs)-1;%200 ms history

[s] = initS(self.nFly,nStates*nKey);
%dSpkRateAll = gradient(self.spk).*self.fs./30;%%%%
dSpkRateAll = self.calcDeltaFR;

dFSmooth = dSpkRateAll(:,1:end-hist);
FSmooth = self.spk(:,1:end-hist);
for i = 2:hist+1
    dFSmooth = dSpkRateAll(:,i:end-(hist+1)+i)+dFSmooth;
    FSmooth = self.spk(:,i:end-(hist+1)+i)+FSmooth;
end
dFSmooth = [zeros(size(dFSmooth,1),hist),dFSmooth]./(hist+1);% average
FSmooth = [zeros(size(FSmooth,1),hist),FSmooth]./(hist+1);% average

for fly = 1:self.nFly
    [startNdx,endNdx,type] = startEndSeq(self.states.ndx(fly,:));
    startNdx(type == 0) = [];endNdx(type == 0) = [];type(type == 0) = [];
    startNdx = reshape(startNdx,[],1);endNdx = reshape(endNdx,[],1);
    for i = 1:nStates
        startNdxTmp = startNdx(type == i);
        endNdxTmp = endNdx(type == i);
        len = endNdxTmp-startNdxTmp+1;
        percContr = zeros(numel(startNdxTmp),nKey);
        for j = 1:numel(startNdxTmp)
            percContr(j,:) = histcounts(condNdx(fly,startNdxTmp(j):endNdxTmp(j)),bins)./len(j);
        end
        badTracks = ~any(percContr>thresh,2) | len<=3;
        startNdxTmp(badTracks) = [];
        endNdxTmp(badTracks) = [];
        percContr(badTracks,:) = [];
        [~,c] = max(percContr,[],2);
        
        kk = ones(1,nStates.*nKey);
        for j = 1:numel(startNdxTmp)
            k = (i-1).*nKey+c(j);
            
            spdTmp = spd(fly,startNdxTmp(j):endNdxTmp(j));
            s.avgSpd{fly,k}(kk(k),1) = mean(spdTmp);
            s.totSpd{fly,k}(kk(k),1) = sum(spdTmp);
            s.stdSpd{fly,k}(kk(k),1) = std(spdTmp);
            curvTmp = self.curv(fly,startNdxTmp(j):endNdxTmp(j));
            s.avgCurv{fly,k}(kk(k),1) = mean(curvTmp);
            s.totCurv{fly,k}(kk(k),1) = sum(curvTmp);
            s.maxCurv{fly,k}(kk(k),1) = max(abs(curvTmp));
            
            s.diffSpd{fly,k}(kk(k),1) = (spdTmp(end)-spdTmp(1))./numel(spdTmp);
            s.diffCurv{fly,k}(kk(k),1) = (curvTmp(end)-curvTmp(1))./numel(spdTmp);
            
            % previous state
            if startNdxTmp(j) == startNdx(1)
                prevStartNdx = max(startNdx(1),1);
                prevEndNdx = prevStartNdx+3;
                
                s.prevState{fly,k}(kk(k),1) = nan;
                s.currState{fly,k}(kk(k),1) = nan;
            else
                prevIndex = find(startNdx==startNdxTmp(j))-1;% previous state
                prevStartNdx = startNdx(prevIndex);
                prevEndNdx = startNdxTmp(j)-1;
                
                s.prevState{fly,k}(kk(k),1) = type(prevIndex);
                s.currState{fly,k}(kk(k),1) = type(prevIndex+1);
            end
            
            % next state
            if endNdxTmp(j) == endNdx(end)
                s.nextState{fly,k}(kk(k),1) = nan;
            else
                nextIndex = find(endNdx==endNdxTmp(j))+1;% previous state
                s.nextState{fly,k}(kk(k),1) = type(nextIndex);
            end
            
            curvTmp_prev = self.curv(fly,prevStartNdx:prevEndNdx);
            s.avgCurv_prev{fly,k}(kk(k),1) = mean(curvTmp_prev);
            s.totCurv_prev{fly,k}(kk(k),1) = sum(curvTmp_prev);
            
            %%%%
            startAng = ang(fly,startNdxTmp(j));
            endAng = ang(fly,endNdxTmp(j));
            startDir = [cosd(startAng);sind(startAng)];
            endDir = [cosd(endAng);sind(endAng)];
            startDirRelCent = -[self.x(fly,startNdxTmp(j)); self.y(fly,startNdxTmp(j))];%xH,yH
            endDirRelCent = -[self.x(fly,endNdxTmp(j)); self.y(fly,endNdxTmp(j))];
            
            n = [0;0;1];
            v = [startDir;0];
            u = [startDirRelCent;0];
            thetaStart = get2VecAngle(v,u,n,false);
            n = [0;0;1];
            v = [endDir;0];
            u = [endDirRelCent;0];
            thetaEnd = get2VecAngle(v,u,n,false);
            thetaDiff = thetaEnd-thetaStart;% negative means turned more inwards
            s.optimalStop{fly,k}(kk(k),1) = thetaDiff;
            %%%%
            
            % calculate current direction of movement relative to odor
            % source (+ is to the right, - is to the left)
            % if -, right is optimal dir; if +, left is optimal dir
            currAng = thetaAll(fly,max(prevEndNdx-hist,1):prevEndNdx);
            currAng = nanmean(currAng);
            currAng = wrapTo360(currAng.*180./pi)-90;
            prevDir = [cosd(currAng);sind(currAng)];
            dirRelCent = -[self.x(fly,prevEndNdx); self.y(fly,prevEndNdx)];
            n = [0;0;1];
            v = [prevDir;0];
            u = [dirRelCent;0];
            theta = get2VecAngle(v,u,n,false);
            
            % calculate direction of movement relative to odor source after
            % turning (+ is to the right, - is to the left)
            % if -, right is optimal dir; if +, left is optimal dir
            nextAng = thetaAll(fly,endNdxTmp(j):min(endNdxTmp(j)+1,self.nPt-1));
            nextAng = nanmean(nextAng);
            nextAng = wrapTo360(nextAng.*180./pi)-90;
            nextDir = [cosd(nextAng);sind(nextAng)];
            dirRelCent = -[self.x(fly,prevEndNdx); self.y(fly,prevEndNdx)];
            n = [0;0;1];
            v = [nextDir;0];
            u = [dirRelCent;0];
            theta2 = get2VecAngle(v,u,n,false);
            
            %--------------------------------------------------------------
            % code for turn optimality based on orientation
            %--------------------------------------------------------------
%             prevDir = [self.xH(fly,max(prevEndNdx-2,1):prevEndNdx)-...
%                 self.x(fly,max(prevEndNdx-2,1):prevEndNdx) ;...
%                     self.yH(fly,max(prevEndNdx-2,1):prevEndNdx)-...
%                     self.y(fly,max(prevEndNdx-2,1):prevEndNdx)];
                
            prevDir = [self.xH(fly,startNdxTmp(j):startNdxTmp(j)+2)-...
                self.x(fly,startNdxTmp(j):startNdxTmp(j)+2) ;...
                    self.yH(fly,startNdxTmp(j):startNdxTmp(j)+2)-...
                    self.y(fly,startNdxTmp(j):startNdxTmp(j)+2)];    
            
            prevDir = nanmean(prevDir,2);
            dirRelCent = -[self.x(fly,startNdxTmp(j)); self.y(fly,startNdxTmp(j))];
            %dirRelCent = -[self.x(fly,prevEndNdx); self.y(fly,prevEndNdx)];
            n = [0;0;1];
            v = [prevDir;0];
            u = [dirRelCent;0];
            theta_orient = get2VecAngle(v,u,n,false);
            %--------------------------------------------------------------
            
%             theta = nan(5,1);
%             for pt = 1:5
%                 try
%                     prevDir = [self.xH(fly,prevEndNdx-pt)-self.xH(fly,max(prevEndNdx-pt-2,1));
%                     self.yH(fly,prevEndNdx-pt)-self.yH(fly,max(prevEndNdx-pt-2,1))];
%                     dirRelCent = -[self.xH(fly,prevEndNdx); self.yH(fly,prevEndNdx)];
%                     n = [0;0;1];
%                     v = [prevDir;0];
%                     u = [dirRelCent;0];
%                     theta(pt) = get2VecAngle(v,u,n,false);
%                 catch
%                     a = 1;
%                 end
%             end
%             u = [-1, 0, 0];
%             v = [0, -1, 0];
%             theta = get2VecAngle(v,u,n,false);
%             theta(theta<0) = theta(theta<0)+360;
            theta3 = abs(theta);
            if theta3>180
                pause(0.1);
            end
            nonOpt = sign(nanmean(theta)) ~= sign(s.avgCurv{fly,k}(kk(k),1));
            theta3(nonOpt) = 360-theta3(nonOpt);
%             if sign(curvTmp)<0
%                 theta3 = 360-theta3;
%             end
            s.theta{fly,k}(kk(k),1) = theta3;


            % curvature (-) = right turn, (+) = left turn ==> optimal is
            % where curvature sign matches the current direction of
            % movement relative to odor
            s.optimal{fly,k}(kk(k),1) = sign(nanmean(theta)) == sign(s.avgCurv{fly,k}(kk(k),1));
            s.turnBias{fly,k}(kk(k),1) = abs(theta2)<abs(theta);
            
            s.optimalPLOS{fly,k}(kk(k),1) = sign(nanmean(theta_orient)) == sign(s.avgCurv{fly,k}(kk(k),1));
            
            %s.prev{fly,k}(kk(k),1:hist) = theta;
            
            
%             spkRate = self.spk(fly,max(startNdxTmp(j)-hist,1):min(endNdxTmp(j)+hist,self.nPt));
%             dSpkRate = dSpkRateAll(fly,max(startNdxTmp(j)-hist,1):min(endNdxTmp(j)+hist,self.nPt));%gradient(spkRate);
%             %midIdx = ceil(length(spkRate)./2)+1;
%             l = min(hist,length(spkRate));
%             s.spkStart{fly,k}(kk(k),1) = mean(spkRate(1:l));
%             %s.spkMid{fly,k}(kk(k),1) = mean(spkRate(midIdx-1:min(midIdx+1,length(spkRate))));
%             s.dSpkStart{fly,k}(kk(k),1) = mean(dSpkRate(1:l));
%             %s.dSpkMid{fly,k}(kk(k),1) = mean(dSpkRate(midIdx-1:min(midIdx+1,length(spkRate))));
%             s.spkEnd{fly,k}(kk(k),1) = mean(spkRate(end-l+1:end));
%             s.dSpkEnd{fly,k}(kk(k),1) = mean(dSpkRate(end-l+1:end));
            
            s.spkStart{fly,k}(kk(k),1) = FSmooth(fly,startNdxTmp(j));
            %s.spkMid{fly,k}(kk(k),1) = mean(spkRate(midIdx-1:min(midIdx+1,length(spkRate))));
            s.dSpkStart{fly,k}(kk(k),1) = dFSmooth(fly,startNdxTmp(j));
            %s.dSpkMid{fly,k}(kk(k),1) = mean(dSpkRate(midIdx-1:min(midIdx+1,length(spkRate))));
            s.spkEnd{fly,k}(kk(k),1) = FSmooth(fly,endNdxTmp(j));
            s.dSpkEnd{fly,k}(kk(k),1) = dFSmooth(fly,endNdxTmp(j));
            
%             if any(spkRate(1:end-hist)==0)
%                 prevTrack = max(prevEndNdx-hist,1):prevEndNdx;
%                 figure(100);
%                 plot(self.x(fly,prevTrack),self.y(fly,prevTrack),'b','Linewidth',2);hold on;
%                 plot(self.x(fly,startNdxTmp(j):endNdxTmp(j)),self.y(fly,startNdxTmp(j):endNdxTmp(j)),'r','Linewidth',2);
%                 q0 = quiver(self.x(fly,prevTrack(end)),self.y(fly,prevTrack(end)),-self.x(fly,prevTrack(end)),-self.y(fly,prevTrack(end)));
%                 q0.Color = 'blue';q0.MaxHeadSize = 1;q0.LineStyle = '--';
%                 q1 = quiver(self.x(fly,prevTrack(end)),self.y(fly,prevTrack(end)),cosd(currAng)./2,sind(currAng)./2);
%                 q2 = quiver(self.x(fly,endNdxTmp(j)),self.y(fly,endNdxTmp(j)),cosd(nextAng)./2,sind(nextAng)./2);
%                 q1.Color = 'blue';q1.MaxHeadSize = 1;q1.LineStyle = '--';
%                 q2.Color = 'red';q2.MaxHeadSize = 1;q2.LineStyle = '--';
%                 plotCircle([0,0],1,100,'r');
%                 plotCircle([0,0],0.01,100,'r');
%                 axis([-2 2 -2 2]);axis square
%                 hold off;
%                 title({['turnBias: ' num2str(s.turnBias{fly,k}(kk(k),1))...
%                     ' optimality: ' num2str(s.optimal{fly,k}(kk(k),1))...
%                     ' curv: ' num2str(s.totCurv{fly,k}(kk(k),1))...
%                     ' theta: ' num2str(s.theta{fly,k}(kk(k),1))], ...
%                     [' spkBefore: ' num2str(s.spkStart{fly,k}(kk(k),1))...
%                     ' spkAft: ' num2str(s.spkEnd{fly,k}(kk(k),1))]})
%                 
%                 pause(0.1);
%             end
            
            
            
            try
                s.locStart{fly,k}(kk(k),1) = mean(self.rH(fly,startNdxTmp(j):startNdxTmp(j)+1));%+l
            catch
                a;
            end
                                        %mod(abs(sum(abs(curvTmp))).*180./pi,360);
            %------------------
            %s.ndx{fly,k}(kk(k),1) = startNdxTmp(j)+find(abs(curvTmp)==max(abs(curvTmp)),1)-1;
            s.ndx{fly,k}(kk(k),1) = startNdxTmp(j)+1;
            %------------------
            s.dur{fly,k}(kk(k),1) = endNdxTmp(j)-startNdxTmp(j)+1;
            kk(k) = kk(k)+1;
        end
    end
end
[A,B] = meshgrid(self.states.key,key);
c=cat(2,A,B);
keyNdx=reshape(c,[],2);

% figure;hold on
% for i = 1:2
%     currRLoc = cell2mat(s.locStart(:,i));
%     xx =[0:0.1:4.1];
%     plot(histcounts(currRLoc,xx)./(pi.*(xx(2:end).^2-xx(1:end-1).^2)));
% end
% 
% figure;hold on
% for i = 1:1
%     currRLoc = cell2mat(s.locStart(:,i));
%     xx =[0:0.1:4.1];
%     tmp = histcounts(currRLoc,xx)./(pi.*(xx(2:end).^2-xx(1:end-1).^2));
%     plot(xx(1:end-1),tmp./sum(tmp));
%     xlabel('radial position');ylabel('Turn Density')
% end
% 
% s2 = s;
% xx2 = xx;
% load('tempRadLocTurns_1230', 's');
% 
% figure;hold on
% for i = 1:1
%     currRLoc = cell2mat(s.locStart(:,i));
%     xx =[0:0.1:4.1];
%     tmp = histcounts(currRLoc,xx)./(pi.*(xx(2:end).^2-xx(1:end-1).^2));
%     plot(xx(1:end-1),tmp./sum(tmp));hold on;
%     
%     currRLoc = cell2mat(s2.locStart(:,i));
%     tmp = histcounts(currRLoc,xx2)./(pi.*(xx2(2:end).^2-xx2(1:end-1).^2));
%     plot(xx2(1:end-1),tmp./sum(tmp));
%     
%     xlabel('radial position');ylabel('Turn Density')
% end

end


function [s] = initS(nFly,nCond)

s.avgSpd = cell(nFly,nCond);
s.totSpd = cell(nFly,nCond);
s.avgCurv = cell(nFly,nCond);
s.totCurv = cell(nFly,nCond);
s.ndx = cell(nFly,nCond);
s.dur = cell(nFly,nCond);

end