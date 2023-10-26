function [self,stateDistBefore,stateDistEarly,stateDistLate] = getBaselineKinematics(self,condNdx,key)
%close all
spdAll = self.calcSpd;

% separate the track into combinations of state and condition
[condState, stateKey] = self.getCond(condNdx,key);
nCondStates = size(stateKey,1);
nState = numel(self.states.key);

% get the index of each state (i.e. sharp turn/curved walk)
stateNdx = cell(nState,1);
for i = 1:nState
    stateNdx{i} = cellfun(@(x) strcmp(x,self.states.key{i}),stateKey(:,1),'UniformOutput',true);
end
% get the index for before and baseline spiking during
before = cellfun(@(x) strcmp(x,'before'),stateKey(:,2),'UniformOutput',true);
baseline = cellfun(@(x) strcmp(x,'baseline'),stateKey(:,2),'UniformOutput',true);

% get the raw values of kinematics for baseline (no spiking) kinematics for
% all conditions and states
% -------------------------------------------------------------------------
dur = cell(nCondStates,1);
spd = cell(nCondStates,1);
avgCurv = cell(nCondStates,1);
totCurv = cell(nCondStates,1);
dur2 = cell(numel(key),1);
for i = 1:self.nFly
    [startNdx,endNdx,type] = startEndSeq(condState(i,:));
    for j = 1:max(type)
        dur{j} = [dur{j} endNdx(type == j)-startNdx(type == j)+1];
        
        currTrack = find(type == j);
        for k = 1:numel(currTrack)
            kk = currTrack(k);
            spd{j} = [spd{j} nanmean(spdAll(i,startNdx(kk):endNdx(kk)))];
            avgCurv{j} = [avgCurv{j} nanmean(self.curv(i,startNdx(kk):endNdx(kk)))];
            totCurv{j} = [totCurv{j} nansum(self.curv(i,startNdx(kk):endNdx(kk)))];
        end
    end
    
    [startNdx,endNdx,type] = startEndSeq(self.states.ndx(i,:));
    % dur2 is just a checking function for when I was debugging
    for j = 1:numel(key)
        dur2{j} = [dur2{j} endNdx(type == j)-startNdx(type == j)+1];
    end
end

% save the raw values of kinematics for baseline (no spiking) kinematics for
% before and baseline spiking during conditions for all states
stateDistBefore = cell(4,1);
stateDistEarly = cell(4,1);
for i = 1:nState
    stateDistBefore{i,1}.dur = dur{stateNdx{i}&before}./self.fs;
    stateDistBefore{i,1}.spd = spd{stateNdx{i}&before};
    stateDistBefore{i,1}.avgCurv = abs(avgCurv{stateNdx{i}&before}).*180./pi./self.fs;
    stateDistBefore{i,1}.totCurv = abs(totCurv{stateNdx{i}&before}).*180./pi;
    
    stateDistEarly{i,1}.dur = dur{stateNdx{i}&baseline}./self.fs;
    stateDistEarly{i,1}.spd = spd{stateNdx{i}&baseline};
    stateDistEarly{i,1}.avgCurv = abs(avgCurv{stateNdx{i}&baseline}).*180./pi./self.fs;
    stateDistEarly{i,1}.totCurv = abs(totCurv{stateNdx{i}&baseline}).*180./pi;
end
% not differentiating between early and late (no adaptation/habituation)
stateDistLate = stateDistEarly;

% check to make sure all conditions are present
assert(all(sum(reshape(cellfun(@(x) numel(x),dur),4,4))' == cellfun(@(x) numel(x),dur2)));

% fit raw values to distribution and save to the current fly object if the
% model model parameters have already been initialized
if ~isempty(self.model)
    for i = 1:numel(self.model.params)
        kin = self.model.params{i}.state.kin;
        state = self.model.params{i}.state.state;
        
        l = cellfun(@(c)strcmp(c,self.states.key),{state},'UniformOutput',false);
        ndx = find(l{1});
        
        % early baseline after first entry
        x = stateDistEarly{ndx}.(kin);
%         try
        parmhat =  lognfit(x+eps);
%         catch
%             a = 1;
%         end
        self.model.params{i}.baselineEarly = parmhat;
        
        % later baseline after first entry
        x = stateDistLate{ndx}.(kin);
        parmhat =  lognfit(x+eps);
        self.model.params{i}.baselineLate = parmhat;
    end
end

end