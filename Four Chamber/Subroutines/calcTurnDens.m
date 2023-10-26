function turnDens = calcTurnDens(self,lab,s,key)
dx = [0:0.025:1];
dx(end) = 1.1;

STNdx = any(strcmpi('sharp turns', key),2);
turns = s.ndx(:,STNdx);
cond = key(STNdx,2);nCond = numel(cond);
befNdx = strcmpi(cond,'before');

tmpBef = [];tmpDur = [];
for fly = 1:self.nFly
    tmpBef = [tmpBef self.(['r' lab])(fly,turns{fly,befNdx})];
    tmpDur = [tmpDur self.(['r' lab])(fly,cell2mat(reshape(turns(fly,~befNdx),[],1)))];
end
tmpBef = tmpBef./self.rBound;tmpDur = tmpDur./self.rBound;
decisionDensBef = histcounts(tmpBef,dx)./(dx(2:end).^2 - dx(1:end-1).^2);
turnDens.Bef = decisionDensBef./sum(decisionDensBef);
decisionDensDur = histcounts(tmpDur,dx)./(dx(2:end).^2 - dx(1:end-1).^2);
turnDens.Dur = decisionDensDur./sum(decisionDensDur);
turnDens.x = dx;

end