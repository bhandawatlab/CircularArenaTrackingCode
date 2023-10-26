function attractNdx = calcAttractionIndex(self,lab,border)
fe = getFirstEntry(self,lab,border);
for i = 1:self.nFly
    totInside = sum(self.(['r' lab])(i,fe(i):end)<border);
    totOutside = sum(self.(['r' lab])(i,fe(i):end)>=border);
    attractNdx.during(i) = totInside./(totInside+totOutside);
    attractNdx.duringIn(i) = totInside;
    attractNdx.duringTot(i) = totInside+totOutside;
    
    totInside = sum(self.(['r' lab])(i,1:fe(i))<border);
    totOutside = sum(self.(['r' lab])(i,1:fe(i))>=border);
    attractNdx.before(i) = totInside./(totInside+totOutside);
    attractNdx.beforeIn(i) = totInside;
    attractNdx.beforeTot(i) = totInside+totOutside;
end

end