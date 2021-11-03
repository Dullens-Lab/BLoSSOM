function [dbn2] = rejectOutlier(dbn,it)


dbn2 = dbn;
for a = 1:it
    mdist = mean(dbn2);
    stddist = std(dbn2);
    dbn2 = dbn2(dbn2 > mdist - 2* stddist & dbn2 <...
        mdist + 2* stddist);
end 
end

