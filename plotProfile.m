function y = plotProfile(trenchHeight,trenchLWall,trenchRWall)

chunksize = 200;
chunknum = trenchHeight/chunksize;
result = zeros(chunknum,1);

for i = 1:chunknum
    chunkL = trenchLWall(((i-1)*chunksize)+1:(i*chunksize));
    chunkR = trenchRWall(((i-1)*chunksize)+1:(i*chunksize));
    chunk = (chunkL + chunkR)/2;
    result(chunknum+1-i,1) = sum(chunk);
end

y = result/chunksize;