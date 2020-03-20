[readframe,nframes,fid,headerinfo] = get_readframe_fcn('180925_1_3_video.ufmf');
[mu,sig] = compute_bg_mean(readframe,1,5400,varargin);
tic
for i = 1:1000
    img = readframe(i);
end
toc