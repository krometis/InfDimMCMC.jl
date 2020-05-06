# Assemble a chain split across multiple files
function assembleChain(file; verbose=true, sampleCols=:, obsCols=:, lpdfCols=:, sampleThin=1, obsThin=sampleThin, lpdfThin=sampleThin, arThin=sampleThin)
  files = [ file ];
  origfile = false;
  
  #get information from first file
  f = h5open(file,"r");
  #nsamp, sampdim = size(f["samples"][1:sampleThin:end,sampleCols]);
  #obsdim         = size(f["obs"][1:obsThin:end,obsCols],2);
  nsampAll, sampdim = size(f["samples"]);
  nsamp   = length(1:sampleThin:nsampAll);
  sampdim = length((1:sampdim)[sampleCols]);
  obsdim  = size(f["obs"],2);
  obsdim  = length((1:obsdim)[obsCols]);
  if any(x->x=="restartfile", names(f))
    file = read(f,"restartfile");
    files = [ file; files ];
  else
      origfile = true;
  end
  close(f);
  
  #find restart files
  while !origfile
    f = h5open(file,"r");
    nsampTmp = size(f["samples"],1);
    nsampAll += nsampTmp;
    nsamp    += length(1:sampleThin:nsampTmp);
    if any(x->x=="restartfile", names(f))
      file = read(f,"restartfile");
      files = [ file; files ];
    else
      origfile = true;
    end
    close(f);
  end
  
  samples = zeros(nsamp,sampdim);
  obs     = zeros(nsamp,obsdim);
  lpdfs   = zeros(nsamp,length((1:3)[lpdfCols]));
  ar      = zeros(Int64,nsamp);
  sampCnt = obsCnt = lpdfCnt = arCnt = 0;
  for file in files
    f = h5open(file,"r");
    #read samples
    n  = size(f["samples"],1);
    ns = length(1:sampleThin:n);
    samples[sampCnt+1:sampCnt+ns,:] = f["samples"][1:sampleThin:end,sampleCols];
    sampCnt += ns;
    #read obs
    nobs = length(1:obsThin:n);
    obs[obsCnt+1:obsCnt+nobs,:] = f["obs"][1:obsThin:end,obsCols];
    obsCnt += nobs;
    #read lpdf
    nlpdf = length(1:lpdfThin:n);
    lpdfs[lpdfCnt+1:lpdfCnt+nlpdf,:] = f["lpdfs"][1:lpdfThin:end,lpdfCols];
    lpdfCnt += nlpdf;
    #read ar 
    nar = length(1:arThin:n);
    ar[arCnt+1:arCnt+nar,:] = f["ar"][1:arThin:end];
    arCnt += nar;
    verbose && @printf("Read %d samples from %s\n",n,file);
  end
  verbose && @printf("Done. %d total samples read out of %d found.\n", nsamp, nsampAll);
  
  return samples, obs, lpdfs, ar, files;
end

# # Assemble a chain split across multiple files
# function assembleChain(file; verbose=true)
#   files = [ file ];
#   origfile = false;
#   
#   #get information from first file
#   f = h5open(file,"r");
#   nsamp, sampdim = size(f["samples"]);
#   obsdim         = size(f["obs"],2);
#   if any(x->x=="restartfile", names(f))
#     file = read(f,"restartfile");
#     files = [ file; files ];
#   else
#       origfile = true;
#   end
#   close(f);
#   
#   #find restart files
#   while !origfile
#     f = h5open(file,"r");
#     nsamp += size(f["samples"],1);
#     if any(x->x=="restartfile", names(f))
#       file = read(f,"restartfile");
#       files = [ file; files ];
#     else
#       origfile = true;
#     end
#     close(f);
#   end
#   
#   samples = zeros(nsamp,sampdim);
#   obs     = zeros(nsamp,obsdim);
#   lpdfs   = zeros(nsamp,3);
#   ar      = zeros(Int64,nsamp);
#   cnt = 0;
#   for file in files
#     sTemp = h5read(file,"samples");
#     samples[cnt+1:cnt+size(sTemp,1),:] = sTemp;
#     obs[cnt+1:cnt+size(sTemp,1),:]     = h5read(file,"obs");
#     lpdfs[cnt+1:cnt+size(sTemp,1),:]   = h5read(file,"lpdfs");
#     ar[cnt+1:cnt+size(sTemp,1),:]      = h5read(file,"ar");
#     cnt += size(sTemp,1);
#     verbose && @printf("Read %d samples from %s\n",size(sTemp,1),file);
#   end
#   verbose && @printf("Done. %d total samples found.\n", nsamp);
#   
#   return samples, obs, lpdfs, ar, files;
# end

