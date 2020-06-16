# Assemble a chain split across multiple files
function assembleChain2(file; verbose=true, sampleCols=:, obsCols=:, lpdfCols=:, sampleThin=1, obsThin=sampleThin, lpdfThin=sampleThin, arThin=sampleThin)
  files = [ file ];
  origfile = false;
  
  #get information from first file
  f = h5open(file,"r");
  #nsamp, sampdim = size(f["samples"][1:sampleThin:end,sampleCols]);
  #obsdim         = size(f["obs"][1:obsThin:end,obsCols],2);
  nsampAll, sampdim = size(f["samples"]);

  #nsamp   = length(1:sampleThin:nsampAll);
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
    #nsamp    += length(1:sampleThin:nsampTmp);
    if any(x->x=="restartfile", names(f))
      file = read(f,"restartfile");
      files = [ file; files ];
    else
      origfile = true;
    end
    close(f);
  end
  
  #handle cases where thinning parameters are ratios rather than a fixed number
  ( sampleThin < 1 ) && ( sampleThin = round(Int,sampleThin*nsampAll) );
  ( obsThin    < 1 ) && ( obsThin    = round(Int,obsThin   *nsampAll) );
  ( lpdfThin   < 1 ) && ( lpdfThin   = round(Int,lpdfThin  *nsampAll) );
  ( arThin     < 1 ) && ( arThin     = round(Int,arThin    *nsampAll) );

  println("sampleThin: $(sampleThin)");
  println("   obsThin: $(obsThin)");
  println("  lpdfThin: $(lpdfThin)");
  println("    arThin: $(arThin)");

  sampIdx  = 1:sampleThin:nsampAll;
  obsIdx   = 1:obsThin:nsampAll;
  lpdfIdx  = 1:lpdfThin:nsampAll;
  arIdx    = 1:arThin:nsampAll;
  
  nsamp = length(sampIdx);
  nobs  = length(obsIdx);
  nlpdf = length(lpdfIdx);
  nar   = length(arIdx);

  samples = zeros(nsamp,sampdim);
  obs     = zeros(nobs,obsdim);
  lpdfs   = zeros(nlpdf,length((1:3)[lpdfCols]));
  ar      = zeros(Int64,nar);
  sampCnt = obsCnt = lpdfCnt = arCnt = 0;
  sampCntAll = obsCntAll = lpdfCntAll = arCntAll = 0;
  for file in files
    f = h5open(file,"r");
    #read samples
    n  = size(f["samples"],1);
    # #nsampFl = length(1:sampleThin:n);
    # #samples[sampCnt+1:sampCnt+nsampFl,:] = f["samples"][1:sampleThin:end,sampleCols];
    # sampIdxFl = sampIdx .- sampCnt;
    # sampIdxFl = sampIdxFl[ (sampIdxFl .> 0) .& (sampIdxFl .<= n) ];
    # sampIdxFl = sampIdx[ (sampIdx .> sampCnt) .& (sampIdx .<= n+sampCnt) ] .- sampCnt;
    #sampIdxFl = (sampIdx .- sampCnt)[ (sampIdx .- sampCnt) .> 0 ][1]:sampleThin:n;
    sampIdxFl = sampIdx[sampCnt+1]-sampCntAll:sampleThin:n;
    nsampFl = length(sampIdxFl);
    samples[sampCnt+1:sampCnt+nsampFl,:] = f["samples"][sampIdxFl,sampleCols];
    sampCnt += nsampFl;
    sampCntAll += n;
    #read obs
    # #nobsFl = length(1:obsThin:n);
    # #obs[obsCnt+1:obsCnt+nobsFl,:] = f["obs"][1:obsThin:end,obsCols];
    # obsIdxFl = obsIdx .- obsCnt;
    # obsIdxFl = obsIdxFl[ (obsIdxFl .> 0) .& (obsIdxFl .<= n) ];
    # obsIdxFl = obsIdx[ (obsIdx .> obsCnt) .& (obsIdx .<= n+obsCnt) ] .- obsCnt;
    # obsIdxFl = (obsIdx .- obsCnt)[ (obsIdx .- obsCnt) .> 0 ][1]:obsThin:n;
    obsIdxFl = obsIdx[obsCnt+1]-obsCntAll:obsThin:n;
    nobsFl = length(obsIdxFl);
    obs[obsCnt+1:obsCnt+nobsFl,:] = f["obs"][obsIdxFl,obsCols];
    obsCnt += nobsFl;
    obsCntAll += n;
    #read lpdf
    # #nlpdfFl = length(1:lpdfThin:n);
    # #lpdfs[lpdfCnt+1:lpdfCnt+nlpdfFl,:] = f["lpdfs"][1:lpdfThin:end,lpdfCols];
    # lpdfIdxFl = lpdfIdx .- lpdfCnt;
    # lpdfIdxFl = lpdfIdxFl[ (lpdfIdxFl .> 0) .& (lpdfIdxFl .<= n) ];
    # lpdfIdxFl = lpdfIdx[ (lpdfIdx .> lpdfCnt) .& (lpdfIdx .<= n+lpdfCnt) ] .- lpdfCnt;
    #lpdfIdxFl = (lpdfIdx .- lpdfCnt)[ (lpdfIdx .- lpdfCnt) .> 0 ][1]:lpdfThin:n;
    lpdfIdxFl = lpdfIdx[lpdfCnt+1]-lpdfCntAll:lpdfThin:n;
    nlpdfFl = length(lpdfIdxFl);
    lpdfs[lpdfCnt+1:lpdfCnt+nlpdfFl,:] = f["lpdfs"][lpdfIdxFl,lpdfCols];
    lpdfCnt += nlpdfFl;
    lpdfCntAll += n;
    #read ar 
    # #narFl = length(1:arThin:n);
    # #ar[arCnt+1:arCnt+narFl,:] = f["ar"][1:arThin:end];
    # arIdxFl = arIdx .- arCnt;
    # arIdxFl = arIdxFl[ (arIdxFl .> 0) .& (arIdxFl .<= n) ];
    # arIdxFl = arIdx[ (arIdx .> arCnt) .& (arIdx .<= n+arCnt) ] .- arCnt;
    #arIdxFl = (arIdx .- arCnt)[ (arIdx .- arCnt) .> 0 ][1]:arThin:n;
    arIdxFl = arIdx[arCnt+1]-arCntAll:arThin:n;
    narFl = length(arIdxFl);
    ar[arCnt+1:arCnt+narFl,:] = f["ar"][arIdxFl];
    arCnt += narFl;
    arCntAll += n;
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

