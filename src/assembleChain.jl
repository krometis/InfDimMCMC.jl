# Assemble a chain split across multiple files
function assembleChain(file; verbose=true, sampleCols=:, obsCols=:, lpdfCols=:, sampleThin=1, obsThin=sampleThin, lpdfThin=sampleThin, arThin=sampleThin)
  files = [ file ];
  origfile = false;
  
  #get information from first file
  f = h5open(file,"r");
  nsampAll, sampdim = size(f["samples"]);

  sampdim = length((1:sampdim)[sampleCols]);
  obsdim  = size(f["obs"],2);
  obsdim  = length((1:obsdim)[obsCols]);
  if haskey(f,"restartfile")
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
    if haskey(f,"restartfile")
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

  #println("sampleThin: $(sampleThin)");
  #println("   obsThin: $(obsThin)");
  #println("  lpdfThin: $(lpdfThin)");
  #println("    arThin: $(arThin)");

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
    sampIdxFl = sampIdx[sampCnt+1]-sampCntAll:sampleThin:n;
    nsampFl = length(sampIdxFl);
    samples[sampCnt+1:sampCnt+nsampFl,:] = f["samples"][sampIdxFl,sampleCols];
    sampCnt += nsampFl;
    sampCntAll += n;

    #read obs
    obsIdxFl = obsIdx[obsCnt+1]-obsCntAll:obsThin:n;
    nobsFl = length(obsIdxFl);
    obs[obsCnt+1:obsCnt+nobsFl,:] = f["obs"][obsIdxFl,obsCols];
    obsCnt += nobsFl;
    obsCntAll += n;

    #read lpdf
    lpdfIdxFl = lpdfIdx[lpdfCnt+1]-lpdfCntAll:lpdfThin:n;
    nlpdfFl = length(lpdfIdxFl);
    lpdfs[lpdfCnt+1:lpdfCnt+nlpdfFl,:] = f["lpdfs"][lpdfIdxFl,lpdfCols];
    lpdfCnt += nlpdfFl;
    lpdfCntAll += n;

    #read ar 
    arIdxFl = arIdx[arCnt+1]-arCntAll:arThin:n;
    narFl = length(arIdxFl);
    ar[arCnt+1:arCnt+narFl,:] = f["ar"][arIdxFl];
    arCnt += narFl;
    arCntAll += n;

    close(f);
    verbose && @printf("Read %d samples from %s\n",n,file);
  end
  verbose && @printf("Done. %d total samples read out of %d found.\n", nsamp, nsampAll);
  
  return samples, obs, lpdfs, ar, files;
end
