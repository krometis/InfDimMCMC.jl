function mcmcRun(mcmcP::mcmcProb, cur::mcmcSample; verbose=0, outFile="none", ncheck=10, samplesFlush=1, targetAR=0.0, sampAdapt=max(1,round(Int64,mcmcP.nsamp/100)),maxAdapt=((x)->(2-x)))

  #number of iterations to keep in memory between checkpoints
  checkpoint = (outFile != "none");
  nSampMem = checkpoint ? ceil(Int,mcmcP.nsamp/ncheck) : mcmcP.nsamp;

  #initialize
  mcmcFillSample(cur,mcmcP);
  
  #run mcmc forward
  #samples = Array{typeof(cur.samp)}(undef,mcmcP.nsamp);
  #ar      = Array{UInt8}(undef,mcmcP.nsamp);
  #obs     = Array{typeof(cur.obs)}(undef,mcmcP.nsamp);
  #lpdfs   = Array{Float64}(undef,mcmcP.nsamp,3);
  samples = Array{eltype(cur.samp)}(undef,nSampMem,length(cur.samp));
  ar      = Array{UInt8}(undef,nSampMem);
  obs     = Array{eltype(cur.obs)}(undef,nSampMem,length(cur.obs));
  lpdfs   = Array{Float64}(undef,nSampMem,3);

  #burnin
  (verbose>0) && @printf("\n-----------------------Starting Burn In-----------------------\n");
  accCnt  = 0;
  sampCnt = 0;
  for i=1:mcmcP.nburn
    (verbose>1) && @printf("\nBurn in Step #%d:\n",i);
    cur,accept,_ = mcmcP.step(cur, mcmcP; verbose=verbose);
    accCnt  += accept;
    sampCnt += 1;
    if (targetAR > 0.0) && (i % sampAdapt == 0)
      mcmcAdapt(mcmcP,accCnt/sampCnt,targetAR;verbose=1);
      accCnt  = 0;
      sampCnt = 0;
    end
  end

  (verbose>0) && @printf("\n----------------------Starting Sampling-----------------------\n");
  ns = 1; #total number of samples
  i  = 1; #samples since last checkpoint
  while ns <= mcmcP.nsamp
    (verbose>1) && @printf("\nSample #%d:\n",ns);
    cur,ar[i],can = mcmcP.step(cur, mcmcP; verbose=verbose);
    samples[i,:] = cur.samp;
    obs[i,:]     = cur.obs;
    lpdfs[i,1]   = cur.prLogPdf;
    lpdfs[i,2]   = cur.llLogPdf;
    lpdfs[i,3]   = cur.postLogPdf;

    #adapt
    if (targetAR > 0.0) && (ns % sampAdapt == 0)
      mcmcAdapt(mcmcP,mean(ar[i-sampAdapt+1:i]),targetAR;verbose=1,maxAdapt=maxAdapt(ns/mcmcP.nsamp));
    end

    #checkpoint
    if checkpoint && ( i == nSampMem || ns == mcmcP.nsamp )
      (verbose>0) && @printf("\n---- Saving checkpoint at Sample %d (%d samples to save) ----\n",ns,i);

      f = h5open(outFile, "cw");  #Open HDF5 file for writing

      ##create datasets on first checkpoint
      #if ( ns == nSampMem || ns == mcmcP.nsamp )
      if !haskey(f, "samples")
        f_samples = create_dataset(f, "samples", datatype(eltype(samples)), (size(samples),(mcmcP.nsamp,size(samples,2))), chunk=(nSampMem,size(samples,2)));
        f_ar      = create_dataset(f, "ar",      datatype(eltype(ar)),      ((nSampMem,),  (mcmcP.nsamp,)), chunk=(nSampMem,));
        f_obs     = create_dataset(f, "obs",     datatype(eltype(obs)),     (size(obs),    (mcmcP.nsamp,size(obs,2))), chunk=(nSampMem,size(obs,2)));
        f_lpdfs   = create_dataset(f, "lpdfs",   datatype(eltype(lpdfs)),   ((nSampMem,3), (mcmcP.nsamp,3)), chunk=(nSampMem,3));
        write(f,"sampComplete",ns);
      
      #resize on subsequent checkpoints
      else
        f_samples = f["samples"];
        f_ar      = f["ar"];
        f_obs     = f["obs"];
        f_lpdfs   = f["lpdfs"];

        #resize
        HDF5.set_extent_dims(f_samples,(ns,size(samples,2)));
        HDF5.set_extent_dims(f_ar,(ns,));
        HDF5.set_extent_dims(f_obs,(ns,size(obs,2)));
        HDF5.set_extent_dims(f_lpdfs,(ns,3));
      end

      #save values
      write(f["sampComplete"],ns);
      f_samples[ns-i+1:ns,:] = samples[1:i,:];
      f_ar[ns-i+1:ns]        = ar[1:i];
      f_obs[ns-i+1:ns,:]     = obs[1:i,:];
      f_lpdfs[ns-i+1:ns,:]   = lpdfs[1:i,:];

      close(f);

      i=0;
    end #end checkpoint

    i +=1;
    ns+=1;

  end #end iteration
  (verbose>0) && @printf("\n------------Sampling Complete------------\n");

  #flush stdout if ns is a multiple of samplesFlush
  ( ns % samplesFlush == 0 ) && flush(stdout);

  return samples, obs, lpdfs, ar;
end

#function mcmcRun(mcmcP::mcmcProb, s0=rand(mcmcP.prior); verbose=0, outFile="none", ncheck=10, samplesFlush=1, targetAR=false, sampAdapt=round(Int64,mcmcP.nsamp/100))
function mcmcRun(mcmcP::mcmcProb, s0=rand(mcmcP.prior); kwargs...)
  #initialize
  cur = mcmcSample();
  cur.samp = s0;

  return mcmcRun(mcmcP, cur; kwargs...);
end
