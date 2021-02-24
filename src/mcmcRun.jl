function mcmcRun(mcmcP::mcmcProb, cur::mcmcSample; verbose=0, outFile="none", ncheck=10, samplesFlush=1)

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
  for i=1:mcmcP.nburn
    (verbose>0) && @printf("\nBurn in Step #%d:\n",i);
    cur,_,_ = mcmcP.step(cur, mcmcP; verbose=verbose);
  end

  ##Write data to HDF5 archive (unless user specified "none")
  #if checkpoint
  #  prot = isfile(outFile) ? "r+" : "w";
  #  f = h5open(outFile, prot);  #Open HDF5 file for writing

  #  f_samples = d_create(f, "samples", datatype(eltype(samples)), (size(samples),(mcmcP.nsamp,size(samples,2))), "chunk", (nSampMem,size(samples,2)));
  #  f_ar = d_create(f, "ar", datatype(eltype(ar)), ((nSampMem,),(mcmcP.nsamp,)), "chunk", (nSampMem,));
  #  f_obs = d_create(f, "obs", datatype(eltype(obs)), (size(obs),(mcmcP.nsamp,size(obs,2))), "chunk", (nSampMem,size(obs,2)));
  #  f_lpdfs = d_create(f, "lpdfs", datatype(eltype(lpdfs)), ((nSampMem,3),(mcmcP.nsamp,3)), "chunk", (nSampMem,3));
  #end

  #(verbose>0) && @printf "\n------------Starting Metropolis-Hastings Sampling------------\n"
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

    #checkpoint
    if checkpoint && ( i == nSampMem || ns == mcmcP.nsamp )
      (verbose>0) && @printf("\n---- Saving checkpoint at Sample %d (%d samples to save) ----\n",ns,i);

      #prot = isfile(outFile) ? "r+" : "w";
      #f = h5open(outFile, prot);  #Open HDF5 file for writing
      f = h5open(outFile, "cw");  #Open HDF5 file for writing

      ##create datasets on first checkpoint
      #if ( ns == nSampMem || ns == mcmcP.nsamp )
      if !exists(f, "samples")
        f_samples = d_create(f, "samples", datatype(eltype(samples)), (size(samples),(mcmcP.nsamp,size(samples,2))), "chunk", (nSampMem,size(samples,2)));
        f_ar      = d_create(f, "ar", datatype(eltype(ar)), ((nSampMem,),(mcmcP.nsamp,)), "chunk", (nSampMem,));
        f_obs     = d_create(f, "obs", datatype(eltype(obs)), (size(obs),(mcmcP.nsamp,size(obs,2))), "chunk", (nSampMem,size(obs,2)));
        f_lpdfs   = d_create(f, "lpdfs", datatype(eltype(lpdfs)), ((nSampMem,3),(mcmcP.nsamp,3)), "chunk", (nSampMem,3));
        write(f,"sampComplete",ns);
      
      #resize on subsequent checkpoints
      else
        f_samples = f["samples"];
        f_ar      = f["ar"];
        f_obs     = f["obs"];
        f_lpdfs   = f["lpdfs"];

        #resize
        set_dims!(f_samples,(ns,size(samples,2)));
        set_dims!(f_ar,(ns,));
        set_dims!(f_obs,(ns,size(obs,2)));
        set_dims!(f_lpdfs,(ns,3));
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

function mcmcRun(mcmcP::mcmcProb, s0=rand(mcmcP.prior); verbose=0, outFile="none", ncheck=10, samplesFlush=1)
  #initialize
  cur = mcmcSample();
  cur.samp = s0;

  return mcmcRun(mcmcP, cur; verbose=verbose, outFile=outFile, ncheck=ncheck, samplesFlush=samplesFlush);
end
