#mcmcAdapt: adjust MCMC parameters (e.g., step size) to target a particular accept/reject ratio
function mcmcAdapt(mcmcP::mcmcProb, acceptRatio::Float64, targetAR::Float64; verbose=0, maxAdapt=2.0, minAdapt=1.0/maxAdapt)

  adaptRatio = acceptRatio / targetAR;
  adaptRatio = min(adaptRatio,maxAdapt);
  adaptRatio = max(adaptRatio,minAdapt);

  d = mcmcP.mcmc;
  if d["mcmc"] == :pcn
    d["beta"] *= adaptRatio;
    d["beta"] = min(1.0,d["beta"]); #beta>1 doesn't make sense
    mcmcSetSampler(mcmcP,d; verbose=0);
    (verbose > 0) && @printf("mcmcAdapt: acceptRatio=%5.3f, target=%5.3f. Adjusting step size beta by factor of %5.3f to %8.6f\n",acceptRatio, targetAR, adaptRatio, d["beta"]);

  elseif d["mcmc"] == :hmc
    d["eps"] *= adaptRatio;
    d["tau"] *= adaptRatio;
    mcmcSetSampler(mcmcP,d; verbose=0);
    (verbose > 0) && @printf("mcmcAdapt: acceptRatio=%5.3f, target=%5.3f. Adjusting step size eps and integration time tau by factor of %5.3f to %8.6f and %8.6f, respectively.\n",acceptRatio, targetAR, adaptRatio, d["eps"], d["tau"]);

  elseif d["mcmc"] == :mala
    d["h"] *= adaptRatio;
    mcmcSetSampler(mcmcP,d; verbose=0);
    (verbose > 0) && @printf("mcmcAdapt: acceptRatio=%5.3f, target=%5.3f. Adjusting step size h by factor of %5.3f to %8.6f\n",acceptRatio, targetAR, adaptRatio, d["h"]);

  elseif d["mcmc"] == :mppcn || d["mcmc"] == :bbpcn
    d["beta"] *= adaptRatio;
    d["beta"] = min(1.0,d["beta"]); #beta>1 doesn't make sense
    mcmcSetSampler(mcmcP,d; verbose=0);
    (verbose > 0) && @printf("mcmcAdapt: acceptRatio=%5.3f, target=%5.3f. Adjusting step size beta by factor of %5.3f to %8.6f\n",acceptRatio, targetAR, adaptRatio, d["beta"]);

  else
    error("Adaptive MCMC not supported for selected method.");
  end
end
