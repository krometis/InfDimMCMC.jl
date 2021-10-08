#mcmcAdapt: adjust MCMC parameters (e.g., step size) to target a particular accept/reject ratio
function mcmcAdapt(mcmcP::mcmcProb, acceptRatio::Float64, targetAR::Float64; verbose=0, minAdapt=0.5, maxAdapt=2.0)

  adaptRatio = targetAR / acceptRatio;
  adaptRatio = min(adaptRatio,maxAdapt);
  adaptRatio = max(adaptRatio,minAdapt);

  d = mcmcP.mcmc;
  if d["mcmc"] == :pcn
    d["beta"] *= adaptRatio;
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

  else
    error("Adaptive MCMC not supported for selected method.");
  end
end
