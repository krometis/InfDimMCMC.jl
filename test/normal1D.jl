using Distributions
using InfDimMCMC
using Printf
using Statistics
using LinearAlgebra
using HDF5

meanPrior=[ 0.0];
meanLlh  =[ 2.0];
nsamp=1e5;
meanThresh = 0.05;
stdThresh  = 0.05;
outFile="normal1D_test.h5";

function mcmcTest(cur,m::mcmcProb,smplr::String; meanLlh=meanLlh, meanThresh=meanThresh, stdThresh=stdThresh,outFile=outFile)
  
  run(`rm -f $(outFile)`);

  @printf("\nSampler: %s\n", smplr);
  #@printf("Initial sample: %6.4f\n", cur.samp);

  mcmcSetSampler(m,smplr);
  #m.step = mcmcStep;
  samples, obs, lpdfs, ar = mcmcRun(m,s0;outFile=outFile,verbose=0);

  @printf("Acceptance ratio: %6.4f\n",mean(ar));
  
  meanTrue = mean([mean(m.prior),meanLlh]);
  meanComp = mean(samples, dims=1)[:];
  meanErr  = norm(meanTrue - meanComp);
  stdTrue  = 1.0/sqrt(2);
  stdComp  = std(samples);
  stdErr   = abs(stdTrue - stdComp);
  @printf("mean: %6.4f %6.4f %6.4f\n",meanComp[1],meanTrue[1],meanErr);
  @printf(" std: %6.4f %6.4f %6.4f\n", stdComp, stdTrue, stdErr);

  @test meanErr < meanThresh;
  @test  stdErr <  stdThresh;

  #check that samples were written correctly
  nSampMem = length(ar);
  @test samples == h5read(outFile,"samples")[end-nSampMem+1:end,:];
  @test obs     == h5read(outFile,"obs"    )[end-nSampMem+1:end,:];
  @test lpdfs   == h5read(outFile,"lpdfs"  )[end-nSampMem+1:end,:];
  @test ar      == h5read(outFile,"ar"     )[end-nSampMem+1:end];

  #return samples, samplesFile
end

#define MCMC problem
m = mcmcProb();
m.nsamp=nsamp;                     #samples
m.nburn=0;                         #burnin
m.prior = MvNormal(meanPrior,1.0); #prior

#define potential (negative log likelihood) and its gradient
InfDimMCMC.mcmcPotMap(s) = -logpdf( MvNormal(meanLlh,1.0), s.samp );
InfDimMCMC.mcmcGradPotMap(s) = (s.samp - meanLlh);


#initial sample
s0 = rand(m.prior);

@printf("\n ----- TEST: 1D NORMAL ----- \n");

#test various MCMC methods
mcmcTest(s0,m,"is");
mcmcTest(s0,m,"pcn|0.5");
mcmcTest(s0,m,"mala|0.5");
mcmcTest(s0,m,"hmc|0.5|2");
mcmcTest(s0,m,"mppcn|0.5|32");
#mcmcTest(s0,m,"bbpcn|0.5|32"); #this will produce the wrong answer, failing the test

@printf("\n -----  END: 1D NORMAL ----- \n");
