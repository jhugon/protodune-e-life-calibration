#!/usr/bin/env python

from helpers import *
#import ROOT as root
import numpy
from matplotlib import pylab as mpl

from toyELifetime import generateCluster, ChargeRatioMethod

if __name__ == "__main__":
  print "Start time: ", datetime.datetime.now().replace(microsecond=0).isoformat(' ')

  nClusters = 10000
  clustersPerCalc = 1000
  qMPV = 300.
  lifetimeTrue = 3000. # us
  doGaus = False
  hitsPerusCases = [40./100.,20./100.]
  #nHitsCases = [20,50,100,150,200,300,400]
  nHitsCases = [100,200,400]
  crms = []
  for i, hitsPerus in enumerate(hitsPerusCases):
    distType = "Landau"
    distTypeLabel = "Landau Charge "
    if doGaus:
      distType = "Gaus"
      distTypeLabel = "Gaussian Charge "
    crms.append([])
    for j, hits in enumerate(nHitsCases):
      caseStr = "_{}_hits{}_hitsPerus{:.1f}".format(distType,hits,hitsPerus).replace('.','p')
      crms[i].append(ChargeRatioMethod(caseStr))

  lifetimesShape = (len(hitsPerusCases),len(nHitsCases),nClusters // clustersPerCalc)
  lifetimes = numpy.zeros(lifetimesShape)
  lifetimeErrs = numpy.zeros(lifetimesShape)
  chi2s = numpy.zeros(lifetimesShape)

  for iCluster in range(nClusters):
    for i, hitsPerus in enumerate(hitsPerusCases):
      ts, qMeass = generateCluster(qMPV,lifetimeTrue,nHitsCases[-1],hitsPerus,doGaus=doGaus)
      for j, hits in enumerate(nHitsCases):
        crms[i][j].processCluster(ts[:hits],qMeass[:hits])
        if (iCluster+1) % clustersPerCalc == 0:
          k = iCluster // clustersPerCalc
          lifetimes[i][j][k], lifetimeErrs[i][j][k], chi2s[i][j][k] = crms[i][j].calculate("_clusters{}".format(iCluster+1))

  pulls = (lifetimes - lifetimeTrue)/(lifetimeErrs)
  
  fig, ax = mpl.subplots()
  nClustersCalcd = numpy.arange(nClusters // clustersPerCalc)*clustersPerCalc
  for i, hitsPerus in enumerate(hitsPerusCases):
    for j, hits in enumerate(nHitsCases):
      ax.plot(nClustersCalcd,pulls[i,j,:],label="{} Hits {:.1f} Hits/us".format(hits,hitsPerus))
  ax.legend()
  ax.set_xlabel("N Clusters")
  ax.set_ylabel(r"$(\tau_{meas}-\tau_{true})/\sigma$")
  fig.text(0.15,0.955,"Toy Study")
  fig.savefig("EvalChargeratioMethod_Pulls.png")
  fig.savefig("EvalChargeratioMethod_Pulls.pdf")
  mpl.close(fig)

  fig, ax = mpl.subplots()
  for i, hitsPerus in enumerate(hitsPerusCases):
    for j, hits in enumerate(nHitsCases):
      ax.plot(nClustersCalcd,lifetimes[i,j,:],label="{} Hits {:.1f} Hits/us".format(hits,hitsPerus))
  ax.legend()
  ax.set_xlabel("N Clusters")
  ax.set_ylabel(r"Electron Lifetime Estimate [us]")
  fig.text(0.15,0.955,"Toy Study")
  fig.savefig("EvalChargeratioMethod_Lifes.png")
  fig.savefig("EvalChargeratioMethod_Lifes.pdf")
  mpl.close(fig)

  print "End time: ", datetime.datetime.now().replace(microsecond=0).isoformat(' ')
