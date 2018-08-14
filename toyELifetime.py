#!/usr/bin/env python

import ROOT as root
import numpy
import sys
from matplotlib import pylab as mpl

RAND = root.TRandom3(7)

def toyCluster(qMPV,lifetimeTrue,nBins=10,pointsPerBin=20,usPerBin=100.,suffix="",doLogFit=False,doPlots=True,doGaus=False):
  """
  qMPV is the true charge deposited
  lifetimeTrue is the true lifetime in us
  trackSlope is points / us
  """
  landauWidth = qMPV*0.22
  #nBins = int(nPoints // (trackSlope*usPerBin))
  nBins = int(nBins)
  nPoints = int(pointsPerBin * nBins)
  trackSlope = pointsPerBin / usPerBin

  ts = numpy.zeros(nPoints) # in us
  qTrues = numpy.zeros(nPoints) # in ADC
  qMeass = numpy.zeros(nPoints) # in ADC
  
  tck = numpy.zeros(nBins)
  ave = numpy.zeros(nBins)
  err = numpy.zeros(nBins)
  cnt = numpy.zeros(nBins)
  
  fig = None
  ax = None
  if doPlots:
    fig, ax = mpl.subplots()
  
  for iPoint in range(nPoints):
      # correct for dumb root MPV!!!!
      param2 = landauWidth
      param1 = qMPV*0.22278+param2
      qTrue = None
      if doGaus:
        qTrue = RAND.Gaus(qMPV,landauWidth)
      else:
        qTrue = RAND.Landau(qMPV,landauWidth)
        #qTrue = RAND.Landau(param1,param2)
      t = iPoint / trackSlope
      qMeas = qTrue*numpy.exp(-t/lifetimeTrue)
  
      ts[iPoint] = t
      qTrues[iPoint] = qTrue
      qMeass[iPoint] = qMeas
  
      ## Now do the fit stuff
      iBin = int(t // usPerBin)
      if doLogFit:
        if iBin < nBins and qMeas < 1500.:
          tck[iBin] += t
          ave[iBin] += numpy.log(qMeas)
          err[iBin] += numpy.log(qMeas)**2
          cnt[iBin] += 1
      else:
        if iBin < nBins and qMeas < 1500.:
          tck[iBin] += t
          ave[iBin] += qMeas
          err[iBin] += qMeas*qMeas
          cnt[iBin] += 1

  tck /= cnt
  ave /= cnt
  
  if doPlots:
    ax.plot(tck,ave,'or')
  
  maxChg = ave * 1.3
  minChg = ave * 0.5
  
  if not doLogFit:
    tck = numpy.zeros(nBins)
    ave = numpy.zeros(nBins)
    err = numpy.zeros(nBins)
    cnt = numpy.zeros(nBins)
    
    # now truncate
    for iPoint in range(nPoints):
        t = ts[iPoint]
        qMeas = qMeass[iPoint]
        iBin = int(t // usPerBin)
        if iBin < nBins and qMeas > minChg[iBin] and qMeas < maxChg[iBin]:
          tck[iBin] += t
          ave[iBin] += qMeas
          err[iBin] += qMeas*qMeas
          cnt[iBin] += 1
    
    tck /= cnt
    ave /= cnt
  
  arg = err - cnt * ave * ave
  err = numpy.sqrt(arg / (cnt - 1))
  err /= numpy.sqrt(cnt)
  
  #for iBin in range(nBins):
  #    print iBin, tck[iBin],ave[iBin],err[iBin],cnt[iBin]
  
  
  # fit to find the 1/lifetime
  Sum = 0.;
  sumx = 0.;
  sumy = 0.;
  sumxy = 0.;
  sumx2 = 0.;
  sumy2 = 0.;
  xx = 0.
  yy = 0. 
  wght = 0. 
  arg = 0.
  fitcnt = 0;
  
  logPlot_xs = []
  logPlot_ys = []
    
  for ihist in range(nBins):
    if(cnt[ihist] < 3):
       continue;
    if(err[ihist] == 0):
      continue;
    xx = (tck[ihist] - tck[0]);
    if doLogFit:
      yy = ave[ihist]
    else:
      yy = numpy.log(ave[ihist]);
    logPlot_xs.append(xx)
    logPlot_ys.append(yy)
    # error on log(x) = dx / x
    arg = ave[ihist] / err[ihist];
    wght = arg * arg;
    Sum += wght;
    sumx += wght * xx;
    sumy += wght * yy;
    sumx2 += wght * xx * xx;
    sumxy += wght * xx * yy;
    sumy2 += wght * yy * yy;
    ++fitcnt;
  # calculate coefficients
  delta = Sum * sumx2 - sumx * sumx;
  #print "delta ",delta
  A = (sumx2 * sumy - sumx * sumxy) / delta;
  #print "A ", A
  # slope = 1 / lifetime
  B = (sumxy * Sum  - sumx * sumy) / delta;
  # calculate the error
  ndof = fitcnt - 2;
  varnce = (sumy2 + A*A*Sum + B*B*sumx2 - 2 * (A*sumy + B*sumxy - A*B*sumx)) / ndof;
  BErr = numpy.sqrt(varnce * Sum / delta);
  #print "B ",B," ",BErr;
  
  lifeInv = -B;
  lifeInvErr = BErr / (B * B) ;
  
  logFun_xs = []
  logFun_ys = []
  # calculate chisq
  chi2 = 0;
  for ihist in range(nBins):
    if(cnt[ihist] < 3):
       continue;
    if(err[ihist] == 0):
       continue;
    xx = (tck[ihist] - tck[0]);
    yy = numpy.exp(A - xx * lifeInv);
    arg = (yy - ave[ihist]) / err[ihist];
    chi2 += arg * arg;
    #print "chk ", ihist, " xx ", xx, " yy ", yy, " ave ", ave[ihist], " arg ", arg, "\n";
    logFun_xs.append(xx)
    logFun_ys.append(A-xx*lifeInv)
  chi2ndof = chi2 / ndof;

  ### Redo numpy

  redoave = numpy.zeros(nBins)
  redoavevariance = numpy.zeros(nBins)
  redoaveerr = numpy.zeros(nBins)
  redologave = numpy.zeros(nBins)
  redologaveerr = numpy.zeros(nBins)
  redot = numpy.zeros(nBins)
  qMeassLog = numpy.log(qMeass)
  for ihist in range(nBins):
    #ts, qMeass
    inBin = numpy.logical_and(ts > (ihist*usPerBin),ts < ((ihist+1)*usPerBin))
    nPoints = len(qMeass[inBin])
    redoave[ihist] = numpy.mean(qMeass[inBin])
    redoavevariance[ihist] = numpy.var(qMeass[inBin])/nPoints
    redot[ihist] = numpy.mean(ts[inBin])
    redologave[ihist] = numpy.mean(qMeassLog[inBin])
    redologaveerr[ihist] = numpy.std(qMeassLog[inBin])/numpy.sqrt(nPoints)
  redoaveerr = numpy.sqrt(redoavevariance)
  redoavelog = numpy.log(redoave)
  redoaveerrlog = numpy.abs(redoaveerr/redoave)
  
  coefs, cov = numpy.polyfit(redot,redoavelog,1,w=1./redoaveerrlog,cov=True)
  coefsLogFirst, covLogFirst = numpy.polyfit(redot,redologave,1,w=1./redologaveerr,cov=True)
  numpyLife = -1/coefs[0]
  numpyLogLife = -1/coefsLogFirst[0]
  numpyLifeVar = numpy.abs(cov[0][0]/coefs[0]**4)
  numpyLogLifeVar = numpy.abs(covLogFirst[0][0]/coefsLogFirst[0]**4)
  #print "Bruce: {:.2f}, Numpy: {:.2f} Numpy Log First: {:.2f}".format(1./lifeInv/1000.,-1/coefs[0]/1000.,-1/coefsLogFirst[0]/1000.)
  #print "Pull: Numpy: {:.2f} Numpy Log First: {:.2f}".format((numpyLife-lifetimeTrue)/numpy.sqrt(numpyLifeVar),(numpyLogLife-lifetimeTrue)/numpy.sqrt(numpyLogLifeVar))
  
  if doPlots:
    ax.set_ylim(0,qMPV*8)
    ax.scatter(ts,qMeass,2.,c='k',lw=0)
    ax.plot(ts,qMPV*numpy.exp(-ts/lifetimeTrue),'-m')
    ax.plot(ts,numpy.exp(ts*B+A),'-g')
    ax.plot(ts,numpy.exp(ts*coefs[0]+coefs[1]),'-c')
    #ax.plot(tck,ave,'og')
    ax.errorbar(numpy.arange(nBins)*usPerBin+0.5*usPerBin,ave,xerr=0.5*usPerBin,yerr=err,fmt="ob")
    ax.set_xlabel("Drift Time [us]")
    ax.set_ylabel("Charge")
    ax.set_ylim(0,1000)
    fig.savefig("Landau{}.png".format(suffix))
    fig.savefig("Landau{}.pdf".format(suffix))
    
    
    fig, ax = mpl.subplots()
    ax.scatter(ts-tck[0],numpy.log(qMeass),2.,c='k',lw=0)
    ax.plot(logPlot_xs,logPlot_ys,"bo")
    ax.plot(logFun_xs,logFun_ys,"g-")
    ax.plot(ts,ts*coefs[0]+coefs[1],'-c')
    ax.plot(ts,numpy.log(qMPV)-ts/lifetimeTrue,"m-")
    ax.set_xlabel("Drift Time [us]")
    ax.set_ylabel("log(Charge)")
    ax.set_ylim(4.5,7)
    fig.savefig("LandauLog{}.png".format(suffix))
    fig.savefig("LandauLog{}.pdf".format(suffix))

  return 1./lifeInv, numpyLife, numpyLogLife, numpyLifeVar, numpyLogLifeVar

if __name__ == "__main__":
  nBins = 10
  pointsPerBin = 100./nBins
  usPerBin = 100.
  qMPV = 300.
  lifetimeTrue = 3000. # us
  doLogFit = False
  doGaus = True

  #landauPoints = numpy.array([RAND.Landau(qMPV,qMPV*0.22) for i in range(100000)])
  #fig, ax = mpl.subplots()
  #ax.hist(landauPoints,bins=100,range=[0,2000],histtype='step')
  #ax.axvline(numpy.mean(landauPoints))
  #ax.set_xlabel("Landau Distributed Points")
  #ax.set_ylabel("Points / Bin")
  #fig.savefig("ToyLandau.png")
  #fig.savefig("ToyLandau.pdf")

  #logLandauPoints = numpy.log(landauPoints)
  #logmean = numpy.mean(logLandauPoints)
  #logrms = numpy.std(logLandauPoints)
  #truncatedPoints = logLandauPoints[logLandauPoints < logmean+1.0*logrms]
  #logtruncmean = numpy.mean(truncatedPoints)
  #logtruncrms = numpy.std(truncatedPoints)
  #logtruncmeanerror = logtruncrms*len(truncatedPoints)**(-0.5)

  #fig, ax = mpl.subplots()
  #ax.axvspan(logmean-logrms,logmean+logrms,fc="r",alpha=0.3)
  #ax.axvline(logmean,c='r')
  ##ax.axvspan(logtruncmean-logtruncrms,logtruncmean+logtruncrms,fc="b",alpha=0.3)
  #ax.axvspan(logtruncmean-logtruncmeanerror,logtruncmean+logtruncmeanerror,fc="b",alpha=0.3)
  #ax.axvline(logtruncmean,c='b')
  #ax.hist(logLandauPoints,bins=100,range=[4,12],histtype='step',color='k')
  #ax.set_xlabel("Log-Landau Distributed Points")
  #ax.set_ylabel("Points / Bin")
  #fig.savefig("ToyLogLandau.png")
  #fig.savefig("ToyLogLandau.pdf")
  

  lifes = []
  lifesNumpy = []
  lifesLogNumpy = []
  pullsNumpy = []
  pullsLogNumpy = []
  for iCluster in range(1000):
    doPlots = (iCluster < 5)
    life, lifeNumpy, lifeLogNumpy, lifeNumpyVar, lifeLogNumpyVar = toyCluster(qMPV,lifetimeTrue,nBins,pointsPerBin,usPerBin,suffix="_{}".format(iCluster),doLogFit=doLogFit,doPlots=doPlots,doGaus=doGaus)
    lifes.append(life/1000.)
    lifesNumpy.append(lifeNumpy/1000.)
    lifesLogNumpy.append(lifeLogNumpy/1000.)
    pullsNumpy.append(lifeNumpy/numpy.sqrt(lifeNumpyVar))
    pullsLogNumpy.append(lifeLogNumpy/numpy.sqrt(lifeLogNumpyVar))

  fig, ax = mpl.subplots()
  ax.hist(lifes,bins=30,range=[0,6],histtype='step')
  ax.hist(lifesNumpy,bins=30,range=[0,6],histtype='step')
  ax.hist(lifesLogNumpy,bins=30,range=[0,6],histtype='step')
  ax.axvline(lifetimeTrue/1000.,c='g')
  ax.set_xlabel("Electron Lifetime [ms]")
  ax.set_ylabel("Toy Clusters / Bin")
  distType = "Landau"
  distTypeLabel = ""
  if doGaus:
    distType = "Gaus"
    distTypeLabel = "Gaussian Charge "
  fig.text(0.15,0.9,"{}Hits: {} Bins: {}, Hits/Bin: {:.1f}".format(distTypeLabel,int(nBins*pointsPerBin),nBins,pointsPerBin),ha='left')
  fig.savefig("ToyLifetime_{}_bins{}_hitpbin{:.0f}.png".format(distType,nBins,pointsPerBin))
  fig.savefig("ToyLifetime_{}_bins{}_hitpbin{:.0f}.pdf".format(distType,nBins,pointsPerBin))

  fig, ax = mpl.subplots()
  ax.hist(pullsNumpy,bins=30,range=[-10,10],histtype='step')
  ax.hist(pullsLogNumpy,bins=30,range=[-10,10],histtype='step')
  ax.set_xlabel("Electron Lifetime Pull")
  ax.set_ylabel("Toy Clusters / Bin")
  fig.text(0.15,0.9,"{}Hits: {} Bins: {}, Hits/Bin: {:.1f}".format(distTypeLabel,int(nBins*pointsPerBin),nBins,pointsPerBin),ha='left')
  fig.savefig("Pulls.png")
