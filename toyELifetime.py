#!/usr/bin/env python

import ROOT as root
import numpy
import sys
from matplotlib import pylab as mpl

RAND = root.TRandom3(7)

def toyCluster(nPoints,qMPV,lifetimeTrue,trackSlope=0.14,usPerBin=100.,suffix="",doLogFit=False,doPlots=True):
  """
  qMPV is the true charge deposited
  lifetimeTrue is the true lifetime in us
  trackSlope is points / us
  """
  landauWidth = qMPV*0.22
  nBins = int(nPoints // (trackSlope*usPerBin))

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
      qTrue = RAND.Gaus(qMPV,landauWidth)
      ##qTrue = RAND.Landau(qMPV,landauWidth)
      ###qTrue = RAND.Landau(param1,param2)
      t = iPoint / trackSlope
      qMeas = qTrue*numpy.exp(-t/lifetimeTrue)
  
      ts[iPoint] = t
      qTrues[iPoint] = qTrue
      qMeass[iPoint] = qMeas
  
      ## Now do the fit stuff
      iBin = int(t // usPerBin)
      if doLogFit:
        if iBin < nBins:
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
  
  maxChg = ave * 5
  minChg = ave * 0.2
  
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
  
  if doPlots:
    ax.set_ylim(0,qMPV*8)
    ax.scatter(ts,qMeass,2.,c='k',lw=0)
    ax.plot(ts,qMPV*numpy.exp(-ts/lifetimeTrue),'-g')
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
    ax.plot(ts,numpy.log(qMPV)-ts/lifetimeTrue,"m-")
    ax.set_xlabel("Drift Time [us]")
    ax.set_ylabel("log(Charge)")
    ax.set_ylim(4.5,7)
    fig.savefig("LandauLog{}.png".format(suffix))
    fig.savefig("LandauLog{}.pdf".format(suffix))

  return 1./lifeInv

if __name__ == "__main__":
  nPoints = 200
  #nPoints = 300
  #nPoints = 400
  trackSlope = 0.2 # points / us 
  #trackSlope = 0.3 # points / us 
  #trackSlope = 0.5 # points / us 
  qMPV = 300.
  lifetimeTrue = 3000. # us
  doLogFit = False

  landauPoints = numpy.array([RAND.Landau(qMPV,qMPV*0.22) for i in range(100000)])
  fig, ax = mpl.subplots()
  ax.hist(landauPoints,bins=100,range=[0,2000],histtype='step')
  ax.axvline(numpy.mean(landauPoints))
  ax.set_xlabel("Landau Distributed Points")
  ax.set_ylabel("Points / Bin")
  fig.savefig("ToyLandau.png")
  fig.savefig("ToyLandau.pdf")

  logLandauPoints = numpy.log(landauPoints)
  logmean = numpy.mean(logLandauPoints)
  logrms = numpy.std(logLandauPoints)
  truncatedPoints = logLandauPoints[logLandauPoints < logmean+1.0*logrms]
  logtruncmean = numpy.mean(truncatedPoints)
  logtruncrms = numpy.std(truncatedPoints)
  logtruncmeanerror = logtruncrms*len(truncatedPoints)**(-0.5)

  fig, ax = mpl.subplots()
  ax.axvspan(logmean-logrms,logmean+logrms,fc="r",alpha=0.3)
  ax.axvline(logmean,c='r')
  #ax.axvspan(logtruncmean-logtruncrms,logtruncmean+logtruncrms,fc="b",alpha=0.3)
  ax.axvspan(logtruncmean-logtruncmeanerror,logtruncmean+logtruncmeanerror,fc="b",alpha=0.3)
  ax.axvline(logtruncmean,c='b')
  ax.hist(logLandauPoints,bins=100,range=[4,12],histtype='step',color='k')
  ax.set_xlabel("Log-Landau Distributed Points")
  ax.set_ylabel("Points / Bin")
  fig.savefig("ToyLogLandau.png")
  fig.savefig("ToyLogLandau.pdf")
  

  lifes = []
  for iCluster in range(1000):
    doPlots = (iCluster < 5)
    life = toyCluster(nPoints,qMPV,lifetimeTrue,trackSlope=trackSlope,suffix="_{}".format(iCluster),doLogFit=doLogFit,doPlots=doPlots)
    lifes.append(life/1000.)

  fig, ax = mpl.subplots()
  ax.hist(lifes,bins=30,range=[0,5],histtype='step')
  ax.axvline(lifetimeTrue/1000.,c='g')
  ax.set_xlabel("Electron Lifetime [ms]")
  ax.set_ylabel("Toy Clusters / Bin")
  fig.text(0.15,0.9,"Hits: {}, Slope: {} us".format(nPoints,1./trackSlope),ha='left')
  fig.savefig("ToyLifetime_{}_{}.png".format(nPoints,1./trackSlope))
  fig.savefig("ToyLifetime_{}_{}.pdf".format(nPoints,1./trackSlope))
