#!/usr/bin/env python

from helpers import *
#import ROOT as root
import numpy
import sys
from matplotlib import pylab as mpl

root.gROOT.SetBatch(True)
RAND = root.TRandom3(7)

def toyCluster(qMPV,lifetimeTrue,nBins=10,pointsPerBin=20,usPerBin=100.,suffix="",doLogFit=False,doPlots=True,doGaus=False,chargeRatioVdtHist=None):
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
  numpyLife = 0.
  numpyLogLife = 0.
  numpyLifeVar = 0.
  numpyLogLifeVar = 0.

  if False:
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
      nGoodPoints = len(qMeass[inBin])
      redoave[ihist] = numpy.mean(qMeass[inBin])
      redoavevariance[ihist] = numpy.var(qMeass[inBin])/nGoodPoints
      redot[ihist] = numpy.mean(ts[inBin])
      redologave[ihist] = numpy.mean(qMeassLog[inBin])
      redologaveerr[ihist] = numpy.std(qMeassLog[inBin])/numpy.sqrt(nGoodPoints)
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


  if chargeRatioVdtHist:
    qMeassRatios = numpy.zeros(nPoints**2)
    qMeassDeltaTimes = numpy.zeros(nPoints**2)
    k = 0
    for i in range(0,nPoints,10):
      for j in range(i,nPoints):
        chargeRatioVdtHist.Fill(ts[j]-ts[i],log(qMeass[j]/qMeass[i]))
        qMeassRatios[k] = qMeass[j]/qMeass[i]
        qMeassDeltaTimes[k] = ts[j]-ts[i]
        k += 1
    qMeassRatios = qMeassRatios[:k]
    qMeassDeltaTimes = qMeassDeltaTimes[:k]
    qMeassRatiosLog = numpy.log(qMeassRatios)

    if False and doPlots:
      fig, ax = mpl.subplots()
      ax.scatter(qMeassDeltaTimes,qMeassRatios,2.,c='k',lw=0)
      ax.set_xlabel("$\Delta t$ [us]")
      ax.set_ylabel("$Q_1/Q_2$")
      ax.set_ylim(0,3.5)
      fig.savefig("Ratios{}.png".format(suffix))
      fig.savefig("Ratios{}.pdf".format(suffix))

      fig, ax = mpl.subplots()
      ax.scatter(qMeassDeltaTimes,qMeassRatiosLog,2.,c='k',lw=0)
      ax.set_xlabel("$\Delta t$ [us]")
      ax.set_ylabel("$\ln(Q_1/Q_2$)")
      ax.set_ylim(-1.5,1.5)
      fig.savefig("RatiosLog{}.png".format(suffix))
      fig.savefig("RatiosLog{}.pdf".format(suffix))

  if doPlots:
    ax.set_ylim(0,qMPV*8)
    ax.scatter(ts,qMeass,2.,c='k',lw=0)
    ax.plot(ts,qMPV*numpy.exp(-ts/lifetimeTrue),'-m')
    ax.plot(ts,numpy.exp(ts*B+A),'-g')
    #ax.plot(ts,numpy.exp(ts*coefs[0]+coefs[1]),'-c')
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
    #ax.plot(ts,ts*coefs[0]+coefs[1],'-c')
    ax.plot(ts,numpy.log(qMPV)-ts/lifetimeTrue,"m-")
    ax.set_xlabel("Drift Time [us]")
    ax.set_ylabel("log(Charge)")
    ax.set_ylim(4.5,7)
    fig.savefig("LandauLog{}.png".format(suffix))
    fig.savefig("LandauLog{}.pdf".format(suffix))

  return 1./lifeInv, numpyLife, numpyLogLife, numpyLifeVar, numpyLogLifeVar

if __name__ == "__main__":
  nBins = 10
  pointsPerBin = 400./nBins
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
  
  #chargeRatioVdt = root.TH2F("chargeRatioVdt","",500,0,2000,500,-1.5,1.5)
  chargeRatioVdt = root.TH2F("chargeRatioVdt","",100,0,1000,100,-1.5,1.5)
  setHistTitles(chargeRatioVdt,"#Delta t [us]","log(Q_{1}/Q_{2})")
  #chargeRatioVdt = None

  lifes = []
  lifesNumpy = []
  lifesLogNumpy = []
  pullsNumpy = []
  pullsLogNumpy = []
  for iCluster in range(10000):
  #for iCluster in range(100):
    doPlots = (iCluster < 5)
    life, lifeNumpy, lifeLogNumpy, lifeNumpyVar, lifeLogNumpyVar = toyCluster(qMPV,lifetimeTrue,nBins,pointsPerBin,usPerBin,suffix="_{}".format(iCluster),doLogFit=doLogFit,doPlots=doPlots,doGaus=doGaus,chargeRatioVdtHist=chargeRatioVdt)
    lifes.append(life/1000.)
    lifesNumpy.append(lifeNumpy/1000.)
    lifesLogNumpy.append(lifeLogNumpy/1000.)
    pullsNumpy.append(lifeNumpy/numpy.sqrt(lifeNumpyVar))
    pullsLogNumpy.append(lifeLogNumpy/numpy.sqrt(lifeLogNumpyVar))

  fig, ax = mpl.subplots()
  ax.hist(lifes,bins=30,range=[0,6],histtype='step')
  ax.hist(lifesNumpy,bins=30,range=[0,6],histtype='step')
  ax.hist(lifesLogNumpy,bins=30,range=[0,6],histtype='step')
  ax.axvline(lifetimeTrue/1000.,c='m')
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

  if chargeRatioVdt:
    #chargeRatioVdt.RebinX(20)
    profileX = chargeRatioVdt.ProfileX()
    projectionX = chargeRatioVdt.ProjectionX()
    projectionY = chargeRatioVdt.ProjectionY()
    setHistTitles(profileX,"#Delta t [us]","<log(Q_{1}/Q_{2})>")
    setHistTitles(projectionX,"#Delta t [us]","Hit Pairs / Bin")
    setHistTitles(projectionY,"log(Q_{1}/Q_{2})","Hit Pairs / Bin")
    canvas = root.TCanvas("c")
    setupCOLZFrame(canvas)
    chargeRatioVdt.Draw("colz")
    profileX.Draw("same")
    drawStandardCaptions(canvas,"Toy Study")
    canvas.SaveAs("ChargeRatioVDeltaT.png")
    setupCOLZFrame(canvas,True)
    projectionX.Draw()
    drawStandardCaptions(canvas,"Toy Study")
    canvas.SaveAs("ChargeRatioVDeltaT_projX.png")
    projectionY.Draw()
    drawStandardCaptions(canvas,"Toy Study")
    canvas.SaveAs("ChargeRatioVDeltaT_projY.png")

    fitfunc = root.TF1("func","pol1",0,5000)
    fitfunc.SetLineColor(root.kBlue)
    #fitrslt = profileX.Fit(fitfunc,"WLFSEM","",50,600)
    #fitrslt = profileX.Fit(fitfunc,"S","",50,600)
    fitrslt = profileX.Fit(fitfunc,"SFM","",50,600)
    profileX.Draw("")
    fitfunc.Draw("same")
    chi2ndf = fitrslt.Chi2()/fitrslt.Ndf()
    slope = fitrslt.Parameter(1)
    slopeErr = fitrslt.ParError(1)
    lifetime = -1.
    lifetimeErr = -1.
    if slope != 0.:
      lifetime = -1./slope/1000.
      lifetimeErr = slopeErr/slope**2/1000.
    drawStandardCaptions(canvas,"Toy Study",
            captionright1="e^{{-}} Lifetime = {:.2f} +/- {:.2f} ms".format(lifetime,lifetimeErr),
            captionright2="#chi^{{2}}/NDF = {:.2f}".format(chi2ndf)
        )
    canvas.SaveAs("ChargeRatioVDeltaT_profX.png")
    
    gausfunc = root.TF1("gausfunc","gaus",-0.75,0.75)
    nBinsX = chargeRatioVdt.GetXaxis().GetNbins()
    muGraph = root.TGraphErrors()
    muList = numpy.zeros(nBinsX)
    dtList = numpy.zeros(nBinsX)
    muErrList = numpy.zeros(nBinsX)
    for iBinX in range(1,nBinsX+1):
      deltaTCenter = chargeRatioVdt.GetXaxis().GetBinCenter(iBinX)
      deltaTLow = chargeRatioVdt.GetXaxis().GetBinLowEdge(iBinX)
      deltaTHigh = chargeRatioVdt.GetXaxis().GetBinUpEdge(iBinX)
      xBinHist = getXBinHist(chargeRatioVdt,iBinX)
      fitrslt = xBinHist.Fit(gausfunc,"S","",-0.75,0.75)
      chi2ndf = fitrslt.Chi2()/fitrslt.Ndf()
      gausmu = fitrslt.Parameter(1)
      gausmuerr = fitrslt.ParError(1)
      gaussigma = fitrslt.Parameter(2)
      gaussigmaerr = fitrslt.ParError(2)
      muGraph.SetPoint(iBinX-1,deltaTCenter,gausmu)
      muGraph.SetPointError(iBinX-1,0.5*(deltaTHigh-deltaTLow),gausmuerr)
      muList[iBinX-1] = gausmu
      muErrList[iBinX-1] = gausmuerr
      dtList[iBinX-1] = deltaTCenter
      if iBinX % (nBinsX // 10) == 0:
        xBinHist.Draw()
        gausfunc.Draw("same")
        setHistTitles(xBinHist,"log(Q_{1}/Q_{2})","Hit Pairs / Bin")
        drawStandardCaptions(canvas,"Toy Study Slice {}, #Delta t: {:.0f} to {:.0f} us".format(iBinX,deltaTLow,deltaTHigh),
                captionright1="#mu = {:.3f} #pm {:.3f}".format(gausmu,gausmuerr),
                captionright2="#sigma = {:.3f} #pm {:.3f}".format(gaussigma,gaussigmaerr),
                captionright3="#chi^{{2}}/NDF = {:.2f}".format(chi2ndf)
            )
        canvas.SaveAs("ChargeRatioVDeltaT_binHist{}.png".format(iBinX))

    axisHist = drawGraphs(canvas,[muGraph],"#Delta t [us]","log(Q_{1}/Q_{2})")
    fitrslt = muGraph.Fit(fitfunc,"S","",500,1000)
    fitfunc.Draw("same")
    chi2ndf = fitrslt.Chi2()/fitrslt.Ndf()
    slope = fitrslt.Parameter(1)
    slopeErr = fitrslt.ParError(1)
    lifetime = -1.
    lifetimeErr = -1.
    if slope != 0.:
      lifetime = -1./slope/1000.
      lifetimeErr = slopeErr/slope**2/1000.
    drawStandardCaptions(canvas,"Toy Study",
            captionright1="e^{{-}} Lifetime = {:.2f} +/- {:.2f} ms".format(lifetime,lifetimeErr),
            captionright2="#chi^{{2}}/NDF = {:.2f}".format(chi2ndf)
        )
    canvas.SaveAs("ChargeRatioVDeltaT_gausFitFit.png")

    import scipy.interpolate
    spline = scipy.interpolate.UnivariateSpline(dtList,muList,w=1/muErrList)
    splineHard = scipy.interpolate.UnivariateSpline(dtList,muList,w=1/muErrList,s=0)
    spline2 = scipy.interpolate.UnivariateSpline(dtList,muList,w=1/muErrList,s=len(muErrList)/3.)
    xs = numpy.linspace(dtList[0],dtList[-1],1000)
    
    fig, ax = mpl.subplots()
    ax.errorbar(dtList,muList,muErrList,fmt='ko',ecolor='k',barsabove=True,markersize=3)
    ax.plot(xs,spline(xs))
    #ax.plot(xs,splineHard(xs))
    ax.plot(xs,spline2(xs))
    ax.set_xlabel("Delta t [us]")
    ax.set_ylabel("Smoothed Fitted $log(Q_1/Q_2)$")
    fig.savefig("ChargeRatioVDeltaT_spline.png")
      
    fig, ax = mpl.subplots()
    ax.plot(xs,-1./spline(xs,1)/1000.)
    #ax.plot(xs,-1./splineHard(xs,1)/1000.)
    ax.plot(xs,-1./spline2(xs,1)/1000.)
    ax.set_xlabel("Delta t [us]")
    ax.set_ylabel("Electron lifetime [ms]")
    ax.set_ylim(0,6)
    fig.savefig("ChargeRatioVDeltaT_splineDeriv.png")
