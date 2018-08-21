#!/usr/bin/env python

from helpers import *
#import ROOT as root
import numpy
import sys
from matplotlib import pylab as mpl

root.gROOT.SetBatch(True)
RAND = root.TRandom3(7)

def generateCluster(qMPV,lifetimeTrue,nHits,hitsPerus,doGaus=False,doLinear=False):
  """
  qMPV is the true charge deposited
  lifetimeTrue is the true lifetime in us
  nHits is the number of hits to generate
  hitsPerus is the number of hits per us
  """
  landauWidth = qMPV*0.22

  ts = numpy.zeros(nHits) # in us
  qTrues = numpy.zeros(nHits) # in ADC
  qMeass = numpy.zeros(nHits) # in ADC
  
  for iHit in range(nHits):
    # correct for dumb root MPV!!!!
    param2 = landauWidth
    param1 = qMPV*0.22278+param2
    qTrue = None
    if doGaus:
      qTrue = RAND.Gaus(qMPV,landauWidth)
    else:
      qTrue = RAND.Landau(qMPV,landauWidth)
      #qTrue = RAND.Landau(param1,param2)
    t = iHit / hitsPerus
    if doLinear:
      qMeas = qTrue-t*qMPV/lifetimeTrue
    else:
      qMeas = qTrue*numpy.exp(-t/lifetimeTrue)
  
    ts[iHit] = t
    qTrues[iHit] = qTrue
    qMeass[iHit] = qMeas

  return ts, qMeass

def rootExpFitPoints(xs,ys,yerrs,fitrange,suffix=None):
  assert(len(xs)==len(ys))
  assert((yerrs is None) or len(xs)==len(yerrs))
  assert(len(fitrange)==2)
  # exp([0]+[1]*x)
  name = uuid.uuid1().hex
  canvas = None
  if suffix:
    canvas = root.TCanvas(name+"canvas")
  f = root.TF1(name,"expo",fitrange[0],fitrange[1])
  graph = root.TGraphErrors()
  for iPoint in range(len(xs)):
    graph.SetPoint(iPoint,xs[iPoint],ys[iPoint])
    if yerrs:
      graph.SetPointError(iPoint,0.,yerrs[iPoint])
  fitrslt = graph.Fit(f,"EX0SQ","",fitrange[0],fitrange[1])
  chi2ndf = -1.
  ndf = fitrslt.Ndf()
  if ndf != 0:
    chi2ndf = fitrslt.Chi2()/ndf
  param = fitrslt.Parameter(1)
  paramErr = fitrslt.ParError(1)
  lifetime = -1./param
  lifetimeErr = paramErr/param**2
  constParam = fitrslt.Parameter(0)
  constParamErr = fitrslt.Parameter(1)
  #print "life: {} +/- {}".format(lifetime,lifetimeErr)

  if suffix:
    axisHist = drawGraphs(canvas,[graph],"Drift Time [us]","Charge",ylims=[0,1000])
    f.Draw("same")
    canvas.SaveAs("ExpFit{}.png".format(suffix))
    del canvas
    del axisHist
  del graph
  return lifetime, lifetimeErr, constParam, constParamErr, chi2ndf

def directFitExpHits(tss,qMeass,suffix="",doPlots=True):
  ts = numpy.array(tss)
  qs = numpy.array(qMeass)
  goodHits = qs < 400.
  ts = ts[goodHits]
  qs = qs[goodHits]
  suffix2=None
  if doPlots:
    suffix2="DirectFit"+suffix
  life, lifeErr, constParam, constParamErr, chi2ndf = rootExpFitPoints(ts,qs,None,[ts[0],ts[-1]],suffix=suffix2)
  return life

def bruceMethod(ts,qMeass,usPerBin=100.,suffix="",doLogFit=False,doPlots=True,chargeRatioVdtHist=None,assumeLinear=False,qMPV=None,lifetimeTrue=None,doRootExpFit=False):

  assert(len(ts)==len(qMeass))
  assert(not (doLogFit and assumeLinear))
  nHits = len(ts)
  nBins = int(numpy.ceil((ts[-1]-ts[0]) / usPerBin))

  tck = numpy.zeros(nBins)
  ave = numpy.zeros(nBins)
  err = numpy.zeros(nBins)
  cnt = numpy.zeros(nBins)

  fig = None
  ax = None
  if doPlots:
    fig, ax = mpl.subplots()
  
  for t, qMeas in zip(ts,qMeass):
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
    for t, qMeas in zip(ts,qMeass):
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
  
  A = None
  B = None
  lifeInv = None
  logPlot_xs = []
  logPlot_ys = []
  logFun_xs = []
  logFun_ys = []

  if doRootExpFit:
    expSuffix = None
    if doPlots:
      expSuffix = suffix
    lifetime, lifetimeErr, constParam, constParamErr, chi2ndf = rootExpFitPoints(tck,ave,err,[tck[0],tck[-1]],suffix=expSuffix)
    lifeInv = 1./lifetime
    A = constParam
    B = -lifeInv
  else:
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
    
    for ihist in range(nBins):
      if(cnt[ihist] < 3):
         continue;
      if(err[ihist] == 0):
        continue;
      xx = (tck[ihist] - tck[0]);
      if doLogFit or assumeLinear:
        yy = ave[ihist]
      else:
        yy = numpy.log(ave[ihist]);
      logPlot_xs.append(xx)
      logPlot_ys.append(yy)
      # error on log(x) = dx / x
      arg = ave[ihist] / err[ihist];
      if assumeLinear:
        arg = 1. / err[ihist];
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
    if assumeLinear:
      lifeInv = -B/A
      if qMPV:
        lifeInv = -B/qMPV
    
    # calculate chisq
    chi2 = 0;
    for ihist in range(nBins):
      if(cnt[ihist] < 3):
         continue;
      if(err[ihist] == 0):
         continue;
      xx = (tck[ihist] - tck[0]);
      yy = numpy.exp(A - xx * lifeInv);
      if assumeLinear:
          yy = A - xx * lifeInv
      arg = (yy - ave[ihist]) / err[ihist];
      chi2 += arg * arg;
      #print "chk ", ihist, " xx ", xx, " yy ", yy, " ave ", ave[ihist], " arg ", arg, "\n";
      logFun_xs.append(xx)
      logFun_ys.append(A-xx*lifeInv)
    chi2ndof = chi2 / ndof;

  if doPlots:
    if qMPV:
      ax.set_ylim(0,qMPV*8)
    ax.scatter(ts,qMeass,2.,c='k',lw=0)
    if qMPV and lifetimeTrue:
      if assumeLinear:
        ax.plot(ts,qMPV*(1-ts/lifetimeTrue),'-m')
      else:
        ax.plot(ts,qMPV*numpy.exp(-ts/lifetimeTrue),'-m')
    if assumeLinear:
      ax.plot(ts,ts*B+A,'-g')
    else:
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
    if qMPV and lifetimeTrue:
      ax.plot(ts,numpy.log(qMPV)-ts/lifetimeTrue,"m-")
    ax.set_xlabel("Drift Time [us]")
    ax.set_ylabel("log(Charge)")
    ax.set_ylim(4.5,7)
    fig.savefig("LandauLog{}.png".format(suffix))
    fig.savefig("LandauLog{}.pdf".format(suffix))

  return 1./lifeInv


def bruceNumpy(ts,qMeass,usPerBin=100.,suffix="",doLogFit=False,doPlots=False,assumeLinear=False,doRootExpFit=False,qMPV=None,lifetimeTrue=None):
  assert(len(ts)==len(qMeass))
  assert( not (doLogFit and assumeLinear))
  nHits = len(ts)
  nBins = int(numpy.ceil((ts[-1]-ts[0]) / usPerBin))

  redoave = numpy.zeros(nBins)
  redoavevariance = numpy.zeros(nBins)
  redoaveerr = numpy.zeros(nBins)
  redot = numpy.zeros(nBins)
  redoavelog = None
  redoaveerrlog = None
  qMeassLog = None
  ltMaxHits = None
  nontruncave = numpy.zeros(nBins)
  if doLogFit:
    qMeassLog = numpy.log(qMeass)
  else:
    ltMaxHits = qMeass < 1500.
  for ihist in range(nBins):
    #ts, qMeass
    inBin = numpy.logical_and(ts > (ihist*usPerBin),ts < ((ihist+1)*usPerBin))
    goodHits = numpy.logical_and(inBin,ltMaxHits)
    
    nGoodPoints = len(qMeass[goodHits])
    redot[ihist] = numpy.mean(ts[goodHits])
    if doLogFit:
      redoavelog[ihist] = numpy.mean(qMeassLog[goodHits])
      redoaveerrlog[ihist] = numpy.std(qMeassLog[goodHits])/numpy.sqrt(nGoodPoints)
    else:
      nontruncave[ihist] = numpy.mean(qMeass[goodHits])
      goodHits2 = numpy.logical_and(qMeass[goodHits] > nontruncave[ihist]*0.5,qMeass[goodHits] < nontruncave[ihist]*1.3)
      redoave[ihist] = numpy.mean(qMeass[goodHits][goodHits2])
      redoavevariance[ihist] = numpy.var(qMeass[goodHits][goodHits2])/nGoodPoints
  if not doLogFit:
    redoavelog = numpy.log(redoave)
    redoaveerr = numpy.sqrt(redoavevariance)
    redoaveerrlog = numpy.abs(redoaveerr/redoave)
  if assumeLinear:
    redoavelog = redoave
    redoaveerr = numpy.sqrt(redoavevariance)
    redoaveerrlog = redoaveerr
  
  coefs = [0.,0.]
  numpyLife = None
  numpyLifeVar = None
  if doRootExpFit:
    expSuffix = None
    if doPlots:
      expSuffix = "Numpy"+suffix
    numpyLife, lifetimeErr, constParam, constParamErr, chi2ndf = rootExpFitPoints(redot,redoave,redoaveerr,[redot[0],redot[-1]],suffix=expSuffix)
    coefs[1] = constParam
    coefs[0] = -1./numpyLife
    numpyLifeVar = lifetimeErr**2
  else:
    coefs, cov = numpy.polyfit(redot,redoavelog,1,w=1./redoaveerrlog,cov=True)
    numpyLife = -1/coefs[0]
    numpyLifeVar = numpy.abs(cov[0][0]/coefs[0]**4)
    if assumeLinear:
      numpyLife = -coefs[1]/coefs[0]
      if qMPV:
        numpyLife = -qMPV/coefs[0]

  if doPlots:
    fig, ax = mpl.subplots()
    if qMPV:
      ax.set_ylim(0,qMPV*8)
    ax.scatter(ts,qMeass,2.,c='k',lw=0)
    if qMPV and lifetimeTrue:
      if assumeLinear:
        ax.plot(ts,qMPV*(1-ts/lifetimeTrue),'-m')
      else:
        ax.plot(ts,qMPV*numpy.exp(-ts/lifetimeTrue),'-m')
    if assumeLinear:
      ax.plot(ts,ts*coefs[0]+coefs[1],'-g')
    else:
      ax.plot(ts,numpy.exp(ts*coefs[0]+coefs[1]),'-g')
    #ax.plot(tck,ave,'og')
    ax.plot(numpy.arange(nBins)*usPerBin+0.5*usPerBin,nontruncave,"or")
    ax.errorbar(numpy.arange(nBins)*usPerBin+0.5*usPerBin,redoave,xerr=0.5*usPerBin,yerr=redoaveerr,fmt="ob")
    ax.set_xlabel("Drift Time [us]")
    ax.set_ylabel("Charge")
    ax.set_ylim(0,1000)
    fig.savefig("NumpyLandau{}.png".format(suffix))
    fig.savefig("NumpyLandau{}.pdf".format(suffix))
    
    fig, ax = mpl.subplots()
    ax.scatter(ts-redot[0],numpy.log(qMeass),2.,c='k',lw=0)
    ax.errorbar(redot,redoavelog,yerr=redoaveerrlog,fmt="ob")
    ax.plot(ts,ts*coefs[0]+coefs[1],'-g')
    if qMPV and lifetimeTrue:
      ax.plot(ts,numpy.log(qMPV)-ts/lifetimeTrue,"m-")
    ax.set_xlabel("Drift Time [us]")
    ax.set_ylabel("log(Charge)")
    ax.set_ylim(4.5,7)
    fig.savefig("NumpyLandauLog{}.png".format(suffix))
    fig.savefig("NumpyLandauLog{}.pdf".format(suffix))

  return numpyLife, numpyLifeVar

class ChargeRatioMethod(object):

  def __init__(self):
    self.chargeRatioVdt = root.TH2F("chargeRatioVdt","",100,0,1000,100,-1.5,1.5)
    setHistTitles(self.chargeRatioVdt,"#Delta t [us]","log(Q_{1}/Q_{2})")

  def processCluster(self,ts,qMeass):
    assert(len(ts)==len(qMeass))
    nHits = len(ts)
  
    for i in range(0,nHits,10):
      for j in range(i,nHits):
        if qMeass[j] > 0 and qMeass[i] > 0:
          try:
            self.chargeRatioVdt.Fill(ts[j]-ts[i],log(qMeass[j]/qMeass[i]))
          except ValueError as e:
            print ts[j], ts[i], qMeass[j], qMeass[i]
            raise e

  def calculate(self):
    chargeRatioVdt = self.chargeRatioVdt
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

    try:
      import scipy.interpolate
    except ImportError:
      print "Couldn't import scipy.interpolate"
    else:
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

if __name__ == "__main__":
  nBins = 10
  pointsPerBin = 400./nBins
  usPerBin = 100.
  qMPV = 300.
  lifetimeTrue = 3000. # us
  doLogFit = False
  doGaus = False
  doLinear = False
  doRootExpFit = False

  crm = ChargeRatioMethod()
  
  lifes = []
  lifesNumpy = []
  lifesLogNumpy = []
  lifesDirect = []
  pullsNumpy = []
  pullsLogNumpy = []
  #for iCluster in range(1000):
  for iCluster in range(1000):
    doPlots = (iCluster < 5)
    #doPlots = False
    ts, qMeass = generateCluster(qMPV,lifetimeTrue,int(nBins*pointsPerBin),pointsPerBin/usPerBin,doGaus,doLinear)
    life = bruceMethod(ts,qMeass,usPerBin,suffix="_{}".format(iCluster),doLogFit=doLogFit,doPlots=doPlots,assumeLinear=doLinear,doRootExpFit=doRootExpFit,qMPV=qMPV,lifetimeTrue=lifetimeTrue)
    lifeNumpy, lifeNumpyVar = bruceNumpy(ts,qMeass,usPerBin,suffix="_{}".format(iCluster),doLogFit=doLogFit,assumeLinear=doLinear,doRootExpFit=doRootExpFit,doPlots=doPlots,qMPV=qMPV,lifetimeTrue=lifetimeTrue)
    lifeDirect = directFitExpHits(ts,qMeass,doPlots=doPlots,suffix="_{}".format(iCluster))
    #crm.processCluster(ts,qMeass)
    lifes.append(life/1000.)
    lifesNumpy.append(lifeNumpy/1000.)
    pullsNumpy.append(lifeNumpy/numpy.sqrt(lifeNumpyVar))
    lifesDirect.append(lifeDirect/1000.)

  fig, ax = mpl.subplots()
  ax.hist(lifes,bins=30,range=[0,6],histtype='step')
  ax.hist(lifesNumpy,bins=30,range=[0,6],histtype='step')
  ax.hist(lifesDirect,bins=30,range=[0,6],histtype='step')
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

  #crm.calculate()
