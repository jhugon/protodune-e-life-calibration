#!/usr/bin/env python

from helpers import *
#import ROOT as root
import numpy
import sys
from matplotlib import pylab as mpl
import matplotlib.patches
import matplotlib.collections

root.gROOT.SetBatch(True)
RAND = root.TRandom3(7)

rooMsgServiceInstance = root.RooMsgService.instance()
rooMsgServiceInstance.setSilentMode(True);
for iStream in range(rooMsgServiceInstance.numStreams()):
  rooMsgServiceInstance.setStreamStatus(iStream,False);

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
    if not (yerrs is None):
      graph.SetPointError(iPoint,0.,yerrs[iPoint])
  fitrslt = graph.Fit(f,"EX0SQ","",fitrange[0],fitrange[1])
  chi2ndof = -1.
  ndf = fitrslt.Ndf()
  if ndf != 0:
    chi2ndof = fitrslt.Chi2()/ndf
  param = fitrslt.Parameter(1)
  paramErr = fitrslt.ParError(1)
  lifetime = -1./param
  lifetimeErr = paramErr/param**2
  constParam = fitrslt.Parameter(0)
  constParamErr = fitrslt.Parameter(0)
  #print "life: {} +/- {}".format(lifetime,lifetimeErr)

  if suffix:
    axisHist = drawGraphs(canvas,[graph],"Drift Time [us]","Charge",ylims=[0,1000])
    f.Draw("same")
    canvas.SaveAs("ExpFit{}.png".format(suffix))
    del canvas
    del axisHist
  del graph
  return lifetime, lifetimeErr, constParam, constParamErr, chi2ndof

def directFitExpHits(tss,qMeass,suffix="",doPlots=True,usPerBin=100.,qMPV=None,lifetimeTrue=None):
  ts = numpy.array(tss)
  qs = numpy.array(qMeass)
  nHits = len(ts)
  nBins = int(numpy.ceil((ts[-1]-ts[0]) / usPerBin))
  ave = numpy.zeros(nBins)
  cnt = numpy.zeros(nBins)
  for t, qMeas in zip(ts,qMeass):
      ## Now do the fit stuff
      iBin = int(t // usPerBin)
      if iBin < nBins and qMeas < 1500.:
        ave[iBin] += qMeas
        cnt[iBin] += 1
  ave /= cnt
  maxChg = ave * 1.3

  goodHits = numpy.full(nHits,False,dtype=bool)
  for ihist in range(nBins):
    #ts, qMeass
    inBin = numpy.logical_and(ts >= (ihist*usPerBin),ts < ((ihist+1)*usPerBin))
    goodInBin = numpy.logical_and(inBin,qs < maxChg[ihist])
    goodHits = numpy.logical_or(goodHits,goodInBin)

  tsGood = ts[goodHits]
  qsGood = qs[goodHits]
  life, lifeErr, constParam, constParamErr, chi2ndof = rootExpFitPoints(tsGood,qsGood,None,[tsGood[0],tsGood[-1]])

  if doPlots:
    fig, ax = mpl.subplots()
    if qMPV:
      ax.set_ylim(0,qMPV*8)
    patchList = []
    for ihist in range(nBins):
      patchList.append(matplotlib.patches.Rectangle((ihist*usPerBin,0),usPerBin,maxChg[ihist]))
    patchCollection = matplotlib.collections.PatchCollection(patchList,facecolor='0.7',edgecolor='0.7')
    ax.add_collection(patchCollection)
    ax.scatter(ts,qs,15.,c=goodHits,lw=0,cmap="PiYG")
    if qMPV and lifetimeTrue:
      ax.plot(ts,qMPV*numpy.exp(-ts/lifetimeTrue),'-m')
    ax.plot(ts,numpy.exp(-ts/life+constParam),'-g')
    ax.set_xlabel("Drift Time [us]")
    ax.set_ylabel("Charge")
    ax.set_ylim(0,2000)
    ax.margins(x=0.,y=0.)
    fig.savefig("LandauDirect{}.png".format(suffix))
    fig.savefig("LandauDirect{}.pdf".format(suffix))
    mpl.close(fig)

  return life, lifeErr

def rooLandauGausFitter(hist,suffix):
   #hist = Hist(1500,0,1500)
   #for chargeVal in chargeValues:
   #  hist.Fill(chargeVal)

   q = root.RooRealVar("q","Charge",0,1500)
   observables = root.RooArgSet(q)
   observableList = root.RooArgList(q)

   mpv = root.RooRealVar("mpv","mpv landau",300.,0,1500)
   xi = root.RooRealVar("xi","xi landau",50.,0,500.)

   landau = root.RooLandau("landau","landau",q,mpv,xi)

   mg = root.RooRealVar("mg","mg",0)
   sg = root.RooRealVar("sg","sg",100,0,100.)
   gauss = root.RooGaussian("gauss","gauss",q,mg,sg)

   q.setBins(10000,"cache")
   langaus = root.RooFFTConvPdf("langaus","landau (X) gauss",q,landau,gauss)

   model = landau
   nParams = 2

   #data = model.generate(observables,10000)
   data = root.RooDataHist("data","",observableList,hist)

   fitrslt = model.fitTo(data,root.RooFit.Save(True),root.RooFit.Verbose(False),root.RooFit.PrintEvalErrors(-1))

   binChi2 = -9999.
   frame = q.frame(root.RooFit.Title(""))
   data.plotOn(frame)
   model.plotOn(frame)
   binChi2 = frame.chiSquare(nParams)
   if not (suffix is None):

     c = root.TCanvas("rf208_convolution","rf208_convolution",600,600)
     root.gPad.SetLeftMargin(0.15)
     frame.GetYaxis().SetTitleOffset(1.4)
     frame.Draw()
     #frame.Draw("same")
     #axisHist = root.TH2F("axisHist","",1,0,50,1,0,1000)
     ##axisHist = root.TH2F("axisHist","",1,-5,5,1,0,1200)
     #axisHist.Draw()
     #frame.Draw("same")
     frame.SetTitle("")
     frame.GetYaxis().SetTitle("Hits/Bin")
     drawStandardCaptions(c,"{}".format(suffix),
                   captionright1="MPV {:.3f} #pm {:.3f}".format(mpv.getVal(),mpv.getError()),
                   captionright2="Min NLL: {:.2f}".format(fitrslt.minNll()),
                   captionright3="#chi^{{2}}/NDF: {:.2f}".format(binChi2)
           )
     c.SaveAs("roofit{}.png".format(suffix))
   return mpv.getVal(), mpv.getError(), binChi2


def bruceMethod(ts,qMeass,usPerBin=100.,suffix="",doLogFit=False,doPlots=True,chargeRatioVdtHist=None,assumeLinear=False,qMPV=None,lifetimeTrue=None,doRootExpFit=False,doFitBinsAndExp=False):

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
  name = uuid.uuid1().hex
  canvas = None
  binHists = [Hist(50,0,1500) for i in range(nBins)]
  if doPlots:
    fig, ax = mpl.subplots()
    canvas = root.TCanvas(name+"canvas")

  for t, qMeas in zip(ts,qMeass):
      ## Now do the fit stuff
      iBin = int(t // usPerBin)
      if doFitBinsAndExp:
        binHists[iBin].Fill(qMeas)
      else:
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

  rooMPVs = numpy.zeros(nBins)
  rooErrs = numpy.zeros(nBins)
  rooBinChi2s = [0]*nBins
  A = None
  B = None
  lifeInv = None
  lifetime = None
  lifetimeErr = None
  chi2ndof = None
  logPlot_xs = []
  logPlot_ys = []
  logFun_xs = []
  logFun_ys = []

  if doFitBinsAndExp:
    for iBin in range(nBins):
      rooSuffix = None
      if doPlots:
        rooSuffix = "{}_bin{}".format(suffix,iBin)
      rooMPVs[iBin], rooErrs[iBin], rooBinChi2s[iBin] = rooLandauGausFitter(binHists[iBin],suffix=rooSuffix)
    lifetime, lifetimeErr, constParam, constParamErr, chi2ndof = rootExpFitPoints(numpy.arange(nBins)*usPerBin+0.5*usPerBin,rooMPVs,rooErrs,[tck[0],tck[-1]])
    A = constParam
    B = -1./lifetime
  else:
    tck /= cnt
    ave /= cnt
    
    maxChg = ave * 1.3
    minChg = ave * 0.5

    if doPlots:
      ax.plot(tck,ave,'or')
      #gausfunc = root.TF1("gausfunc","gaus",0,1500)
      #for iBin in range(nBins):
      #  fitrslt = binHists[iBin].Fit(gausfunc,"S","",0,ave[iBin]*1.5)
      #  chi2ndof = fitrslt.Chi2()/fitrslt.Ndf()
      #  gausmu = fitrslt.Parameter(1)
      #  gausmuerr = fitrslt.ParError(1)
      #  gaussigma = fitrslt.Parameter(2)
      #  gaussigmaerr = fitrslt.ParError(2)
      #  setHistTitles(binHists[iBin],"Charge","Hits")
      #  binHists[iBin].Sumw2()
      #  binHists[iBin].Draw()
      #  gausfunc.Draw("same")
      #  drawStandardCaptions(canvas,"Cluster {} Bin {}".format(suffix,iBin),
      #                captionright1="Truncated: {:4.0f}-{:4.0f}".format(minChg[iBin],maxChg[iBin]),
      #                captionright2="Ave: {:4.0f}".format(ave[iBin]),
      #                captionright3="T: {:4.0f} Cnt: {:4.0f}".format(tck[iBin],cnt[iBin])
      #        )
      #  canvas.SaveAs("BinHist_{}_bin{}.png".format(suffix,iBin))

    
    
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
    
    if doRootExpFit:
      expSuffix = None
      if doPlots:
        expSuffix = suffix
      lifetime, lifetimeErr, constParam, constParamErr, chi2ndof = rootExpFitPoints(tck,ave,err,[tck[0],tck[-1]],suffix=expSuffix)
      A = constParam
      B = -1./lifetime
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
        fitcnt += 1
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
      #print "B ",B," +/- ",BErr;
      
      lifeInv = -B;
      lifeInvErr = BErr / (B * B) ;
      lifetime = 1./lifeInv
      lifetimeErr = lifeInvErr
      if assumeLinear:
        lifeInv = -B/A
        lifetime = -A/B
        if qMPV:
          lifeInv = -B/qMPV
          lifetime = -qMPV/B
      
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
    if rooMPVs[0] != 0.:
      ax.errorbar(numpy.arange(nBins)*usPerBin+0.5*usPerBin,rooMPVs,yerr=rooErrs,fmt="ob")
    else:
      ax.errorbar(numpy.arange(nBins)*usPerBin+0.5*usPerBin,ave,xerr=0.5*usPerBin,yerr=err,fmt="ob")
    #ax.plot(ts,numpy.exp(ts*coefs[0]+coefs[1]),'-c')
    #ax.plot(tck,ave,'og')
    ax.set_xlabel("Drift Time [us]")
    ax.set_ylabel("Charge")
    ax.set_ylim(0,1000)
    fig.savefig("Landau{}.png".format(suffix))
    fig.savefig("Landau{}.pdf".format(suffix))
    mpl.close(fig)
    
    
    if rooMPVs[0] == 0.:
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
      mpl.close(fig)

  return lifetime, lifetimeErr, chi2ndof, rooBinChi2s


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
    inBin = numpy.logical_and(ts >= (ihist*usPerBin),ts < ((ihist+1)*usPerBin))
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
    numpyLife, lifetimeErr, constParam, constParamErr, chi2ndof = rootExpFitPoints(redot,redoave,redoaveerr,[redot[0],redot[-1]],suffix=expSuffix)
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
    mpl.close(fig)
    
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
    mpl.close(fig)

  return numpyLife, numpy.sqrt(numpyLifeVar)

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
    chi2ndof = fitrslt.Chi2()/fitrslt.Ndf()
    slope = fitrslt.Parameter(1)
    slopeErr = fitrslt.ParError(1)
    lifetime = -1.
    lifetimeErr = -1.
    if slope != 0.:
      lifetime = -1./slope/1000.
      lifetimeErr = slopeErr/slope**2/1000.
    drawStandardCaptions(canvas,"Toy Study",
            captionright1="e^{{-}} Lifetime = {:.2f} +/- {:.2f} ms".format(lifetime,lifetimeErr),
            captionright2="#chi^{{2}}/NDF = {:.2f}".format(chi2ndof)
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
      chi2ndof = fitrslt.Chi2()/fitrslt.Ndf()
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
                captionright3="#chi^{{2}}/NDF = {:.2f}".format(chi2ndof)
            )
        canvas.SaveAs("ChargeRatioVDeltaT_binHist{}.png".format(iBinX))

    axisHist = drawGraphs(canvas,[muGraph],"#Delta t [us]","log(Q_{1}/Q_{2})")
    fitrslt = muGraph.Fit(fitfunc,"QS","",500,1000)
    fitfunc.Draw("same")
    chi2ndof = fitrslt.Chi2()/fitrslt.Ndf()
    slope = fitrslt.Parameter(1)
    slopeErr = fitrslt.ParError(1)
    lifetime = -1.
    lifetimeErr = -1.
    if slope != 0.:
      lifetime = -1./slope/1000.
      lifetimeErr = slopeErr/slope**2/1000.
    drawStandardCaptions(canvas,"Toy Study",
            captionright1="e^{{-}} Lifetime = {:.2f} +/- {:.2f} ms".format(lifetime,lifetimeErr),
            captionright2="#chi^{{2}}/NDF = {:.2f}".format(chi2ndof)
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
      mpl.close(fig)
        
      fig, ax = mpl.subplots()
      ax.plot(xs,-1./spline(xs,1)/1000.)
      #ax.plot(xs,-1./splineHard(xs,1)/1000.)
      ax.plot(xs,-1./spline2(xs,1)/1000.)
      ax.set_xlabel("Delta t [us]")
      ax.set_ylabel("Electron lifetime [ms]")
      ax.set_ylim(0,6)
      fig.savefig("ChargeRatioVDeltaT_splineDeriv.png")
      mpl.close(fig)

if __name__ == "__main__":
  print "Start time: ", datetime.datetime.now().replace(microsecond=0).isoformat(' ')

  nClusters = 1000
  nBins = 10
  pointsPerBin = 400./nBins
  #nBins = 5
  #pointsPerBin = 100./nBins
  usPerBin = 100.
  qMPV = 300.
  lifetimeTrue = 3000. # us
  doLogFit = False
  doGaus = False
  doLinear = False
  doRootExpFit = False
  distType = "Landau"
  distTypeLabel = "Landau Charge "
  if doGaus:
    distType = "Gaus"
    distTypeLabel = "Gaussian Charge "

  caseStr = "{}_bins{}_hitpbin{:.0f}".format(distType,nBins,pointsPerBin)

  crm = ChargeRatioMethod()
  
  lifes = -1e6*numpy.ones(nClusters)
  lifesErr = -1e6*numpy.ones(nClusters)
  lifesChi2 = -1e6*numpy.ones(nClusters) 
  lifesNumpy = -1e6*numpy.ones(nClusters)
  lifesNumpyErr = -1e6*numpy.ones(nClusters)
  lifesMPV = -1e6*numpy.ones(nClusters)
  lifesMPVErr = -1e6*numpy.ones(nClusters)
  lifesMPVChi2 = -1e6*numpy.ones(nClusters) 
  allBinChi2s = []
  lifesDirect = -1e6*numpy.ones(nClusters) 
  lifesDirectErr = -1e6*numpy.ones(nClusters) 
  for iCluster in range(nClusters):
    doPlots = (iCluster < 5)
    #doPlots = False
    ts, qMeass = generateCluster(qMPV,lifetimeTrue,int(nBins*pointsPerBin),pointsPerBin/usPerBin,doGaus,doLinear)
    lifes[iCluster], lifesErr[iCluster], lifesChi2[iCluster], DUMMY = bruceMethod(ts,qMeass,usPerBin,suffix="_"+caseStr+"_{}".format(iCluster),doLogFit=doLogFit,doPlots=doPlots,assumeLinear=doLinear,doRootExpFit=doRootExpFit,qMPV=qMPV,lifetimeTrue=lifetimeTrue)
    lifesNumpy[iCluster], lifesNumpyErr[iCluster] = bruceNumpy(ts,qMeass,usPerBin,suffix="_"+caseStr+"_{}".format(iCluster),doLogFit=doLogFit,assumeLinear=doLinear,doRootExpFit=doRootExpFit,doPlots=doPlots,qMPV=qMPV,lifetimeTrue=lifetimeTrue)
    lifesMPV[iCluster], lifesMPVErr[iCluster], lifesMPVChi2[iCluster], binChi2s = bruceMethod(ts,qMeass,usPerBin,suffix="_"+caseStr+"_MPV_{}".format(iCluster),doPlots=doPlots,assumeLinear=doLinear,doFitBinsAndExp=True,qMPV=qMPV,lifetimeTrue=lifetimeTrue)
    allBinChi2s += binChi2s
    crm.processCluster(ts,qMeass)
    lifesDirect[iCluster], lifesDirectErr[iCluster] = directFitExpHits(ts,qMeass,suffix="_"+caseStr+"_Direct_{}".format(iCluster),doPlots=doPlots)

  fig, ax = mpl.subplots()
  ax.hist(lifes/1000.,bins=50,range=[0,10],histtype='step',label="Bruce Method")
  ax.hist(lifesNumpy/1000.,bins=50,range=[0,10],histtype='step',label="Bruce Numpy")
  ax.hist(lifesDirect/1000.,bins=50,range=[0,10],histtype='step',label="Direct Exp Fit")
  ax.hist(lifesMPV/1000.,bins=50,range=[0,10],histtype='step',label="Fit Bins and Exp")
  ax.axvline(lifetimeTrue/1000.,c='m')
  ax.set_xlabel("Electron Lifetime [ms]")
  ax.set_ylabel("Toy Clusters / Bin")
  ax.legend()
  fig.text(0.15,0.955,"{}Hits: {} Bins: {}, Hits/Bin: {:.1f}".format(distTypeLabel,int(nBins*pointsPerBin),nBins,pointsPerBin),ha='left',va='bottom')
  fig.savefig("ToyLifetime_{}_bins{}_hitpbin{:.0f}.png".format(distType,nBins,pointsPerBin))
  fig.savefig("ToyLifetime_{}_bins{}_hitpbin{:.0f}.pdf".format(distType,nBins,pointsPerBin))
  mpl.close(fig)

  crm.calculate()

  fig, ax = mpl.subplots()
  ax.hist(lifesChi2,bins=50,range=[0,10],histtype='step',label="Bruce")
  ax.hist(lifesMPVChi2,bins=50,range=[0,10],histtype='step',label="Fit Bins and Exp")
  ax.set_xlabel("$\chi^2/NDF$")
  ax.set_ylabel("Toy Clusters / Bin")
  ax.legend()
  fig.text(0.15,0.955,"{}Hits: {} Bins: {}, Hits/Bin: {:.1f}".format(distTypeLabel,int(nBins*pointsPerBin),nBins,pointsPerBin),ha='left',va='bottom')
  fig.savefig("ToyLifetimeChi2_{}_bins{}_hitpbin{:.0f}.png".format(distType,nBins,pointsPerBin))
  fig.savefig("ToyLifetimeChi2_{}_bins{}_hitpbin{:.0f}.pdf".format(distType,nBins,pointsPerBin))
  mpl.close(fig)

  fig, ax = mpl.subplots()
  #ax.hist(lifes/lifesErr,bins=50,range=[-50,50],histtype='step',label="Bruce")
  ax.hist(lifesNumpy/lifesNumpyErr,bins=50,range=[-10,10],histtype='step',label="Bruce Numpy")
  ax.hist(lifesDirect/lifesDirectErr,bins=50,range=[-10,10],histtype='step',label="Direct Exp Fit")
  ax.hist(lifesMPV/lifesMPVErr,bins=50,range=[-10,10],histtype='step',label="Fit Bins and Exp")
  ax.set_xlabel("Pull on Electron Lifetime")
  ax.set_ylabel("Toy Clusters / Bin")
  ax.legend()
  fig.text(0.15,0.955,"{}Hits: {} Bins: {}, Hits/Bin: {:.1f}".format(distTypeLabel,int(nBins*pointsPerBin),nBins,pointsPerBin),ha='left',va='bottom')
  fig.savefig("ToyLifetimePulls_{}_bins{}_hitpbin{:.0f}.png".format(distType,nBins,pointsPerBin))
  fig.savefig("ToyLifetimePulls_{}_bins{}_hitpbin{:.0f}.pdf".format(distType,nBins,pointsPerBin))
  mpl.close(fig)

  #fig, ax = mpl.subplots()
  #hist = ax.hist2d(lifes,lifesErr,bins=30,range=[[0,6],[0,6]])
  #ax.axvline(lifetimeTrue/1000.,c='m')
  #cbar = fig.colorbar(hist[3])
  #fig.text(0.15,0.9,"{}Hits: {} Bins: {}, Hits/Bin: {:.1f}".format(distTypeLabel,int(nBins*pointsPerBin),nBins,pointsPerBin),ha='left')
  #fig.savefig("ToyBruceLifetimeErrVLife_{}_bins{}_hitpbin{:.0f}.png".format(distType,nBins,pointsPerBin))
  #fig.savefig("ToyBruceLifetimeErrVLife_{}_bins{}_hitpbin{:.0f}.pdf".format(distType,nBins,pointsPerBin))
  #mpl.close(fig)

  fig, ax = mpl.subplots()
  hist = ax.hist2d(lifes/1000.,lifesChi2,bins=30,range=[[0,6],[0,15]])
  ax.axvline(lifetimeTrue/1000.,c='m')
  cbar = fig.colorbar(hist[3])
  fig.text(0.15,0.955,"{}Hits: {} Bins: {}, Hits/Bin: {:.1f}".format(distTypeLabel,int(nBins*pointsPerBin),nBins,pointsPerBin),ha='left',va='bottom')
  ax.set_xlabel("Electron Lifetime [ms]")
  ax.set_ylabel("$\chi^2/NDF$")
  fig.savefig("ToyBruceChi2VLife_{}_bins{}_hitpbin{:.0f}.png".format(distType,nBins,pointsPerBin))
  fig.savefig("ToyBruceChi2VLife_{}_bins{}_hitpbin{:.0f}.pdf".format(distType,nBins,pointsPerBin))
  mpl.close(fig)

  fig, ax = mpl.subplots()
  hist = ax.hist2d(lifesMPV/1000.,lifesMPVChi2,bins=30,range=[[0,6],[0,15]])
  ax.axvline(lifetimeTrue/1000.,c='m')
  cbar = fig.colorbar(hist[3])
  fig.text(0.15,0.955,"{}Hits: {} Bins: {}, Hits/Bin: {:.1f}".format(distTypeLabel,int(nBins*pointsPerBin),nBins,pointsPerBin),ha='left',va='bottom')
  ax.set_xlabel("Electron Lifetime [ms]")
  ax.set_ylabel("$\chi^2/NDF$")
  fig.savefig("ToyMPVChi2VLife_{}_bins{}_hitpbin{:.0f}.png".format(distType,nBins,pointsPerBin))
  fig.savefig("ToyMPVChi2VLife_{}_bins{}_hitpbin{:.0f}.pdf".format(distType,nBins,pointsPerBin))
  mpl.close(fig)

  fig, ax = mpl.subplots()
  hist = ax.hist2d(lifesMPV/1000.,lifesMPVErr/1000.,bins=30,range=[[0,6],[0,6]])
  ax.axvline(lifetimeTrue/1000.,c='m')
  cbar = fig.colorbar(hist[3])
  fig.text(0.15,0.955,"{}Hits: {} Bins: {}, Hits/Bin: {:.1f}".format(distTypeLabel,int(nBins*pointsPerBin),nBins,pointsPerBin),ha='left',va='bottom')
  ax.set_xlabel("Electron Lifetime [ms]")
  ax.set_ylabel("Error on Electron Lifetime [ms]")
  fig.savefig("ToyMPVLifetimeErrVLife_{}_bins{}_hitpbin{:.0f}.png".format(distType,nBins,pointsPerBin))
  fig.savefig("ToyLifetimeErrVLife_{}_bins{}_hitpbin{:.0f}.pdf".format(distType,nBins,pointsPerBin))
  mpl.close(fig)

  fig, ax = mpl.subplots()
  hist = ax.hist2d(lifesNumpy/1000.,lifesNumpyErr/1000.,bins=30,range=[[0,6],[0,6]])
  ax.axvline(lifetimeTrue/1000.,c='m')
  cbar = fig.colorbar(hist[3])
  fig.text(0.15,0.955,"{}Hits: {} Bins: {}, Hits/Bin: {:.1f}".format(distTypeLabel,int(nBins*pointsPerBin),nBins,pointsPerBin),ha='left',va='bottom')
  ax.set_xlabel("Electron Lifetime [ms]")
  ax.set_ylabel("Error on Electron Lifetime [ms]")
  fig.savefig("ToyNumpyLifetimeErrVLife_{}_bins{}_hitpbin{:.0f}.png".format(distType,nBins,pointsPerBin))
  fig.savefig("ToyNumpyLifetimeErrVLife_{}_bins{}_hitpbin{:.0f}.pdf".format(distType,nBins,pointsPerBin))
  mpl.close(fig)

  fig, ax = mpl.subplots()
  ax.hist(allBinChi2s,bins=100,range=[0,5],histtype='step')
  ax.set_xlabel("Bin Landau Fit $\chi^2/NDF$")
  ax.set_ylabel("Toy Clusters Bins / Bin")
  fig.text(0.15,0.955,"{}Hits: {} Bins: {}, Hits/Bin: {:.1f}".format(distTypeLabel,int(nBins*pointsPerBin),nBins,pointsPerBin),ha='left',va='bottom')
  fig.savefig("ToyMPVBinChi2_{}_bins{}_hitpbin{:.0f}.png".format(distType,nBins,pointsPerBin))
  fig.savefig("ToyMPVBinChi2_{}_bins{}_hitpbin{:.0f}.pdf".format(distType,nBins,pointsPerBin))
  mpl.close(fig)

  print "End time: ", datetime.datetime.now().replace(microsecond=0).isoformat(' ')
