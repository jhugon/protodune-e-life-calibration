#!/usr/bin/env python

import ROOT as root
from helpers import *
root.gROOT.SetBatch(True)
import sys

import glob

if __name__ == "__main__":

  #for fn in glob.glob("*.root"):
  #  ftest = root.TFile(fn)
  #  ftest.cd("lifetime")
  #  for key in ftest.GetListOfKeys():
  #    obj = key.ReadObj()
  #    print fn, len(obj.GetListOfKeys())
  #sys.exit()

  c = root.TCanvas()

  #f = root.TFile("Lifetime_hist.root")
  #title = "MCC10"
  #f = root.TFile("Lifetime_10_noSCE.root")
  #title = "MCC10, 10 Events, SCE "
  #f = root.TFile("Lifetime_10_SCE.root")
  #f = root.TFile("Lifetime_1k_SCE.root")
  #title = "MCC10, 1k Events, SCE "
  #f = root.TFile("Lifetime_1k_SCE.root")
  #title = "MCC10, 1k Events, No SCE "
  #f = root.TFile("Lifetime_1k_noSCE.root")
  f = root.TFile("Lifetime_100_noSCE.root")
  title = "MCC10, 100 Events, No SCE"
  #f = root.TFile("Lifetime_100_SCE.root")
  #title = "MCC10, 100 Events, SCE"
  for i in range(1,10):
    basename = "lifetime/hitChargeVTick"
    scattername = basename+"_"+str(i)
    firstname = basename+"Fit_"+str(i)+"_0"
    secondname = basename+"Fit_"+str(i)+"_1"
    checkname = basename+"Check_"+str(i)
    scatter = f.Get(scattername)
    if not scatter:
        continue
    firstfit = f.Get(firstname)
    secondfit = f.Get(secondname)
    check = f.Get(checkname)
    basename = "lifetime/fit"
    lifefitname = basename+"_"+str(i)
    lifefitfunname = basename+"Fun_"+str(i)
    lifefit = f.Get(lifefitname)
    lifefitfun = f.Get(lifefitfunname)

    setHistTitles(scatter,"Hit Time Tick","Hit Charge [ADC]")
    scatter.SetTitle("Cluster # "+str(i))
    scatter.SetMarkerStyle(8)
    scatter.SetMarkerSize(0.4)

    firstfit.SetLineWidth(2)
    secondfit.SetLineWidth(2)
    check.SetLineWidth(2)

    firstfit.SetLineColor(root.kRed+1)
    secondfit.SetLineColor(root.kBlue+1)
    check.SetLineColor(root.kGreen+1)

    scatter.Draw("AP")
    firstfit.Draw("PEZ")
    secondfit.Draw("PEZ")
    check.Draw("L")
    c.SaveAs("Cluster_"+str(i)+".png")

    setHistTitles(lifefit,"Hit Time [ms]","ln(Hit Charge) [ln(ADC)]")
    lifefit.SetTitle("Cluster # "+str(i))

    lifefit.SetLineWidth(2)
    lifefitfun.SetLineWidth(2)

    lifefit.SetLineColor(root.kBlue+1)
    lifefitfun.SetLineColor(root.kGreen+1)

    lifefit.Draw("AP")
    lifefitfun.Draw("L")
    c.SaveAs("ClusterFit_"+str(i)+".png")

    firstfitYs = firstfit.GetY()
    secondfitYs = secondfit.GetY()
    secondfitErrYs = secondfit.GetEY()
    scatterXs = scatter.GetX()
    scatterYs = scatter.GetY()
    for iBin in range(firstfit.GetN()):
      xMin = iBin*100.
      xMax = (iBin+1)*100.
      binHist = Hist(75,0,1500)
      binHist.Sumw2()
      nHits = 0
      nHitsKept = 0
      sumY = 0.
      for iPoint in range(scatter.GetN()):
        xPoint = scatterXs[iPoint]
        if xPoint >= xMin and xPoint < xMax:
          yPoint = scatterYs[iPoint]
          binHist.Fill(yPoint)
          nHits += 1
          sumY += yPoint
          if yPoint > firstfitYs[iBin]*0.5 and yPoint < firstfitYs[iBin]*1.3:
            nHitsKept += 1
      axisHist = makeStdAxisHist([binHist],freeTopSpace=0.1,includeErrorBar=True)
      axisHist.Draw()
      vtruncate = drawVSpan(axisHist,firstfitYs[iBin]*0.5,firstfitYs[iBin]*1.3)
      vtruncate.SetFillColor(root.kGray)
      #verror = drawVSpan(axisHist,secondfitYs[iBin]-secondfit.GetErrorY(iBin),secondfitYs[iBin]+secondfit.GetErrorY(iBin))
      #verror.SetFillColor(root.kCyan)
      vline = drawVline(axisHist,firstfitYs[iBin])
      vline.SetLineColor(root.kRed+1)
      vline2 = drawVline(axisHist,secondfitYs[iBin])
      vline2.SetLineColor(root.kBlue+1)
      vline3 = drawVline(axisHist,sumY/nHits)
      vline3.SetLineColor(root.kOrange+1)
      binHist.Draw("same")
      setHistTitles(axisHist,"Hit Charge [ADC]","Hits / Bin")
      drawStandardCaptions(c,"Cluster {}, Time Bin {}, {} to {} us".format(i,iBin,xMin,xMax),
                            captionright1="Cluster Hits: {}".format(scatter.GetN()),
                            captionright2="Bin Hits: {}".format(nHits),
                            captionright3="Bin Fraction Kept: {:.3f}".format(float(nHitsKept)/nHits)
        )
      c.SaveAs("BinHist_cluster"+str(i)+"_bin"+str(iBin)+".png")
      c.Clear()

  f.cd("lifetime")
  f.ls()
  plotHistsSimple([f.Get("lifetime/Life")],None,"Electron Lifetime [ms]","Clusters / Bin",c,"LifetimeAll",captionArgs=[title],rebin=None)

  plotHistsSimple([f.Get("lifetime/ChiDOFZoom")],None,"#chi^2/DOF","Clusters / Bin",c,"Chi2",captionArgs=[title])

  plotHistsSimple([f.Get("lifetime/LifeInvCA"),f.Get("lifetime/LifeInvAC")],["Cathode to Anode","Anode to Cathode"],"1/Lifetime [ms^{-1}]","Normalized Clusters / Bin",c,"LifeInvCAAC",captionArgs=[title],normalize=True)

  plotHistsSimple([f.Get("lifetime/MaxDriftTime")],None,"Electron Drift Time [ms]","Clusters / Bin",c,"MaxDriftTime",captionArgs=[title])

  #hists = [f.Get("lifetime/MaxDriftTime_TPC"+str(i)) for i in range(12)]
  #goodhists = [hist for hist in hists if hist.GetEntries() > 0]
  #labels = ["TPC "+str(i) for i in range(12) if hists[i].GetEntries() > 0]
  #plotHistsSimple(goodhists,labels,"Electron Drift Time [ms]","Clusters / Bin",c,"MaxDriftTimePerTPC",captionArgs=[title])

  plotHist2DSimple(f.Get("lifetime/LifeVChargeEff"),"Cluster Charge Efficiency","Electron Lifetime [ms]",c,"LifeVChargeEff",captionArgs=[title])
  plotHist2DSimple(f.Get("lifetime/LifeVZenith"),"True Particle Zenith Angle [deg]","Electron Lifetime [ms]",c,"LifeVZenithAngle",captionArgs=[title])
  plotHist2DSimple(f.Get("lifetime/LifeVAzimuth"),"True Particle Azimuth Angle [deg]","Electron Lifetime [ms]",c,"LifeVAzimuthAngle",captionArgs=[title])

  plotHist2DSimple(f.Get("lifetime/LifeVNHits"),"Number of Hits in Cluster","Electron Lifetime [ms]",c,"LifeVNHits",captionArgs=[title])
  plotHist2DSimple(f.Get("lifetime/LifeVNBins"),"Number of Histogram Bins","Electron Lifetime [ms]",c,"LifeVNBins",captionArgs=[title])
  plotHist2DSimple(f.Get("lifetime/LifeVMaxDriftTime"),"Electron Drift Time [ms]","Electron Lifetime [ms]",c,"LifeVMaxDriftTime",captionArgs=[title])

  plotHist2DSimple(f.Get("lifetime/FracSelHitsVDriftTime"),"Electron Drift Time [ms]","Fraction of Hits Selected",c,"FracSelHitsVDriftTime",captionArgs=[title],profileX=True)

  plotHist2DSimple(f.Get("lifetime/NHitsPerBinVNHits"),None,None,c,"NHitsPerBinVNHits",captionArgs=[title],profileX=True)
  plotHist2DSimple(f.Get("lifetime/NHitsPerBinVNBins"),None,None,c,"NHitsPerBinVNBins",captionArgs=[title],profileX=True)
  plotHist2DSimple(f.Get("lifetime/LifeVMinNHitsPerBin"),None,None,c,"LifeVMinNHitsPerBin",captionArgs=[title],profileX=True)
