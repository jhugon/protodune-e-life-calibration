#!/usr/bin/env python

import ROOT as root
from helpers import *
root.gROOT.SetBatch(True)
import sys

if __name__ == "__main__":

  c = root.TCanvas()

  #f = root.TFile("Lifetime_hist.root")
  #f = root.TFile("Lifetime_10_noSCE.root")
  f = root.TFile("Lifetime_1k_SCE.root")
  title = "1k Events, SCE "
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


  plotHistsSimple([f.Get("lifetime/Life")],None,"Electron Lifetime [ms]","Clusters / Bin",c,"LifetimeAll",captionArgs=[title])

  plotHistsSimple([f.Get("lifetime/ChiDOFZoom")],None,"#chi^2/DOF","Clusters / Bin",c,"Chi2",captionArgs=[title])

  plotHistsSimple([f.Get("lifetime/LifeInvCA"),f.Get("lifetime/LifeInvAC")],["Cathode to Anode","Anode to Cathode"],"1/Lifetime [ms^{-1}]","Normalized Clusters / Bin",c,"LifeInvCAAC",captionArgs=[title],normalize=True)

  plotHistsSimple([f.Get("lifetime/MaxDriftTime")],None,"Electron Drift Time [ms]","Clusters / Bin",c,"MaxDriftTime",captionArgs=[title])

  hists = [f.Get("lifetime/MaxDriftTime_TPC"+str(i)) for i in range(12)]
  goodhists = [hist for hist in hists if hist.GetEntries() > 0]
  labels = ["TPC "+str(i) for i in range(12) if hists[i].GetEntries() > 0]
  plotHistsSimple(goodhists,labels,"Electron Drift Time [ms]","Clusters / Bin",c,"MaxDriftTimePerTPC",captionArgs=[title])

  plotHist2DSimple(f.Get("lifetime/LifeVChargeEff"),"Cluster Charge Efficiency","Electron Lifetime [ms]",c,"LifeVChargeEff",captionArgs=[title])
  plotHist2DSimple(f.Get("lifetime/LifeVZenith"),"True Particle Zenith Angle [deg]","Electron Lifetime [ms]",c,"LifeVZenithAngle",captionArgs=[title])
  plotHist2DSimple(f.Get("lifetime/LifeVAzimuth"),"True Particle Azimuth Angle [deg]","Electron Lifetime [ms]",c,"LifeVAzimuthAngle",captionArgs=[title])
