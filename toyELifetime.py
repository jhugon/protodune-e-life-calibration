#!/usr/bin/env python

import ROOT as root
import numpy
from matplotlib import pylab as mpl

nPoints = 1000
qMPV = 300.
lifetimeTrue = 2000. # us
driftSpeed = 1 # us / tick
trackSlope = 1 # tick / point
landauWidth = qMPV*0.22
ticksPerBin = 40
nBins = nPoints // ticksPerBin

ts = numpy.zeros(nPoints)
qTrues = numpy.zeros(nPoints)
qMeass = numpy.zeros(nPoints)

rand = root.TRandom3(7)

tck = numpy.zeros(nBins)
ave = numpy.zeros(nBins)
err = numpy.zeros(nBins)
cnt = numpy.zeros(nBins)

fig, ax = mpl.subplots()

for iPoint in range(nPoints):
    # correct for dumb root MPV!!!!
    param2 = 4.*landauWidth
    param1 = qMPV*0.22278+param2
    qTrue = rand.Landau(qMPV,landauWidth)
    #qTrue = rand.Landau(param1,param2)
    t = iPoint * trackSlope # ticks
    qMeas = qTrue*numpy.exp(-t/lifetimeTrue)

    ts[iPoint] = t
    qTrues[iPoint] = qTrue
    qMeass[iPoint] = qMeas

    ## Now do the fit stuff
    iBin = t // ticksPerBin
    if iBin < nBins:
      tck[iBin] += t
      ave[iBin] += qMeas
      err[iBin] += qMeas*qMeas
      cnt[iBin] += 1

tck /= cnt
ave /= cnt

ax.plot(tck,ave,'or')

maxChg = ave * 1.3
minChg = ave * 0.5

tck = numpy.zeros(nBins)
ave = numpy.zeros(nBins)
err = numpy.zeros(nBins)
cnt = numpy.zeros(nBins)

# now truncate
for iPoint in range(nPoints):
    t = ts[iPoint]
    qMeas = qMeass[iPoint]
    iBin = t / ticksPerBin
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

for iBin in range(nBins):
    print iBin, tck[iBin],ave[iBin],err[iBin],cnt[iBin]

ax.set_ylim(0,qMPV*8)
ax.scatter(ts,qMeass,2.,c='k',lw=0)
ax.plot(ts,qMPV*numpy.exp(-ts/lifetimeTrue),'-g')
#ax.plot(tck,ave,'og')
ax.errorbar(numpy.arange(nBins)*ticksPerBin+0.5*ticksPerBin,ave,xerr=0.5*ticksPerBin,yerr=err,fmt="ob")
ax.set_xlabel("Drift Time [us]")
ax.set_ylabel("Charge")
fig.savefig("Landau.png")
fig.savefig("Landau.pdf")
