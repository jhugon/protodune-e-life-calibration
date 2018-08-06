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
ticksPerBin = 60
nBins = nPoints // ticksPerBin
msPerTick = 1.

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
    if iBin < nBins and qMeas < 1500.:
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
  xx = (tck[ihist] - tck[0]) * msPerTick;
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
print "delta ",delta
A = (sumx2 * sumy - sumx * sumxy) / delta;
print "A ", A
# slope = 1 / lifetime
B = (sumxy * Sum  - sumx * sumy) / delta;
# calculate the error
ndof = fitcnt - 2;
varnce = (sumy2 + A*A*Sum + B*B*sumx2 - 2 * (A*sumy + B*sumxy - A*B*sumx)) / ndof;
BErr = numpy.sqrt(varnce * Sum / delta);
print "B ",B," ",BErr;

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
  xx = (tck[ihist] - tck[0]) * msPerTick;
  yy = numpy.exp(A - xx * lifeInv);
  arg = (yy - ave[ihist]) / err[ihist];
  chi2 += arg * arg;
  print "chk ", ihist, " xx ", xx, " yy ", yy, " ave ", ave[ihist], " arg ", arg, "\n";
  logFun_xs.append(xx)
  logFun_ys.append(A-xx*lifeInv)
chi2ndof = chi2 / ndof;

ax.set_ylim(0,qMPV*8)
ax.scatter(ts,qMeass,2.,c='k',lw=0)
ax.plot(ts,qMPV*numpy.exp(-ts/lifetimeTrue),'-g')
#ax.plot(tck,ave,'og')
ax.errorbar(numpy.arange(nBins)*ticksPerBin+0.5*ticksPerBin,ave,xerr=0.5*ticksPerBin,yerr=err,fmt="ob")
ax.set_xlabel("Drift Time [us]")
ax.set_ylabel("Charge")
fig.savefig("Landau.png")
fig.savefig("Landau.pdf")


fig, ax = mpl.subplots()
ax.plot(logPlot_xs,logPlot_ys,"bo")
ax.plot(logFun_xs,logFun_ys,"g-")
ax.plot(ts,numpy.log(qMPV)-ts/lifetimeTrue,"m-")
ax.set_xlabel("Drift Time [us]")
ax.set_ylabel("log(Charge)")
fig.savefig("LandauLog.png")
fig.savefig("LandauLog.pdf")
