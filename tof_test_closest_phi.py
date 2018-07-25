#! /usr/bin/env python

import ROOT as r
import os
import sys
# import glob
import math
import datetime
import array
from random import gauss
import random
import numpy as np
import tof as mytof


clight = 299792458.#  m / s

PDGID = {'e+' : -11, -11 : 'e+',
         'e-' : 11, 11 : 'e-',
         'nu_e' : 12, 12 : 'nu_e',
         'gamma' : 22, 22 : 'gamma',
         'mu+' : -13, -13 : 'mu+',
         'mu-' : 13, 13 : 'mu-',
         'pi0' : 111, 111 : 'pi0',
         'pi+' : 211, 211 : 'pi+',
         'pi-' : -211, -211 : 'pi-',
         'KL0' : 130, 130 : 'KL0',
         'KS0' : 310, 310 : 'KS0',
         'K+' : 321, 321 : 'K+',
         'K-' : -321, -321 : 'K-',
         'n0' : 2112, 2112 : 'n0',
         'p+' : 2212, 2212 : 'p+',
         'p-' : -2212, -2212 : 'p-'}


MASS = {-11 : 0.5109989461, 11 : 0.5109989461,
        -12 : 0., 12 : 0.,
        -13 : .1056583745, 13 : .1056583745,
        -14 : 0., 14 : 0.,
        -15 : 1776.86, 17 : 1776.86,
        -16 : 0., 16 : 0.,
        22 : 0.,
        -211 : 139.57061, 211 : 139.57061,
        111 : 134.9770,
        221 : 547.862,
        113 : 775.26,
        -213 : 775.11, 213 : 775.11,
        223 : 782.65,
        333 : 1019.461,
        -321 : 493.677, 321 : 493.677,
        130 : 0.5293,
        310 : 0.5293}
def pdg2m( pdg ):
    if pdg in MASS:
        return MASS[pdg]
    else:
        print 'pdg2m :  Have not such mass (pdg = {}).'.format(pdg)
        return 0.


def main():
    c1 = r.TCanvas('c1', 'c1', 800, 800)
    frame1 = c1.DrawFrame(-7, -4, 13, 4)
    #
    circle = r.TEllipse()
    line = r.TLine()
    point = r.TEllipse()
    lat = r.TLatex()
    lat.SetTextFont(12)
    lat.SetTextAlign(22)
    def drpoint(x, y, col=r.kRed, radius = 0.05):
        point.SetFillColor(col)
        point.SetLineColor(col)
        point.DrawEllipse(x, y, radius, radius, 0, 360, 0)
    vec = r.TArrow()
    def drmom(x, y, px, py, col=r.kRed-6):
        vec.SetAngle(45)
        vec.SetLineWidth(2)
        vec.SetLineColor(col)
        vec.SetFillColor(col)
        vec.DrawArrow(x, y, x+px, y+py, .010)
    #
    line.SetLineColor(r.kGray+1)
    for i in xrange(0, 4):
        line.DrawLine(-6.5, i, 12.75, i)
        line.DrawLine(-6.5, -i, 12.75, -i)
        for j in xrange(-1, 3):
            line.DrawLine(math.pi*2*j, i-.1, math.pi*2*j, i+.1)
            line.DrawLine(math.pi*2*j, -i-.1, math.pi*2*j, -i+.1)
    #
    #
    def drphi(phi1=.2, phi2=5.8, phi=1., wd=1, y0=0.):
        clo_phi = mytof.closest_phi(phi1, phi2, phi, wd)
        drpoint(phi1, y0, r.kBlue)
        drpoint(phi2, y0, r.kRed)
        drpoint(phi1-math.pi*2, y0, r.kBlue-10)
        drpoint(phi2-math.pi*2, y0, r.kRed-10)
        drpoint(phi1+math.pi*2, y0, r.kBlue-10)
        drpoint(phi2+math.pi*2, y0, r.kRed-10)
        drpoint(phi, y0, r.kBlack)
        drpoint(clo_phi[1], y0+0.1, r.kCyan)
        drmom(phi, y0-.1, wd, 0.)
    #
    drphi(1, 5.0, 2.5, 1, 0.)
    drphi(1, 5.0, 2.5, -1, 1.)
    drphi(1, 5.0, 0.5, 1, 2.)
    drphi(1, 5.0, 0.5, -1, 3.)
    drphi(1, 5.0, 6.0, 1, -1.)
    drphi(1, 5.0, 6.0, -1, -2.)
    #
    c1.Update()
    raw_input()


if __name__ == '__main__':
    main()
