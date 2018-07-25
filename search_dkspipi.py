#! /usr/bin/env python

import ROOT as r
# import os
import sys
# import glob
# import math
import numpy as np
import array
# import tof as mytof


def find_resonance_2(a, b, m):
    '''
    Find the pair of particles with masses closest to m.
    One particle should be from the input list a, another from the input list b.
    Returns the list:
        [ index of the particle from a,
            index of the particle from b,
            invariant mass of the best pair,
            TLorentzVector of the best pair ]
    '''
    dm = 1.e3
    i, j = -1, -1
    m1 = -1.
    lv = -1
    for ia, va in enumerate(a):
        for ib, vb in enumerate(b):
            lv = va + vb
            m1 = lv.Mag()
            if dm > abs(m1-m):
                dm = abs(m1-m)
                i, j = ia, ib
    return [i, j, m1, lv]


def find_resonance_2_by_2(a, b, m):
    dm = 1.e3
    i, j, k, l = -1, -1
    m1 = -1.
    lv = -1
    for ia, va in enumerate(a):
        for ib, vb in enumerate(b):
            lv = va + vb
            m1 = lv.Mag()
            if dm > abs(m1-m):
                dm = abs(m1-m)
                i, j = ia, ib
    return [i, j, m1, lv]


def find_resonance_ks_pip_pim(ks, pip, pim, m):
    '''
    Find the trio of particles with masses closest to m.
    One particle should be from the input list ks, another from the input list pip,
    and third one should be from the list pim.
    Returns the list:
        [ index of the particle from ks,
            index of the particle from pip,
            index of the particle from pim,
            invariant mass of the best pair,
            TLorentzVector of the best pair ]
    '''
    dm = 1.e3
    i, j, k = -1, -1, -1
    m1 = -1.
    lv = -1
    for ia, va in enumerate(ks):
        for ib, vb in enumerate(pip):
            for ic, vc in enumerate(pim):
                lv = va + vb + vc
                m1 = lv.Mag()
                if dm > abs(m1-m):
                    dm = abs(m1-m)
                    i, j, k = ia, ib, ic
    return [i, j, k, m1, lv]


def main():
    finname = 'output/txt2root.root'#  Input file name.
    Fin = r.TFile( finname )#  Open the input file.
    t = Fin.Get( 't' )#  Get the tree with data from the input file.
    #  Declare meson masses.
    mks = .497611#  GeV / c^2
    md0 = 1.86483#  GeV / c^2
    mpip = 0.13957061#  GeV / c^2
    #  Define output histograms.
    t_mks = 'K_{S} mass M(#pi^{+}#pi^{#minus}), GeV/c^{2}'#  An axis title.
    t_md0 = 'D^{0} mass M(2(#pi^{+}#pi^{#minus})), GeV/c^{2}'
    t_pd0 = 'D^{0} momentum P(2(#pi^{+}#pi^{#minus})), GeV/c'
    hmks = r.TH1D( 'hmks', 'hmks;'+t_mks, 1000, 0, 4 )#  An 1-dimensional histogram.
    hmd0 = r.TH1D( 'hmd0', 'hmd0;'+t_md0, 1000, 0, 4 )
    hmksmd0 = r.TH2D( 'hmksmd0', 'hmksmd0;{};{}'.format(t_mks, t_md0), 140, 0.2, 1.6, 200, 0.5, 3.5 )#  A 2-dimensional histogram.
    hmd0pd0 = r.TH2D( 'hmd0pd0', r'hmd0pd0;{};{}'.format(t_md0, t_pd0), 200, 0.5, 3.5, 200, 0, 2 )
    #  Go throught the tree event by event.
    for eve in t:
        if eve.n < 4:#  Check the total number of reconstucted particles in the current event.
            continue
        #  Check the number of charge particles.
        nq = 0
        for q in eve.q:
            nq += 1
        if nq < 4:
            continue
        #  Create lists with pi+ and pi- mesons lorentz vectors.
        #  https://root.cern.ch/doc/master/classTLorentzVector.html
        lvpip, lvpim = [], []
        for i in xrange( eve.n ):
            if eve.q[i] == -1:
                lvpim.append( r.TLorentzVector() )
                lvpim[-1].SetXYZM( eve.px[i], eve.py[i], eve.pz[i], mpip )
            if eve.q[i] == 1:
                lvpip.append( r.TLorentzVector() )
                lvpip[-1].SetXYZM( eve.px[i], eve.py[i], eve.pz[i], mpip )
        if len( lvpip ) < 2:
            continue
        if len( lvpim ) < 2:
            continue
        #  Find resonances and fill histograms.
        ks = find_resonance_2( lvpip, lvpim, mks )#  Construct a K_S meson from pi+ and pi-.
        hmks.Fill( ks[2] )#  Fill an 1d histogram.
        #
        lvpip1 = lvpip
        del lvpip1[ ks[0] ]#  Drop the used pi+ in K_S.
        lvpim1 = lvpim
        del lvpim1[ ks[1] ]#  Drop the used pi- in K_S.
        d0 = find_resonance_ks_pip_pim( [ks[3]], lvpip1, lvpim1, md0 )#  Construct a D0 meson from K_S, pi+ and pi-.
        hmd0.Fill( d0[3] )
        #
        hmksmd0.Fill( ks[2], d0[3] )#  Fill a 2d histogram.
        hmd0pd0.Fill( d0[3], d0[4].P() )
    ################################
    #  Draw histograms with obtained results.
    r.gStyle.SetPalette(56)
    r.gStyle.SetNumberContours(80)
    r.gStyle.SetOptStat(0)
    cx, cy = 804, 828
    #
    line = r.TLine()
    line.SetLineColor( r.kCyan )
    #
    c1 = r.TCanvas( 'c1', 'c1', cx, cy )
    hmks.Draw()
    c1.Update()
    #
    c2 = r.TCanvas( 'c2', 'c2', cx, cy )
    hmd0.Draw()
    c2.Update()
    #
    c3 = r.TCanvas( 'c3', 'c3', cx, cy )
    c3.SetGrid(1)
    # c3.SetLogz(1)
    hmksmd0.Rebin2D( 4, 4 )
    hmksmd0.GetXaxis().SetTitleOffset(1.35)
    hmksmd0.GetYaxis().SetTitleOffset(1.35)
    gopt = 'colz'
    # gopt = 'lego20z'
    hmksmd0.Draw( gopt )
    line.DrawLine( mks, .5, mks, 3.5 )
    line.DrawLine( .2, md0, 1.6, md0 )
    c3.Update()
    #
    c4 = r.TCanvas( 'c4', 'c4', cx, cy )
    c4.SetGrid(1)
    hmd0pd0.Rebin2D( 4, 4 )
    hmd0pd0.GetXaxis().SetTitleOffset(1.35)
    hmd0pd0.GetYaxis().SetTitleOffset(1.35)
    hmd0pd0.Draw( gopt )
    line.DrawLine( md0, 0, md0, 2 )
    c4.Update()
    #
    if not '-b' in sys.argv:
        raw_input()
    #
    if 'save' in sys.argv:
        c1.SaveAs( 'output/dkspipi_mks.png' )
        c2.SaveAs( 'output/dkspipi_md0.png' )
        c3.SaveAs( 'output/dkspipi_mksmd0.png' )
        c4.SaveAs( 'output/dkspipi_md0pd0.png' )


if __name__ == '__main__':
    main()
