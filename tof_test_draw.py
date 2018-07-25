#! /usr/bin/env python

import ROOT as r
# import os
import sys
import glob
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


def readfile(finname):
    f = open(finname)
    res = {}
    for line in f:
        lin = line.split()
        if len(lin) == 2:
            if lin[0] in ['tof.get_time_v2.x0to2', 'tof.get_time_v2.y0to2', 'tof.get_time_v2.z0to2', 'tof.get_time_v2.rho0to2']:
                isArr = True
            else:
                isArr = False
            try:
                if isArr:
                    res[lin[0]] = [ float(lin[1]) ]
                else:
                    res[lin[0]] = float(lin[1])
            except ValueError:
                res[lin[0]] = lin[1]
        elif len(lin) > 2:
            res[lin[0]] = [float(x) for x in lin[1:]]
    f.close()
    for v in res:
        print '{:<30}'.format(v), res[v]
    return res


def test(finname = 'tof_test_1.txt'):
    k = 2
    if '-b' in sys.argv:
        k = 2
    c1 = r.TCanvas('c1', 'c1', int(k*1200), int(k*800))
    c1.Divide(3, 2)
    #
    circle = r.TEllipse()
    line = r.TLine()
    point = r.TEllipse()
    lat = r.TLatex()
    lat.SetTextFont(12)
    lat.SetTextAlign(22)
    def drpoint(x, y, col=r.kRed, radius = 0.025):
        point.SetFillColor(col)
        point.SetLineColor(col)
        point.DrawEllipse(x, y, radius, radius, 0, 360, 0)
    vec = r.TArrow()
    def drmom(x, y, px, py, col=r.kRed-6):
        vec.SetAngle(45)
        vec.SetLineWidth(2)
        vec.SetLineColor(col)
        vec.SetFillColor(col)
        if (px**2 + py**2) > .5**2:
            vec.DrawArrow(x, y, x+px/2., y+py/2., .005)
        elif (px**2 + py**2) > .2**2:
            vec.DrawArrow(x, y, x+px, y+py, .005)
        else:
            vec.DrawArrow(x, y, x+px*5, y+py*5, .005)
    #
    #
    pad1 = c1.cd(1)
    a, b = 3.1, 3.1
    frame1 = pad1.DrawFrame(-a, -a, a, a)
    #
    par = readfile( finname )
    if 'particle.get_lar.r_lar' in par:
        r_lar = par['particle.get_lar.r_lar']
        lar_cen = par['particle.get_lar.centre']
        h_lar = par['particle.get_lar.h_lar']
    rfield = par['field.init.r']
    poi = par['tof.get_time_v2.poi']
    if 'tof.get_time_v2.cropoi' in par:
        cropoi = par['tof.get_time_v2.cropoi']
    else:
        cropoi = None
    tofrin = 0.5
    tofrout = 2.0
    tofz = 2.0
    vertex = par['particle.init.vertex']
    ################################################
    #
    frame1.Draw('a')
    lat.DrawLatexNDC(.5, .95, 'x-y')
    circle.DrawEllipse(0., 0., tofrout, tofrout, 0, 360, 0)
    circle.SetFillColor(r.kYellow - 10)
    circle.DrawEllipse(0., 0., rfield, rfield, 0, 360, 0)
    circle.SetFillStyle(0)
    if 'particle.get_lar.r_lar' in par:
        circle.DrawEllipse(lar_cen[0], lar_cen[1], r_lar, r_lar, 0, 360, 0)
    if cropoi:
        line.DrawLine(cropoi[0], cropoi[1], cropoi[2], cropoi[3])
    drpoint(0., 0., r.kBlack)# detector centre point
    if 'particle.get_lar.r_lar' in par:
        drpoint(lar_cen[0], lar_cen[1], r.kBlack)# Larmor centre point
    #
    if 'tof.get_time_v2.trk_x' in par:
        trkx = par['tof.get_time_v2.trk_x']
        trky = par['tof.get_time_v2.trk_y']
        for i, v in enumerate(trkx):
            drpoint(v, trky[i], r.kGreen, 0.01)
    #
    line.SetLineColor(r.kBlack)
    line.SetLineWidth(1)
    line.DrawLine(poi[0], poi[1], par['tof.get_time_v2.cropoi_true'][0], par['tof.get_time_v2.cropoi_true'][1])
    #
    drpoint(vertex[0], vertex[1])
    drpoint(poi[0], poi[1])
    if 'tof.get_time_v2.x2' in par:
        drpoint(par['tof.get_time_v2.x2'], par['tof.get_time_v2.y2'])
    drmom(vertex[0], vertex[1], par['particle.init.p'][0], par['particle.init.p'][1])
    drpoint(par['tof.get_time_v2.cropoi_true'][0], par['tof.get_time_v2.cropoi_true'][1])
    if 'tof.get_time_v2.px1_prelim' in par:
        drmom(par['tof.get_time_v2.poi_prelim'][0], par['tof.get_time_v2.poi_prelim'][1],
                par['tof.get_time_v2.px1_prelim'], par['tof.get_time_v2.py1_prelim'])
    if 'tof.get_time_v2.px1' in par:
        drmom(poi[0], poi[1], par['tof.get_time_v2.px1'], par['tof.get_time_v2.py1'])
    #
    if 'particle.get_lar.x' in par:
        line.SetLineColor(r.kOrange)
        line.SetLineWidth(2)
        line.DrawLine(vertex[0], vertex[1],
                vertex[0] + par['particle.get_lar.x']*par['particle.get_lar.x.sign'], vertex[1])
        line.DrawLine(vertex[0], vertex[1],
                vertex[0], vertex[1] + par['particle.get_lar.y']*par['particle.get_lar.y.sign'])
    #
    if 'tof.get_time_v2.x0to2' in par:
        xx =  par['tof.get_time_v2.x0to2']
        yy =  par['tof.get_time_v2.y0to2']
        for i, v in enumerate(xx):
            drpoint(v, yy[i], r.kCyan, 0.01)
    #
    if 'tof.get_time_v2.velx0' in par and 'tof.get_time_v2.vely0' in par:
        vely0 = par['tof.get_time_v2.vely0']
        velx0 = par['tof.get_time_v2.velx0']
        velxy0 = (vely0**2+velx0**2)**.5
        drmom(vertex[0], vertex[1], velx0/velxy0, vely0/velxy0, r.kBlue)
    ################################################
    pad2 = c1.cd(2)
    frame2 = pad2.DrawFrame(-b, -a, b, a)
    frame2.Draw('a')
    lat.DrawLatexNDC(.5, .95, '#rho-z')
    box = r.TBox()
    box.DrawBox(-par['tof.z'], -par['tof.r_out'], par['tof.z'], par['tof.r_out'])
    box.SetFillColor(r.kYellow - 10)
    box.DrawBox(-par['field.init.z'], -par['field.init.r'], par['field.init.z'], par['field.init.r'])
    drpoint(0., 0., r.kBlack)# centre point
    line.SetLineColor(r.kGray+1)
    line.SetLineWidth(2)
    line.DrawLine(-par['tof.z'], par['tof.r_in'], -par['tof.z'], par['tof.r_out'])
    line.DrawLine(par['tof.z'], par['tof.r_in'], par['tof.z'], par['tof.r_out'])
    line.DrawLine(-par['tof.z'], -par['tof.r_in'], -par['tof.z'], -par['tof.r_out'])
    line.DrawLine(par['tof.z'], -par['tof.r_in'], par['tof.z'], -par['tof.r_out'])
    line.DrawLine(-par['tof.z'], -par['tof.r_out'], par['tof.z'], -par['tof.r_out'])
    line.DrawLine(-par['tof.z'], par['tof.r_out'], par['tof.z'], par['tof.r_out'])
    #
    if 'tof.get_time_v2.rho2' in par:
        rho2 = par['tof.get_time_v2.rho2']
    else:
        rho2 = (par['tof.get_time_v2.cropoi_true'][0]**2 + par['tof.get_time_v2.cropoi_true'][1]**2)**.5
    line.DrawLine(vertex[2], (vertex[0]**2 + vertex[1]**2)**.5, par['tof.get_time_v2.z1'], (poi[0]**2 + poi[1]**2)**.5)
    line.DrawLine(par['tof.get_time_v2.z1'], (poi[0]**2 + poi[1]**2)**.5, par['tof.get_time_v2.z2'], rho2)
    drmom(vertex[2], (vertex[0]**2 + vertex[1]**2)**.5, par['particle.init.p'][2], (par['particle.init.p'][0]**2 + par['particle.init.p'][1]**2)**.5)
    drpoint(vertex[2], (vertex[0]**2 + vertex[1]**2)**.5)# initial point
    drpoint(par['tof.get_time_v2.z1'], (poi[0]**2 + poi[1]**2)**.5)# Solenoid hit
    drpoint(par['tof.get_time_v2.z2'], rho2)# ToF hit
    if 'tof.get_time_v2.prho1' in par:
        drmom(par['tof.get_time_v2.z1'], (poi[0]**2 + poi[1]**2)**.5, 0., par['tof.get_time_v2.prho1'])
    #
    if 'tof.get_time_v2.z0to2' in par:
        zz = par['tof.get_time_v2.z0to2']
        rhorho =  par['tof.get_time_v2.rho0to2']
        for i, v in enumerate(zz):
            drpoint(v, rhorho[i], r.kCyan, 0.01)
    ################################################
    pad3 = c1.cd(3)
    frame3 = pad3.DrawFrame(-b, -a, b, a)
    frame3.Draw('a')
    #
    lat.DrawLatexNDC(.5, .95, 'y-z')
    box = r.TBox()
    box.DrawBox(-par['tof.z'], -par['tof.r_out'], par['tof.z'], par['tof.r_out'])
    box.SetFillColor(r.kYellow - 10)
    box.DrawBox(-par['field.init.z'], -par['field.init.r'], par['field.init.z'], par['field.init.r'])
    drpoint(0., 0., r.kBlack)# centre point
    line.SetLineColor(r.kGray+1)
    line.SetLineWidth(2)
    line.DrawLine(-par['tof.z'], par['tof.r_in'], -par['tof.z'], par['tof.r_out'])
    line.DrawLine(par['tof.z'], par['tof.r_in'], par['tof.z'], par['tof.r_out'])
    line.DrawLine(-par['tof.z'], -par['tof.r_in'], -par['tof.z'], -par['tof.r_out'])
    line.DrawLine(par['tof.z'], -par['tof.r_in'], par['tof.z'], -par['tof.r_out'])
    line.DrawLine(-par['tof.z'], -par['tof.r_out'], par['tof.z'], -par['tof.r_out'])
    line.DrawLine(-par['tof.z'], par['tof.r_out'], par['tof.z'], par['tof.r_out'])
    #
    if 'tof.get_time_v2.z0to2' in par:
        zz = par['tof.get_time_v2.z0to2']
        yy =  par['tof.get_time_v2.y0to2']
        for i, v in enumerate(zz):
            drpoint(v, yy[i], r.kCyan, 0.01)
    #
    drpoint(vertex[2], vertex[1])# initial point
    drpoint(par['tof.get_time_v2.z1'], poi[1])# Solenoid hit
    drpoint(par['tof.get_time_v2.z2'], par['tof.get_time_v2.cropoi_true'][1])# ToF hit
    drmom(vertex[2], vertex[1], par['particle.init.p'][2], par['particle.init.p'][1])
    if 'tof.get_time_v2.velz1' in par and 'tof.get_time_v2.vely1' in par:
        vely1 = par['tof.get_time_v2.vely1']
        velz1 = par['tof.get_time_v2.velz1']
        velyz1 = (vely1**2+velz1**2)**.5
        drmom(par['tof.get_time_v2.z1'], poi[1], velz1/velyz1, vely1/velyz1, r.kBlue)# Solenoid hit
    if 'tof.get_time_v2.velz0' in par and 'tof.get_time_v2.vely0' in par:
        vely0 = par['tof.get_time_v2.vely0']
        velz0 = par['tof.get_time_v2.velz0']
        velyz0 = (vely0**2+velz0**2)**.5
        drmom(vertex[2], vertex[1], velz0/velyz0, vely0/velyz0, r.kBlue)# Solenoid hit
    ################################################
    pad4 = c1.cd(4)
    # frame4 = pad4.DrawFrame(-b, -a, b, a)
    frame4 = pad4.DrawFrame(-5, -5, 5, 5)
    frame4.Draw('a')
    co = math.cos( math.pi / 6 )
    si = math.sin( math.pi / 6 )
    #
    line.SetLineColor(r.kGray)
    line.SetLineWidth(2)
    q11 = (-5)*co
    q12 = (5)*co
    w11 = (-5)*si
    w12 = (5)*si
    line.DrawLine( q11, w11, q12, w12 )
    q21 = q22 = 0.
    w21 = -5.
    w22 = 5.
    line.DrawLine( q21, w21, q22, w22 )
    q31 = (-5)*co
    q32 = (5)*co
    w31 = -(-5)*si
    w32 = -(5)*si
    line.DrawLine( q31, w31, q32, w32 )
    #
    if 'tof.get_time_v2.z0to2' in par:
        zz = par['tof.get_time_v2.z0to2']
        yy = par['tof.get_time_v2.y0to2']
        xx = par['tof.get_time_v2.x0to2']
        for i, v in enumerate(zz):
            q = xx[i]*co + zz[i]*co
            w = yy[i] + xx[i]*si - zz[i]*si
            drpoint(q, w, r.kCyan, 0.01)
    if 'tof.get_time_v2.velz0' in par and 'tof.get_time_v2.vely0' in par:
        velx0 = par['tof.get_time_v2.velx0']
        vely0 = par['tof.get_time_v2.vely0']
        velz0 = par['tof.get_time_v2.velz0']
        q0 = vertex[0]*co + vertex[2]*co
        w0 = vertex[1] + vertex[0]*si - vertex[2]*si
        velq0 = velx0*co + velz0*co
        velw0 = vely0 + velx0*si - velz0*si
        velqw0 = (velq0**2 + velw0**2)**.5
        drmom(q0, w0, velq0/velqw0, velw0/velqw0, r.kBlue)
    drpoint( vertex[0]*co + vertex[2]*co, vertex[1] + vertex[0]*si - vertex[2]*si )
    drpoint( vertex[0]*co + vertex[2]*co, vertex[1] + vertex[0]*si - vertex[2]*si )
    x1, y1, z1 = poi[0], poi[1], par['tof.get_time_v2.z1']
    drpoint( (x1+z1)*co, y1+(x1-z1)*si )
    x2, y2, z2 = par['tof.get_time_v2.cropoi_true'][0], par['tof.get_time_v2.cropoi_true'][1], par['tof.get_time_v2.z2']
    drpoint( (x2+z2)*co, y2+(x2-z2)*si )
    ################################################
    pad5 = c1.cd(5)
    frame5 = pad5.DrawFrame(-b, -a, b, a)
    frame5.Draw('a')
    #
    lat.DrawLatexNDC(.5, .95, 'x-z')
    box = r.TBox()
    box.DrawBox(-par['tof.z'], -par['tof.r_out'], par['tof.z'], par['tof.r_out'])
    box.SetFillColor(r.kYellow - 10)
    box.DrawBox(-par['field.init.z'], -par['field.init.r'], par['field.init.z'], par['field.init.r'])
    drpoint(0., 0., r.kBlack)# centre point
    line.SetLineColor(r.kGray+1)
    line.SetLineWidth(2)
    line.DrawLine(-par['tof.z'], par['tof.r_in'], -par['tof.z'], par['tof.r_out'])
    line.DrawLine(par['tof.z'], par['tof.r_in'], par['tof.z'], par['tof.r_out'])
    line.DrawLine(-par['tof.z'], -par['tof.r_in'], -par['tof.z'], -par['tof.r_out'])
    line.DrawLine(par['tof.z'], -par['tof.r_in'], par['tof.z'], -par['tof.r_out'])
    line.DrawLine(-par['tof.z'], -par['tof.r_out'], par['tof.z'], -par['tof.r_out'])
    line.DrawLine(-par['tof.z'], par['tof.r_out'], par['tof.z'], par['tof.r_out'])
    #
    if 'tof.get_time_v2.z0to2' in par:
        zz = par['tof.get_time_v2.z0to2']
        xx =  par['tof.get_time_v2.x0to2']
        for i, v in enumerate(zz):
            drpoint(v, xx[i], r.kCyan, 0.01)
    #
    drpoint(vertex[2], vertex[0])# initial point
    drpoint(par['tof.get_time_v2.z1'], poi[0])# Solenoid hit
    drpoint(par['tof.get_time_v2.z2'], par['tof.get_time_v2.cropoi_true'][0])# ToF hit
    drmom(vertex[2], vertex[0], par['particle.init.p'][2], par['particle.init.p'][0])
    if 'tof.get_time_v2.velz1' in par and 'tof.get_time_v2.velx1' in par:
        velx1 = par['tof.get_time_v2.velx1']
        velz1 = par['tof.get_time_v2.velz1']
        velxz1 = (velx1**2+velz1**2)**.5
        drmom(par['tof.get_time_v2.z1'], poi[0], velz1/velxz1, velx1/velxz1, r.kBlue)
    if 'tof.get_time_v2.velz0' in par and 'tof.get_time_v2.velx0' in par:
        velx0 = par['tof.get_time_v2.velx0']
        velz0 = par['tof.get_time_v2.velz0']
        velxz0 = (velx0**2+velz0**2)**.5
        drmom(vertex[2], vertex[0], velz0/velxz0, velx0/velxz0, r.kBlue)
    ################################################
    c1.cd(6)
    lat.SetTextSize(.04)
    if 'tof.get_time_v2.case' in par:
        lat.DrawLatexNDC(.5, .95, par['tof.get_time_v2.case'])
    lat.DrawLatexNDC(.5, .90, 't_{0} = %.4g sec' % par['tof.get_time_v2.time'][0])
    lat.DrawLatexNDC(.5, .85, 't_{1} = %.4g sec' % par['tof.get_time_v2.time'][1])
    lat.DrawLatexNDC(.5, .80, 'charge = %.1f' % par['particle.init.charge'])
    lat.DrawLatexNDC(.5, .75, 'mass = %.5f GeV/c^{2}' % par['particle.init.mass'])
    ################################################
    c1.Update()
    if not '-b' in sys.argv:
        raw_input()
    if 'save' in sys.argv:
        # c1.SaveAs( '~/Downloads/test_img/{}.png'.format(finname[:-4]) )
        c1.SaveAs( '~/Downloads/{}.png'.format(finname[:-4]) )


def test1():
    test('tof_test_1.txt')

def test2():
    test('tof_test_2.txt')

def test3():
    test('tof_test_3.txt')

def test4():
    test('tof_test_4.txt')

def main():
    '''
    test1()
    test('tof_test_1a.txt')
    test2()
    test3()
    test4()
    test('tof_test_5.txt')
    test('tof_test_6.txt')
    test('tof_test_7.txt')
    test('tof_test_8.txt')
    test('tof_test_9.txt')
    test('tof_test_9b.txt')
    test('tof_test_9c.txt')
    test('tof_test_9d.txt')
    test('tof_test_9e.txt')
    test('tof_test_9f.txt')
    '''
    files = glob.glob('tof_test_*.txt')
    for name in files:
        test(name)


if __name__ == '__main__':
    doMain = True
    for i, v in enumerate(sys.argv):
        if v == 'fin':
            test( sys.argv[i+1] )
            doMain = False
    if doMain:
        main()
