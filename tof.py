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

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    header = '\033[95m'
    okblue = '\033[94m'
    okgreen = '\033[92m'
    warning = '\033[93m'
    fail = '\033[91m'
    endc = '\033[0m'
    bold = '\033[1m'
    underline = '\033[4m'

PRINTLEVEL = 0

def saveinfo(s):
    f = open('tof_save.txt', 'a')
    f.write(s + '\n')
    f.close()

def printgreen(s):
    print bcolors.okgreen + s + bcolors.endc

clight = 299792458.#  m / s

PDGID = {'e+' : -11, -11 : 'e+',
         'e-' : 11, 11 : 'e-',
         'nu_e' : 12, 12 : 'nu_e',
         'gamma' : 22, 22 : 'gamma',
         'mu+' : -13, -13 : 'mu+',
         'mu-' : -13, -13 : 'mu-',
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
        -13 : 105.6583745, 13 : 105.6583745,
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


def dot(a, b):
    if len(a) != len(b):
        print '!'*80
        print "tof.py : dot : Vectors' lengths don't match."
        print '!'*80
        raise SystemExit
    return sum( [x*b[i] for i, x in enumerate(a)] )


def cross_two_circles(r1, r2, x2, y2):
    '''
    The first cricle with a radius r1 is placed in the coordinate centre.
    While the second circle with a radius r2 is centred at the point (x2, y2).
    '''
    a = -2 * x2
    b = -2 * y2
    c = x2**2 + y2**2 + r1**2 - r2**2
    if PRINTLEVEL > 0:
        print 'Cross two circles'
        print '\tInput parameters'
        print '\t\tr1', r1
        print '\t\tr2', r2
        print '\t\tx2', x2
        print '\t\ty2', y2
        print '\tIntermediate parameters'
        print '\t\ta', a
        print '\t\tb', b
        print '\t\tc', c
        saveinfo('cross_two_circles.r1 {}'.format(r1))
        saveinfo('cross_two_circles.r2 {}'.format(r2))
        saveinfo('cross_two_circles.x2 {}'.format(x2))
        saveinfo('cross_two_circles.y2 {}'.format(y2))
        saveinfo('cross_two_circles.a {}'.format(a))
        saveinfo('cross_two_circles.b {}'.format(b))
        saveinfo('cross_two_circles.c {}'.format(c))
    if a == 0 and b == 0:
        res = []
    else:
        res = cross_circle_and_line(r1, a, b, c)
    return res

def cross_circle_and_line(r1, a, b, c):
    '''   a * x + b * y + c = 0   '''
    eps = 1.e-6
    x0 = -a*c/(a**2+b**2)
    y0 = -b*c/(a**2+b**2)
    if (c**2) > (r1**2*(a**2+b**2) + eps):
        return []
    elif abs(c**2 - r1**2*(a**2+b**2)) < eps:
        return [[x0, y0]]
    else:
        d = r1**2 - c**2/(a**2+b**2)
        mult = (d/(a**2+b**2))**.5
        ax = x0 + b*mult
        bx = x0 - b*mult
        ay = y0 - a*mult
        by = y0 + a*mult
        return [[ax, ay], [bx, by]]

def closest_phi(phi1, phi2, phi, wd):
    if phi1 < phi2:
        inde = [1, 2]
        a, b = phi1, phi2
    else:
        inde = [2, 1]
        a, b = phi2, phi1
    if a <= phi and phi <= b:
        if wd > 0:
            return [inde[1], b]
        else:
            return [inde[0], a]
    if phi < a:
        if wd > 0:
            return [inde[0], a]
        else:
            return [inde[1], b-2*math.pi]
    if phi > b:
        if wd > 0:
            return [inde[0], a+2*math.pi]
        else:
            return [inde[1], b]


class Tof:

    def __init__(self, r_in, r_out, z, sigma, eff):
        self.r_in = r_in
        self.r_out = r_out
        self.z = z
        self.th_0 = [abs(math.atan(r_in/z))]
        self.th_0.append( math.pi-self.th_0[0] )
        self.th_1 = [abs(math.atan(r_out/z))]
        self.th_1.append( math.pi-self.th_1[0] )
        self.sigma = sigma
        self.eff = eff
        if PRINTLEVEL > 0:
            print 'ToF init'
            print '\tr_in', self.r_in
            print '\tr_out', self.r_out
            print '\tz', self.z
            saveinfo('tof.r_in {}'.format(self.r_in))
            saveinfo('tof.r_out {}'.format(self.r_out))
            saveinfo('tof.z {}'.format(self.z))


    def get_time_v2_neutral(self, ptcl, field):
        if PRINTLEVEL > 0:
            print bcolors.header + 'The neutral way.' + bcolors.endc
        px0, py0 = ptcl.px, ptcl.py
        if ptcl.px != 0:
            a = - py0/px0
            b = 1.
            c = py0/px0*ptcl.vx - ptcl.vy
        else:#  The line is parralel to the vertical axis.
            a = 1.
            b = 0.
            c = - poi[0]
        cropoi_field = cross_circle_and_line(field.r, a, b, c)
        if len(cropoi_field) == 2:
            poi1 = cropoi_field[0]
            if ((poi1[0]-ptcl.vx)*px0 + (poi1[1]-ptcl.vy)*py0) < 0:
                poi1 = cropoi_field[1]
        cropoi_tof = cross_circle_and_line(self.r_out, a, b, c)
        if len(cropoi_tof) == 2:
            poi2 = cropoi_tof[0]
            if ((poi2[0]-ptcl.vx)*px0 + (poi2[1]-ptcl.vy)*py0) < 0:
                poi2 = cropoi_tof[1]
        #
        time = [ ((poi1[0]-ptcl.vx)**2 + (poi1[1]-ptcl.vy)**2)**.5 / ptcl.vt ]
        time.append( ((poi1[0]-poi2[0])**2 + (poi1[1]-poi2[1])**2)**.5 / ptcl.vt )
        #
        velz = ptcl.v * math.cos(ptcl.theta)
        rho0 = (ptcl.vx**2 + ptcl.vy**2)**.5
        z1 = ptcl.vz + velz * time[0]
        if abs(z1) > field.z:
            z1 = field.z * np.sign(velz)
            time[0] = ( (z1 - ptcl.vz) / velz )
        z2 = z1 + velz * time[1]
        if abs(z1 - z2) > 1.e-6:
            drho0to1 = abs( (z1-ptcl.vz) * math.tan( ptcl.theta ) )
            rho1 = rho0 + drho0to1
            poi1 = [ptcl.vx+drho0to1*math.cos(ptcl.phi), ptcl.vy+drho0to1*math.sin(ptcl.phi)]
            if abs(z2) < self.z:
                rho2 = z2 * math.tan( ptcl.theta )
            else:
                z2 = self.z * np.sign(velz)
                drho1to2 = abs( (z2-z1) * math.tan( ptcl.theta ) )
                rho2 = rho1 + drho1to2
                time[1] = (z2-z1)/abs(velz)
                if rho2 < self.r_in:
                    rho2 = rho2 + .1 * math.sin(ptcl.theta)
                    z2 = z2 + .1 * math.cos(ptcl.theta)
                    time[1] = -666.
                poi2[0] = poi1[0] + (rho2-rho1)*math.cos(ptcl.phi)
                poi2[1] = poi1[1] + (rho2-rho1)*math.sin(ptcl.phi)
        else:
            rho1 = field.r
            rho2 = self.r_out
        #
        px1 = px2 = ptcl.px
        py1 = py2 = ptcl.py
        pz1 = pz2 = ptcl.pz
        velx1 = velx2 = ptcl.v * math.sin( ptcl.theta ) * math.cos( ptcl.phi )
        vely1 = vely2 = ptcl.v * math.sin( ptcl.theta ) * math.sin( ptcl.phi )
        #
        if PRINTLEVEL > 0:
            saveinfo( 'tof.get_time_v2.mode neutral' )
            saveinfo( 'tof.get_time_v2.poi {} {}'.format(poi1[0], poi1[1]) )
            saveinfo( 'tof.get_time_v2.cropoi {} {} {} {}'.format(cropoi_tof[0][0],
                cropoi_tof[0][1], cropoi_tof[1][0], cropoi_tof[1][1]) )
            saveinfo( 'tof.get_time_v2.cropoi_true {} {}'.format(poi2[0], poi2[1]) )
            saveinfo( 'tof.get_time_v2.z1 {}'.format(z1) )
            saveinfo( 'tof.get_time_v2.z2 {}'.format(z2) )
            saveinfo( 'tof.get_time_v2.rho1 {}'.format(rho1) )
            saveinfo( 'tof.get_time_v2.rho2 {}'.format(rho2) )
            saveinfo( 'tof.get_time_v2.time {} {}'.format(time[0], time[1]) )
            saveinfo( 'tof.get_time_v2.px1 {}'.format(px1) )
            saveinfo( 'tof.get_time_v2.py1 {}'.format(py1) )
            saveinfo( 'tof.get_time_v2.pz1 {}'.format(pz1) )
            saveinfo( 'tof.get_time_v2.px2 {}'.format(px2) )
            saveinfo( 'tof.get_time_v2.py2 {}'.format(py2) )
            saveinfo( 'tof.get_time_v2.pz2 {}'.format(pz2) )
            saveinfo( 'tof.get_time_v2.velx1 {}'.format(velx1) )
            saveinfo( 'tof.get_time_v2.vely1 {}'.format(vely1) )
            saveinfo( 'tof.get_time_v2.velz1 {}'.format(velz) )
            saveinfo( 'tof.get_time_v2.velx2 {}'.format(velx2) )
            saveinfo( 'tof.get_time_v2.vely2 {}'.format(vely2) )
            saveinfo( 'tof.get_time_v2.velz2 {}'.format(velz) )
        res = { 'x1' : poi1[0], 'y1' : poi1[1], 'z1' : z1,
                'velx1' : velx1, 'vely1' : vely1, 'velz1' : velz,
                'x2' : poi2[0], 'y2' : poi2[1], 'z2' : z2,
                'velx2' : velx2, 'vely2' : vely2, 'velz2' : velz,
                'time0to1' : time[0], 'time1to2' : time[1] }
        return res


    def tracking(self, ptcl, field, x1, y1, z1, z2, x_lar, y_lar, r_lar, lar_phi0, wd, time):
        eps = 1.e-5
        print '\ttracking'
        omega = wd * abs( ptcl.vt / r_lar )
        print '\t\tomega = {}'.format( omega )
        if time[0] >= 0:
            dphi0to1 = omega * time[0]
        else:
            print '!' * 80
            raise SystemExit
        # saveinfo( 'tof.get_time_v2.velx0 {}'.format( ptcl.vt * math.cos(ptcl.phi) ) )
        # saveinfo( 'tof.get_time_v2.vely0 {}'.format( ptcl.vt * math.sin(ptcl.phi) ) )
        saveinfo( 'tof.get_time_v2.velx0 {}'.format( ptcl.v * math.sin(ptcl.theta) * math.cos(ptcl.phi) ) )
        saveinfo( 'tof.get_time_v2.vely0 {}'.format( ptcl.v * math.sin(ptcl.theta) * math.sin(ptcl.phi) ) )
        saveinfo( 'tof.get_time_v2.velz0 {}'.format( ptcl.v * math.cos(ptcl.theta) ) )
        # velx1 = ptcl.vt * math.cos(ptcl.phi + dphi0to1)
        # vely1 = ptcl.vt * math.sin(ptcl.phi + dphi0to1)
        velx1 = ptcl.v * math.sin(ptcl.theta) * math.cos(ptcl.phi + dphi0to1)
        vely1 = ptcl.v * math.sin(ptcl.theta) * math.sin(ptcl.phi + dphi0to1)
        velz =  ptcl.v * math.cos(ptcl.theta)
        x0to2 = [ptcl.vx]
        y0to2 = [ptcl.vy]
        z0to2 = [ptcl.vz]
        rho0to2 = [(ptcl.vx**2 + ptcl.vy**2)**.5]
        # r_lar, h_lar, x_lar, y_lar, wd, lar_phi0
        x, y = ptcl.vx, ptcl.vy
        t, dt = 0., .25e-9
        if time[0] < 5.e-9:
            dt = time[0] / 10.
        if time[0] / dt > 499:
            dt = time[0] / 499.
        # print '\t\t| z0to2[-1]={} |  <  | z1={} |'.format(z0to2[-1], z1)
        print '\t\t| z0to2[-1]={}  -  z1={} |  <>  eps={}'.format(z0to2[-1], z1, eps)
        # while abs(z0to2[-1]) < abs(z1):
        count = 0
        dzmin = abs(z0to2[-1]-z1)
        ddzmin = 0.
        while abs(z0to2[-1]-z1) > eps:
            t += dt
            x = x_lar + r_lar * math.cos(lar_phi0 + omega*t)
            y = y_lar + r_lar * math.sin(lar_phi0 + omega*t)
            z0to2.append( z0to2[-1] + velz*dt )
            rho0to2.append( (x**2 + y**2)**.5 )
            x0to2.append( x )
            y0to2.append( y )
            print '\t\t{:>4})  x = {:.2f}  y = {:.2f}  z = {:.2f}    omega * t  = {:.5g} * {:.5g} = {:.5g}'.format(count, x, y, z0to2[-1], omega, t, omega*t)
            #
            ddzmin = abs(z0to2[-1] - z1) - dzmin
            dzmin = abs(z0to2[-1] - z1)
            if ddzmin > 0.:
                printgreen( '\t\tddzmin has started increasing.' )
                break
            #
            count += 1
            if count > 500:
                break
        z0to2 = z0to2[:-1]
        z0to2.append(z1)
        rho0to2 = rho0to2[:-1]
        rho0to2.append( (x1**2 + y1**2)**.5 )
        t = 0.
        if abs(time[1]) < 5.e-9:
            dt = time[1] / 10.
        # while abs(z0to2[-1]) < abs(z2):
        count = 0
        dzmin = abs(z0to2[-1]-z2)
        ddzmin = 0.
        print '\t\t| z0to2[-1]={}  -  z2={} |  <>  eps={}'.format(z0to2[-1], z2, eps)
        while abs(z0to2[-1] - z2) > eps:
            t += dt
            x = x1 + velx1 * t
            y = y1 + vely1 * t
            rho0to2.append( (x**2 + y**2)**.5 )
            z0to2.append( z0to2[-1] + velz*dt )
            x0to2.append( x )
            y0to2.append( y )
            print '\t\t{:>4})  x = {:.2f}  y = {:.2f}  z = {:.2f}   dzmin = {:.2f}'.format(count, x, y, z0to2[-1], dzmin)
            #
            ddzmin = abs(z0to2[-1] - z2) - dzmin
            dzmin = abs(z0to2[-1] - z2)
            if ddzmin > 0.:
                printgreen( '\t\tddzmin has started increasing.  dzmin = {:.2f}'.format( dzmin ) )
                break
            #
            count += 1
            if count > 100:
                break
        saveinfo( 'tof.get_time_v2.x0to2 {}'.format(' '.join([str(x) for x in x0to2])) )
        saveinfo( 'tof.get_time_v2.y0to2 {}'.format(' '.join([str(x) for x in y0to2])) )
        saveinfo( 'tof.get_time_v2.z0to2 {}'.format(' '.join([str(x) for x in z0to2])) )
        saveinfo( 'tof.get_time_v2.rho0to2 {}'.format(' '.join([str(x) for x in rho0to2])) )


    def get_time_v2(self, ptcl, field):
        if ptcl.charge == 0.:
            return self.get_time_v2_neutral(ptcl, field)
        res = ptcl.get_Lar( field )
        r_lar, h_lar, x_lar, y_lar, wd, lar_phi0 = res[0], res[1], res[2][0], res[2][1], res[3], res[4]
        omega = wd * abs( ptcl.vt / r_lar )
        if self.r_out < field.r:
            # Effectively set the solenoid radius to the ToF one - 1 um
            rfield = self.r_out - 1.e-6
            zfield = self.z - 1.e-6
        else:
            rfield = field.r
            zfield = field.z
        cross_points = cross_two_circles(rfield, r_lar, x_lar, y_lar)
        if len(cross_points) == 2:
            poi1 = cross_points[0]
            poi2 = cross_points[1]
            phi1 = math.atan2(poi1[1]-y_lar, poi1[0]-x_lar) + math.pi
            phi2 = math.atan2(poi2[1]-y_lar, poi2[0]-x_lar) + math.pi
            phi = math.atan2(ptcl.vy-y_lar, ptcl.vx-x_lar) + math.pi
            clo_phi = closest_phi(phi1, phi2, phi, wd)
            if clo_phi[0] == 1:
                poi = poi1
                dphi = clo_phi[1] - phi
            else:
                poi = poi2
                dphi = clo_phi[1] - phi
            if PRINTLEVEL > 0:
                print 'get_time_v2'
                print '\tsolenoid leave point  poi ', poi
            time = [ abs(dphi * r_lar / ptcl.vt)]#  Time to reach the solenoid cylinder.
            #
            px = ptcl.pt * math.cos( ptcl.phi + dphi )
            py = ptcl.pt * math.sin( ptcl.phi + dphi )
            if PRINTLEVEL > 0:
                saveinfo( 'tof.get_time_v2.poi_prelim {} {}'.format(poi[0], poi[1]) )
                saveinfo( 'tof.get_time_v2.px1_prelim {}'.format(px) )
                saveinfo( 'tof.get_time_v2.py1_prelim {}'.format(py) )
            if px != 0:
                a = - py/px
                b = 1.
                c = py/px*poi[0] - poi[1]
            else:#  The line is parralel to the vertical axis.
                a = 1.
                b = 0.
                c = - poi[0]
            cropoi = cross_circle_and_line(self.r_out, a, b, c)
            if len(cropoi) == 2:
                poi1 = cropoi[0]
                poi2 = cropoi[1]
                d1 = ((poi1[0]-poi[0])**2 + (poi1[1]-poi[1])**2)**.5
                d2 = ((poi2[0]-poi[0])**2 + (poi2[1]-poi[1])**2)**.5
                #
                velx1 = ptcl.vt * math.cos( ptcl.phi + dphi )
                vely1 = ptcl.vt * math.sin( ptcl.phi + dphi )
                #
                if PRINTLEVEL > 0:
                    print '\tToF cross points  cropoi ', cropoi
                    saveinfo('tof.get_time_v2.cropoi {} {} {} {}'.format(poi1[0], poi1[1], poi2[0], poi2[1]))
                if d1 < d2:
                    time.append( d1 / ptcl.vt )
                    velrho1 = (velx1*poi1[0] + vely1*poi1[1]) / d1
                    x2, y2 = poi1[0], poi1[1]
                    # time.append( (self.r_out - field.r) / velrho1 )
                    if PRINTLEVEL > 0:
                        print '\t\tvelx1 {:.2g}   vely1 {:.2g}   velrho1 {:.2g}'.format( velx1, vely1, velrho1 )
                        saveinfo( 'tof.get_time_v2.cropoi_true {} {}'.format(poi1[0], poi1[1]) )
                else:
                    time.append( d2 / ptcl.vt )
                    velrho1 = (velx1*poi2[0] + vely1*poi2[1]) / d2
                    x2, y2 = poi2[0], poi2[1]
                    # time.append( (self.r_out - field.r) / velrho1 )
                    if PRINTLEVEL > 0:
                        print '\t\tvelx1 {:.2g}   vely1 {:.2g}   velrho1 {:.2g}'.format( velx1, vely1, velrho1 )
                        saveinfo( 'tof.get_time_v2.cropoi_true {} {}'.format(poi2[0], poi2[1]) )
            elif len(cropoi) == 1:
                poi1 = cropoi[0]
                d1 = (poi1[0]-poi[0])**2 + (poi1[1]-poi[1])**2
                if PRINTLEVEL > 0:
                    saveinfo( 'tof.get_time_v2.cropoi_true {} {}'.format(poi1[0], poi1[1]) )
                time.append( d1 / ptcl.vt )
            else:
                time.append( 0. )
            #  Now time needed to reach the ToF cylinder is added.
            #
            rho0 = (ptcl.vx**2 + ptcl.vy**2)**.5
            z0 = ptcl.vz
            velz =  ptcl.v * math.cos(ptcl.theta)
            if velz * ptcl.pz < 0:
                print "{}Momentum's and velocity's z direction aren't the same.{}".format(bcolors.warning, ptcl.pz, velz, bcolors.endc)
            z1 = z0 + velz * time[0]
            if abs(z1) <= zfield:#  Hits the solenoid cylinder.
                if PRINTLEVEL > 0:
                    printgreen('\tHits the solenoid cylinder:   abs(z1) = abs({:.3f})  <  zfield = {:.3f}'.format(z1, zfield))
                z2 = z1 + velz * time[1]
                if abs(z2) <= self.z:#  Hits the ToF cylinder.
                    velx2 = ptcl.v * math.sin( ptcl.theta ) * math.cos( ptcl.phi + omega * sum(time) )
                    vely2 = ptcl.v * math.sin( ptcl.theta ) * math.cos( ptcl.phi + omega * sum(time) )
                    if PRINTLEVEL > 0:
                        saveinfo( 'tof.get_time_v2.z1 {}'.format(z1) )
                        saveinfo( 'tof.get_time_v2.z2 {}'.format(z2) )
                        saveinfo( 'tof.get_time_v2.poi {} {}'.format(poi[0], poi[1]))
                        saveinfo( 'tof.get_time_v2.time {} {}'.format(time[0], time[1]) )
                        self.tracking(ptcl, field, poi[0], poi[1], z1, z2, x_lar, y_lar, r_lar, lar_phi0, wd, time)
                        saveinfo( 'tof.get_time_v2.case Hits_the_solenoid_and_ToF_cylinders.' )
                    res = { 'x1' : poi1[0], 'y1' : poi1[1], 'z1' : z1,
                            'velx1' : velx1, 'vely1' : vely1, 'velz1' : velz,
                            'x2' : x2, 'y2' : y2, 'z2' : z2,
                            'velx2' : velx2, 'vely2' : vely2, 'velz2' : velz,
                            'time0to1' : time[0], 'time1to2' : time[1] }
                    return res
                else:#  Towards a ToF cap.
                    time[1] = (self.z - abs(z1)) / abs(velz)
                    # rhit = rfield + ptcl.vt * time[1]
                    rhit = rfield + (self.z - abs(z1)) * math.tan(ptcl.theta)
                    velx2 = ptcl.v * math.sin( ptcl.theta ) * math.cos( ptcl.phi + omega * sum(time) )
                    vely2 = ptcl.v * math.sin( ptcl.theta ) * math.cos( ptcl.phi + omega * sum(time) )
                    if PRINTLEVEL > 0:
                        saveinfo( 'tof.get_time_v2.z1 {}'.format(z1) )
                        saveinfo('tof.get_time_v2.poi {} {}'.format(poi[0], poi[1]))
                    if rhit > self.r_in:#  Hits a ToF cap.
                        if PRINTLEVEL > 0:
                            saveinfo( 'tof.get_time_v2.z2 {}'.format(self.z) )
                            saveinfo('tof.get_time_v2.poi {} {}'.format(poi[0], poi[1]))
                            saveinfo( 'tof.get_time_v2.time {} {}'.format(time[0], time[1]) )
                            self.tracking(ptcl, field, poi[0], poi[1], z1, z2, x_lar, y_lar, r_lar, lar_phi0, wd, time)
                        # return sum(time)
                    else:#  Goes through the hole.
                        if PRINTLEVEL > 0:
                            saveinfo( 'tof.get_time_v2.z2 {}'.format(self.z*2) )
                            saveinfo('tof.get_time_v2.poi {} {}'.format(poi[0], poi[1]))
                            saveinfo( 'tof.get_time_v2.time {} {}'.format(time[0], -666.) )
                            self.tracking(ptcl, field, poi[0], poi[1], z1, z2, x_lar, y_lar, r_lar, lar_phi0, wd, time)
                        # return -666.
                        time[1] = -666.
                    res = { 'x1' : poi1[0], 'y1' : poi1[1], 'z1' : z1,
                            'velx1' : velx1, 'vely1' : vely1, 'velz1' : velz,
                            'x2' : x2, 'y2' : y2, 'z2' : z2,
                            'velx2' : velx2, 'vely2' : vely2, 'velz2' : velz,
                            'time0to1' : time[0], 'time1to2' : time[1] }
                    return res
            else:#  Outs from the solenoid through the cap.
                if PRINTLEVEL > 0:
                    printgreen('\tOuts from the solenoid through the cap:   abs(z1) = abs({:.3f})  >  zfield = {:.3f}'.format(z1, zfield))
                t0new = (np.sign(velz)*zfield - z0) / velz
                if t0new < 0:
                    print bcolors.warning + 't0new = {} is less than zero!'.format(t0new) + bcolors.endc
                t1new = (self.z - zfield) / abs(velz)
                if t1new < 0:
                    print bcolors.warning + 't1new = {} is less than zero!'.format(t1new) + bcolors.endc
                dphinew = dphi*t0new/time[0]
                if PRINTLEVEL > 0:
                    print '\tNew dphi angle :  {:.3g}  -->  {:.3g} = {:.3g} * {:.3g} / {:.3g}'.format(dphi, dphinew, dphi, t0new, time[0])
                    print '\tNew dphi angle :  {:.3g}  -->  {:.3g} = {:.3g} * {:.3g} / {:.3g}'.format(dphi*180/math.pi,
                            dphinew*180/math.pi, dphi*180/math.pi, t0new, time[0])
                # x1 = ptcl.vx + x_lar + r_lar * math.cos(phi+dphinew)
                # y1 = ptcl.vx + y_lar + r_lar * math.sin(phi+dphinew)
                x1 = x_lar + r_lar * math.cos(lar_phi0+dphinew)
                y1 = y_lar + r_lar * math.sin(lar_phi0+dphinew)
                if PRINTLEVEL > 0:
                    saveinfo('tof.get_time_v2.poi {} {}'.format(x1, y1))
                    # trk_x_tmp = [ptcl.vx + x_lar + r_lar * math.cos(lar_phi0+.1*x) for x in xrange(1, 25)]
                    # trk_y_tmp = [ptcl.vy + y_lar + r_lar * math.sin(lar_phi0+.1*x) for x in xrange(1, 25)]
                    trk_x_tmp = [x_lar + r_lar * math.cos(lar_phi0+.1*x*wd) for x in xrange(1, 25)]
                    trk_y_tmp = [y_lar + r_lar * math.sin(lar_phi0+.1*x*wd) for x in xrange(1, 25)]
                    saveinfo('tof.get_time_v2.trk_x {}'.format(' '.join([str(x) for x in trk_x_tmp])))
                    saveinfo('tof.get_time_v2.trk_y {}'.format(' '.join([str(x) for x in trk_y_tmp])))
                velx1 = ptcl.vt * math.cos( ptcl.phi + dphi*t0new/time[0] )
                vely1 = ptcl.vt * math.sin( ptcl.phi + dphi*t0new/time[0] )
                rho1 = rho0 + (zfield-z0) * math.tan(ptcl.theta)
                velrho1 = (velx1*x1 + vely1*y1) / (x1**2 + y1**2)**.5
                # rho2 = rho1 + velrho1 * t1new
                x2 = x1 + velx1 * t1new
                y2 = y1 + vely1 * t1new
                rho2 = (x2**2 + y2**2)**.5
                # z1 = field.z * np.sign( math.cos(ptcl.theta) )
                z1 = field.z * np.sign( ptcl.pz )
                px1 = ptcl.pt * math.cos( ptcl.phi + dphi*t0new/time[0] )
                py1 = ptcl.pt * math.sin( ptcl.phi + dphi*t0new/time[0] )
                prho1 = (px1*x1 + py1*y1) / (x1**2 + y1**2)**.5
                time = [t0new, t1new]
                if PRINTLEVEL > 0:
                    saveinfo( 'tof.get_time_v2.z1 {}'.format(z1) )
                    saveinfo( 'tof.get_time_v2.velrho1 {}'.format(velrho1) )
                    saveinfo( 'tof.get_time_v2.prho1 {}'.format(prho1) )
                    saveinfo( 'tof.get_time_v2.px1 {}'.format(px1) )
                    saveinfo( 'tof.get_time_v2.py1 {}'.format(py1) )
                    saveinfo( 'tof.get_time_v2.rho2 {}'.format(rho2) )
                    saveinfo( 'tof.get_time_v2.x2 {}'.format(x2) )
                    saveinfo( 'tof.get_time_v2.y2 {}'.format(y2) )
                if self.r_in <= rho2 and rho2 <= self.r_out:#  Hits a ToF cap.
                    z2 = self.z*np.sign(ptcl.pz)
                    if PRINTLEVEL > 0:
                        printgreen('\tHits a ToF cap:   r_in = {:.3f}  <  rho2 = {:.3f}  <  r_out = {:.3f}'.format(self.r_in, rho2, self.r_out))
                        saveinfo( 'tof.get_time_v2.z2 {}'.format(self.z*np.sign(ptcl.pz)) )
                        saveinfo( 'tof.get_time_v2.time {} {}'.format(t0new, t1new) )
                        self.tracking(ptcl, field, x1, y1, z1, z2, x_lar, y_lar, r_lar, lar_phi0, wd, time)
                        saveinfo( 'tof.get_time_v2.case Outs_from_the_solenoid_through_the_cap_and_hit_a_ToF_cap.' )
                    # return t0new + t1new
                    time = [t0new, t1new]
                elif rho2 > self.r_out:#  Hits the ToF cylinder.
                    if PRINTLEVEL > 0:
                        printgreen('\tHits the ToF cylinder:   rho2 = {:.3f}  >  r_out = {:.3f}'.format(rho2, self.r_out))
                    t1newnew = (self.r_out - rho1) / velrho1
                    z2 = z1 + (self.r_out - rho1) / math.tan( ptcl.theta )
                    if PRINTLEVEL > 0:
                        saveinfo( 'tof.get_time_v2.z2 {}'.format(z2) )
                        saveinfo( 'tof.get_time_v2.time {} {}'.format(t0new, t1newnew) )
                        self.tracking(ptcl, field, x1, y1, z1, z2, x_lar, y_lar, r_lar, lar_phi0, wd, time)
                        saveinfo( 'tof.get_time_v2.case Outs_from_the_solenoid_through_the_cap_and_hit_the_ToF_cylinder.' )
                    # return t0new + t1newnew
                    time = [t0new, t1newnew]
                else:#  Goes through the ToF hole.
                    z2 = self.z*np.sign(ptcl.pz)
                    if PRINTLEVEL > 0:
                        printgreen('\tGoes through the ToF hole:   rho2 = {:.3f}  <  r_in = {:.3f}'.format(rho2, self.r_in))
                        saveinfo( 'tof.get_time_v2.z2 {}'.format((self.z+.25)*np.sign(ptcl.pz)) )
                        saveinfo( 'tof.get_time_v2.time {} {}'.format(t0new, -666.) )
                        self.tracking(ptcl, field, x1, y1, z1, z2, x_lar, y_lar, r_lar, lar_phi0, wd, time)
                    # return -666.
                    time = [t0new, -666.]
                velx2 = ptcl.v * math.sin( ptcl.theta ) * math.cos( ptcl.phi + omega * sum(time) )
                vely2 = ptcl.v * math.sin( ptcl.theta ) * math.cos( ptcl.phi + omega * sum(time) )
                res = { 'x1' : poi1[0], 'y1' : poi1[1], 'z1' : z1,
                        'velx1' : velx1, 'vely1' : vely1, 'velz1' : velz,
                        'x2' : x2, 'y2' : y2, 'z2' : z2,
                        'velx2' : velx2, 'vely2' : vely2, 'velz2' : velz,
                        'time0to1' : time[0], 'time1to2' : time[1] }
                return res
        elif len(cross_points) == 1:
            pass
        else:
            if PRINTLEVEL > 0:
                print bcolors.header + '\nThe helix way.' + bcolors.endc
            velz = ptcl.v * math.cos( ptcl.theta )
            if velz != 0:
                z1 = np.sign(velz) * field.z
                # time = [ (z1 - ptcl.vz) / abs(velz) ]
                time = [ (z1 - ptcl.vz) / velz ]
                z2 = np.sign(velz) * self.z
                # time.append( (z2 - z1) / abs(velz) )
                time.append( (z2 - z1) / velz )
                omega = wd * ptcl.vt / r_lar
                dphi0to1 = omega * time[0]
                # dphi1to2 = omega * time[1]
                x1 = x_lar + r_lar * math.cos(lar_phi0 + dphi0to1)
                y1 = y_lar + r_lar * math.sin(lar_phi0 + dphi0to1)
                velx1 = ptcl.vt * math.cos(ptcl.phi + dphi0to1)
                vely1 = ptcl.vt * math.sin(ptcl.phi + dphi0to1)
                velrho1 = (velx1*x1 + vely1*y1) / (x1**2 + y1**2)**.5
                rho1 = (x1**2 + y1**2)**.5
                x2 = x1 + velx1 * time[1]
                y2 = y1 + vely1 * time[1]
                rho2 = (x2**2 + y2**2)**.5
                velz1 = velz
                if rho2 < self.r_in:
                    time[1] = -666.
                if rho2 > self.r_out:
                    #
                    if velx1 != 0:
                        a = - vely1 / velx1
                        b = 1.
                        c = - a * x1 - y1
                    else:
                        a = 1.
                        b = 0.
                        c = - x1
                    poi2prelim = cross_circle_and_line(self.r_out, a, b, c)
                    if len(poi2prelim) == 2:
                        d1 = [poi2prelim[0][0]-x1, poi2prelim[0][1]-y1]
                        if dot( [velx1, vely1], d1 ) > 0:
                            x2 = poi2prelim[0][0]
                            y2 = poi2prelim[0][1]
                        else:
                            x2 = poi2prelim[1][0]
                            y2 = poi2prelim[1][1]
                    time[1] = ( (x2-x1)**2 + (y2-y1)**2 )**.5 / ptcl.vt
                    #
                    rho2 = (x2**2 + y2**2)**.5
                    z2 = z1 + velz * time[1]
                if PRINTLEVEL > 0:
                    print '\tomega', omega
                    print '\tr_lar', r_lar
                    print '\tlar_phi0', lar_phi0
                    print '\tdphi0to1', dphi0to1
                    # print '\tdphi1to2', dphi1to2
                    saveinfo( 'tof.get_time_v2.velx1 {}'.format(velx1) )
                    saveinfo( 'tof.get_time_v2.vely1 {}'.format(vely1) )
                    saveinfo( 'tof.get_time_v2.velz1 {}'.format(velz1) )
                    saveinfo( 'tof.get_time_v2.omega {}'.format(omega) )
                    saveinfo( 'tof.get_time_v2.poi {} {}'.format(x1, y1) )
                    saveinfo( 'tof.get_time_v2.cropoi_true {} {}'.format(x2, y2) )
                    saveinfo( 'tof.get_time_v2.z1 {}'.format(z1) )
                    saveinfo( 'tof.get_time_v2.z2 {}'.format(z2) )
                    saveinfo( 'tof.get_time_v2.time {} {}'.format(time[0], time[1]) )
                    self.tracking(ptcl, field, x1, y1, z1, z2, x_lar, y_lar, r_lar, lar_phi0, wd, time)
                # return sum(time)
                velx2 = ptcl.v * math.sin( ptcl.theta ) * math.cos( ptcl.phi + omega * sum(time) )
                vely2 = ptcl.v * math.sin( ptcl.theta ) * math.cos( ptcl.phi + omega * sum(time) )
                res = { 'x1' : x1, 'y1' : y1, 'z1' : z1,
                        'velx1' : velx1, 'vely1' : vely1, 'velz1' : velz,
                        'x2' : x2, 'y2' : y2, 'z2' : z2,
                        'velx2' : velx2, 'vely2' : vely2, 'velz2' : velz,
                        'time0to1' : time[0], 'time1to2' : time[1] }
                return res
            else:
                return -666.

    def get_time(self, ptcl, field):
        if self.th_1[0] < ptcl.theta and ptcl.theta < self.th_1[1]:
            if ptcl.charge != 0:
                r_lar = 3.3 * (ptcl.p2+ptcl.m**2)**.5 * ptcl.vt / ptcl.charge / field.mag
                if self.r_out < field.r:
                    if 2*r_lar > self.r_out:
                        phi_curve = 2.*math.asin(self.r_out/2./r_lar)
                        h_lar = r_lar * ptcl.pz / ptcl.pt
                        time = ((h_lar/2./math.pi)**2+r_lar**2)**.5 * phi_curve / ptcl.v
                    else:
                        time = field.z / (ptcl.v * abs(math.cos(ptcl.theta)) )
                else:
                    if 2*r_lar > field.r:
                        phi_curve = 2.*math.asin(field.r/2./r_lar)
                        h_lar = r_lar * ptcl.pz / ptcl.pt
                        time = ((h_lar/2./math.pi)**2+r_lar**2)**.5 * phi_curve / ptcl.v
                        if self.r_out > field.r:
                            time += (self.r_out-field.r) / (ptcl.vt * abs(math.cos(phi_curve/2)) )
                    else:
                        time = field.z / (ptcl.v * abs(math.cos(ptcl.theta)) )
            else:
                time = self.r_out / ptcl.vt
        elif self.th_0[0] < ptcl.theta and ptcl.theta < self.th_0[1]:
            time = self.z / (ptcl.v * abs(math.cos(ptcl.theta)) )
        else:
            time = -600.
        return time


class Field:

    def __init__(self, mag, r, z):
        self.mag = mag
        self.r = r
        self.z = z
        if PRINTLEVEL > 0:
            saveinfo( 'field.init.mag {}'.format(self.mag) )
            saveinfo( 'field.init.r {}'.format(self.r) )
            saveinfo( 'field.init.z {}'.format(self.z) )


class Particle:

    def __init__(self, pdgid, charge, m, px, py, pz, vx=0., vy=0., vz=0.):
        self.pdgid = pdgid
        self.charge = charge
        self.m = m
        self.px = px
        self.py = py
        self.pz = pz
        self.px2 = px**2
        self.py2 = py**2
        self.pz2 = pz**2
        self.pt2 = self.px2+self.py2
        self.pt = self.pt2**.5
        self.p2 = self.pt2+self.pz2
        self.p = self.p2**.5
        self.beta = (1+self.m**2/self.p2)**-.5
        self.v = self.beta * clight
        self.theta = math.acos(self.pz/self.p)
        self.phi = math.atan2(self.py, self.px)# + math.pi
        self.betat = self.beta * math.sin(self.theta)
        self.vt = abs( self.v * math.sin(self.theta) )
        # vertex
        self.vx = vx
        self.vy = vy
        self.vz = vz
        if PRINTLEVEL > 0:
            print 'Particle'
            print '\t{:>8}{:>8}{:>8}{:>8}{:>8}'.format('px', 'py', 'pz', 'pt', 'p')
            print '\t{:8.4f}{:8.4f}{:8.4f}{:8.4f}{:8.4f}'.format(self.px, self.py, self.pz, self.pt, self.p)
            print '\ttheta {}   phi {}'.format(self.theta, self.phi)
            saveinfo( 'particle.init.vertex {} {} {}'.format(self.vx, self.vy, self.vz) )
            saveinfo( 'particle.init.p {} {} {}'.format(self.px, self.py, self.pz) )
            saveinfo( 'particle.init.charge {}'.format(self.charge) )
            saveinfo( 'particle.init.mass {}'.format(self.m) )

    def get_Lar(self, field):
        bz = field.mag
        if self.charge == 0:
            if PRINTLEVEL > 0:
                print 'Larmor radius'
                print '\tThe neutral particle has been passed.'
            return 0., 0., [0., 0.], 0.
        if bz == 0:
            if PRINTLEVEL > 0:
                print 'Larmor radius'
                print '\tThe zero magnet field.'
            return 0., 0., [0., 0.], 0.
        # r_lar = 3.3 * (self.p2+self.m**2)**.5 * self.betat / self.charge / bz
        r_lar = 3.3 * self.pt / self.charge / bz
        r_lar = abs(r_lar)
        h_lar = r_lar * self.pz / self.pt
        # Find the helix centre.
        x = abs( r_lar * self.py / self.pt )
        y = abs( r_lar * self.px / self.pt )
        # x = r_lar * self.py / self.pt
        # y = r_lar * self.px / self.pt
        if x == 0:
            x = 1.e-6
        if y == 0:
            y = 1.e-6
        if self.px > 0 and self.py > 0:
            # xc = self.vx + x*np.sign(field.mag)*np.sign(self.charge)
            # yc = self.vy - y*np.sign(field.mag)*np.sign(self.charge)
            xc = self.vx - x*np.sign(field.mag)*np.sign(self.charge)
            yc = self.vy + y*np.sign(field.mag)*np.sign(self.charge)
            if PRINTLEVEL > 0:
                saveinfo('particle.get_lar.x.sign -1')
                saveinfo('particle.get_lar.y.sign 1')
        elif self.px > 0 and self.py < 0:
            # xc = self.vx - x*np.sign(field.mag)*np.sign(self.charge)
            # yc = self.vy - y*np.sign(field.mag)*np.sign(self.charge)
            xc = self.vx + x*np.sign(field.mag)*np.sign(self.charge)
            yc = self.vy + y*np.sign(field.mag)*np.sign(self.charge)
            if PRINTLEVEL > 0:
                saveinfo('particle.get_lar.x.sign 1')
                saveinfo('particle.get_lar.y.sign 1')
        elif self.px < 0 and self.py < 0:
            # xc = self.vx - x*np.sign(field.mag)*np.sign(self.charge)
            # yc = self.vy + y*np.sign(field.mag)*np.sign(self.charge)
            xc = self.vx + x*np.sign(field.mag)*np.sign(self.charge)
            yc = self.vy - y*np.sign(field.mag)*np.sign(self.charge)
            if PRINTLEVEL > 0:
                saveinfo('particle.get_lar.x.sign 1')
                saveinfo('particle.get_lar.y.sign -1')
        elif self.px < 0 and self.py > 0:
            # xc = self.vx + x*np.sign(field.mag)*np.sign(self.charge)
            # yc = self.vy + y*np.sign(field.mag)*np.sign(self.charge)
            xc = self.vx - x*np.sign(field.mag)*np.sign(self.charge)
            yc = self.vy - y*np.sign(field.mag)*np.sign(self.charge)
            if PRINTLEVEL > 0:
                saveinfo('particle.get_lar.x.sign -1')
                saveinfo('particle.get_lar.y.sign -1')
        elif self.px == 0 and self.py != 0:
            # xc = self.vx + r_lar * np.sign(field.mag)*np.sign(self.charge)*np.sign(self.py)
            xc = self.vx - r_lar * np.sign(field.mag)*np.sign(self.charge)*np.sign(self.py)
            yc = self.vy
            if PRINTLEVEL > 0:
                saveinfo('particle.get_lar.x.sign -1')
                saveinfo('particle.get_lar.y.sign 0')
        elif self.px != 0 and self.py == 0:
            xc = self.vx
            # yc = self.vy + r_lar * np.sign(field.mag)*np.sign(self.charge)*np.sign(self.px)
            yc = self.vy + r_lar * np.sign(field.mag)*np.sign(self.charge)*np.sign(self.px)
            if PRINTLEVEL > 0:
                saveinfo('particle.get_lar.x.sign 0')
                saveinfo('particle.get_lar.y.sign -1')
        lar_phi0 = math.atan2( self.vy-yc, self.vx-xc)
        #
        # Find the angula speed direction.
        # -1 --- an azimuthal angle decreases with time.
        #  1 --- an azimuthal angle increases with time.
        wd = np.sign(field.mag) * np.sign(self.charge)
        #
        # return r_lar, h_lar, [xc, yc], wd
        if PRINTLEVEL > 0:
            print 'Larmor radius'
            print '\tBz', bz
            print '\tpt', self.pt
            print '\tr_lar', r_lar
            print '\th_lar', h_lar
            print '\tx = {}    y = {}'.format(x, y)
            print '\tcentre', xc, yc
            print '\twd', wd
            saveinfo('particle.get_lar.r_lar {}'.format(r_lar))
            saveinfo('particle.get_lar.h_lar {}'.format(h_lar))
            saveinfo('particle.get_lar.centre {} {}'.format(xc, yc))
            saveinfo('particle.get_lar.lar_phi0 {}'.format(lar_phi0))
            saveinfo('particle.get_lar.x {}'.format(x))
            saveinfo('particle.get_lar.y {}'.format(y))
        return [r_lar, h_lar, [xc, yc], wd, lar_phi0]


def main():
    Finname = 'output.root'
    Fin = r.TFile(Finname)
    t = Fin.Get('events')
    nentries = t.GetEntries()
    print '{} entries in the tree.'.format(nentries)
    ################
    tof_r_out = 1.090# m,  Time-of-Flight radius
    tof_r_in = .250# m,  Time-of-Flight radius around the beam line
    tof_z = 2.550/2# m
    tof_th_0 = [abs(math.atan(tof_r_in/tof_z))]
    tof_th_0.append( math.pi-tof_th_0[0] )
    tof_th_1 = [abs(math.atan(tof_r_out/tof_z))]
    tof_th_1.append( math.pi-tof_th_1[0] )
    tof_sigma_avr = 1.3# ns,  STD of average time
    tof_sigma_0 = 1.7# ns,  STD of a single counter
    tof_eff = 0.99
    tof_code_along_pipe = -100.
    tof_code_lost_signal = -101.
    B_mag = 1.# T, magnetic field magnitude
    B_r = 1.575# m, magnetic field radius
    B_z = 3.780/2# m, magnetic field half length
    ################
    field = Field(1., 1.575, 3.780/2)
    tof = Tof(.250, 1.090, 2.550/2, .1, .99)
    mu0 = Tof(.950, 1.930, 3.900/2, .1, .99)
    ################
    htime = r.TH1D('htime', 'htime;TOF time, ns;Entries', 1000, 0, 100)
    hv = r.TH1D('hv', 'hv;velocity;Entries', 1000, 0, 1)
    htimep = r.TH2D('htimep', 'htimep;Momentum p, MeV/c;TOF time, ns;',
            100, 0, 2, 100, 0, 20)
    ################
    Fout = r.TFile('eos.root', 'recreate')
    tout = r.TTree('t', 't')
    npat = array.array('i', [0])
    toftime = array.array('d', [0.]*100)
    toftime0 = array.array('d', [0.]*100)
    mutime = array.array('d', [0.]*100)
    mutime0 = array.array('d', [0.]*100)
    charge = array.array('i', [0]*100)
    pdgid = array.array('i', [0]*100)
    th = array.array('d', [0.]*100)
    phi = array.array('d', [0.]*100)
    mom = array.array('d', [0.]*100)
    tout.Branch('n', npat, 'n/I')
    tout.Branch('charge', charge, 'charge[n]/I')
    tout.Branch('pdgid', pdgid, 'pdgid[n]/I')
    tout.Branch('toftime0', toftime0, 'toftime0[n]/D')
    tout.Branch('toftime', toftime, 'toftime[n]/D')
    tout.Branch('mutime0', mutime0, 'mutime0[n]/D')
    tout.Branch('mutime', mutime, 'mutime[n]/D')
    tout.Branch('p', mom, 'p[n]/D')
    tout.Branch('th', th, 'th[n]/D')
    tout.Branch('phi', phi, 'phi[n]/D')
    ################
    i = 0
    if 'full' in sys.argv:
        imax = -1
    else:
        imax = 50000
    for eve in t:
        if i % 10000 == 0:
            print  '{}   {:%Y-%m-%d %H:%M:%S}'.format( i, datetime.datetime.now() )
        if imax > 0 and i > imax:
            break
        gps = eve.GenParticle
        # print gps, 'len', len(gps)
        npat[0] = len(gps)
        for j in xrange(len(gps)):
            px = gps[j].core.p4.px
            py = gps[j].core.p4.py
            pz = gps[j].core.p4.pz
            mass = gps[j].core.p4.mass
            charge[j] = gps[j].core.charge
            pdgid[j] = gps[j].core.pdgId
            #
            ptcl = Particle(pdgid[j], charge[j], mass, px, py, pz)
            mom[j] = ptcl.p
            th[j] = ptcl.theta
            phi[j] = ptcl.phi
            #
            toftime0[j] = tof.get_time( ptcl, field )
            time = toftime0[j]
            mutime0[j] = mu0.get_time( ptcl, field )
            #
            if time > 0:
                dice = random.random()
                if dice > tof_eff:
                    time = tof_code_lost_signal
                else:
                    # toftime[j] = time+gauss(0., tof_sigma_0)
                    toftime[j] = time+gauss(0., tof_sigma_avr)
                    toftime0[j] = time
            htime.Fill(time)
            htimep.Fill(ptcl.p, time)
        tout.Fill()
        i += 1
    c1 = r.TCanvas('c1', 'c1', 600, 600)
    '''
    c1.Divide(2, 2)
    c1.cd(1)
    htime.Draw()
    c1.cd(2)
    hv.Draw()
    c1.cd(3)
    '''
    htimep.Draw('colz')
    c1.Update()
    raw_input()
    tout.Write()
    Fout.Close()


if __name__ == '__main__':
    main()
