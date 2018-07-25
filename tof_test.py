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


def test1():
    os.system('rm -rf tof_save.txt')
    field = mytof.Field(1., 1., 2.)
    # def __init__(self, pdgid, charge, m, px, py, pz, vx=0., vy=0., vz=0.):
    pat = mytof.Particle(13, -1., 0.1056583745, .5/2**.5, .5/2**.5, 0., 0., 0., 0.)
    # def __init__(self, r_in, r_out, z, sigma, eff):
    tof = mytof.Tof(.5, 2., 3., 25.e-9, .999)
    print tof.get_time_v2(pat, field)
    os.system('cp -f tof_save.txt tof_test_1.txt')

def test1b():
    os.system('rm -rf tof_save.txt')
    field = mytof.Field(1., 1., 2.)
    pat = mytof.Particle(-13, 1., 0.1056583745, .5/2**.5, .5/2**.5, 0., 0., 0., 0.)
    tof = mytof.Tof(.5, 2., 3., 25.e-9, .999)
    print tof.get_time_v2(pat, field)
    os.system('cp -f tof_save.txt tof_test_1a.txt')


def test2():
    os.system('rm -rf tof_save.txt')
    field = mytof.Field(1., 1., 2.)
    pat = mytof.Particle(13, -1., 0.1056583745, .5/2**.5, .5/2**.5, 0, 0.5, 0.25, 0.)
    tof = mytof.Tof(.5, 2., 3., 25.e-9, .999)
    print tof.get_time_v2(pat, field)
    os.system('cp -f tof_save.txt tof_test_2.txt')


def test3():
    os.system('rm -rf tof_save.txt')
    field = mytof.Field(1., 1., 2.)
    pat = mytof.Particle(13, -1., 0.1056583745, .5/3**.5, .5/3**.5, .5/3**.5, 0.5, 0.25, -.75)
    tof = mytof.Tof(.5, 2., 3., 25.e-9, .999)
    print tof.get_time_v2(pat, field)
    os.system('cp -f tof_save.txt tof_test_3.txt')


def test4():
    os.system('rm -rf tof_save.txt')
    field = mytof.Field(1., 1., 2.)
    pat = mytof.Particle(13, -1., 0.1056583745, .5/2**.5, .5/2**.5, 1., 0.5, 0.25, 1.75)
    tof = mytof.Tof(.5, 2., 3., 25.e-9, .999)
    print tof.get_time_v2(pat, field)
    os.system('cp -f tof_save.txt tof_test_4.txt')


def test5():
    os.system('rm -rf tof_save.txt')
    field = mytof.Field(1., 1., 2.)
    pat = mytof.Particle(13, -1., 0.1056583745, .5/2**.5, .5/2**.5, 1., 0.5, 0.25, 1.5)
    tof = mytof.Tof(1.5, 2., 3., 25.e-9, .999)
    print tof.get_time_v2(pat, field)
    os.system('cp -f tof_save.txt tof_test_5.txt')


def test6():
    os.system('rm -rf tof_save.txt')
    field = mytof.Field(1., 1., 2.)
    pat = mytof.Particle(-13, 1., 0.1056583745, .5/2**.5, .5/2**.5, 1., 0.5, 0.25, 1.5)
    tof = mytof.Tof(1.5, 2., 3., 25.e-9, .999)
    print tof.get_time_v2(pat, field)
    os.system('cp -f tof_save.txt tof_test_6.txt')

def test6b():
    os.system('rm -rf tof_save.txt')
    field = mytof.Field(1., 1., 2.)
    pat = mytof.Particle(-13, 1., 0.1056583745, .3/2**.5, .3/2**.5, -1., -.90/2**.5, .90/2**.5, 1.5)
    tof = mytof.Tof(1.5, 2., 3., 25.e-9, .999)
    print tof.get_time_v2(pat, field)
    os.system('cp -f tof_save.txt tof_test_6b.txt')

def test6c():
    os.system('rm -rf tof_save.txt')
    field = mytof.Field(1., 1., 2.)
    phi = math.pi / 180 * 88.
    pat = mytof.Particle(-13, 1., 0.1056583745, .3*math.cos(phi), .3*math.sin(phi), -.1, -.99, 0.001, 1.5)
    tof = mytof.Tof(1.5, 2., 3., 25.e-9, .999)
    print tof.get_time_v2(pat, field)
    os.system('cp -f tof_save.txt tof_test_6c.txt')


def test7():
    os.system('rm -rf tof_save.txt')
    field = mytof.Field(1., 1., 2.)
    pat = mytof.Particle(130, 0., pdg2m(130), .5/2**.5, .5/2**.5, 0., 0.5, 0.25, 1.5)
    tof = mytof.Tof(1.5, 2., 3., 25.e-9, .999)
    print tof.get_time_v2(pat, field)
    os.system('cp -f tof_save.txt tof_test_7.txt')


def test8():
    os.system('rm -rf tof_save.txt')
    field = mytof.Field(1., 1., 2.)
    pat = mytof.Particle(130, 0., pdg2m(130), .5/2**.5, .5/2**.5, 1., 0.5, 0.25, 1.5)
    tof = mytof.Tof(1.5, 2., 3., 25.e-9, .999)
    print tof.get_time_v2(pat, field)
    os.system('cp -f tof_save.txt tof_test_8.txt')

def test8a():
    os.system('rm -rf tof_save.txt')
    field = mytof.Field(1., 1., 2.)
    pat = mytof.Particle(130, 0., pdg2m(130), .5/2**.5, .5/2**.5, 1., 0.5, 0.25, -1.5)
    tof = mytof.Tof(1.5, 2., 3., 25.e-9, .999)
    print tof.get_time_v2(pat, field)
    os.system('cp -f tof_save.txt tof_test_8a.txt')

def test8b():
    os.system('rm -rf tof_save.txt')
    field = mytof.Field(1., 1., 2.)
    pat = mytof.Particle(130, 0., pdg2m(130), .5/2**.5, .5/2**.5, 1., 0.5, 0.25, 1.5)
    tof = mytof.Tof(1.5, 2., 3., 25.e-9, .999)
    print tof.get_time_v2(pat, field)
    os.system('cp -f tof_save.txt tof_test_8.txt')


def test9():
    os.system('rm -rf tof_save.txt')
    field = mytof.Field(1., 1., 2.)
    pat = mytof.Particle(-13, 1., 0.1056583745, .10/2**.5, .10/2**.5, .1, 0.1, 0.05, -1.5)
    tof = mytof.Tof(1.5, 2., 3., 25.e-9, .999)
    print tof.get_time_v2(pat, field)
    os.system('cp -f tof_save.txt tof_test_9.txt')

def test9b():
    os.system('rm -rf tof_save.txt')
    field = mytof.Field(1., 1., 2.)
    pat = mytof.Particle(-13, 1., 0.1056583745, .10/2**.5, .10/2**.5, .1, 0.0, 0.00, -1.5)
    tof = mytof.Tof(1.5, 2., 3., 25.e-9, .999)
    print tof.get_time_v2(pat, field)
    os.system('cp -f tof_save.txt tof_test_9b.txt')

def test9c():
    os.system('rm -rf tof_save.txt')
    field = mytof.Field(1., 1., 2.)
    pat = mytof.Particle(13, -1., 0.1056583745, .10/2**.5, .10/2**.5, .1, 0.0, 0.00, -1.5)
    tof = mytof.Tof(1.5, 2., 3., 25.e-9, .999)
    print tof.get_time_v2(pat, field)
    os.system('cp -f tof_save.txt tof_test_9c.txt')

def test9d():
    os.system('rm -rf tof_save.txt')
    field = mytof.Field(1., 1., 2.)
    pat = mytof.Particle(-13, 1., 0.1056583745, .0, -.1, .1, 0.33, 0., -1.5)
    tof = mytof.Tof(1.5, 2., 3., 25.e-9, .999)
    print tof.get_time_v2(pat, field)
    os.system('cp -f tof_save.txt tof_test_9d.txt')

def test9e():
    os.system('rm -rf tof_save.txt')
    field = mytof.Field(1., 1., 2.)
    pat = mytof.Particle(-13, 1., 0.1056583745, .0, -.1, .1, 0.33, 0., -1.5)
    tof = mytof.Tof(.5, 2., 3., 25.e-9, .999)
    print tof.get_time_v2(pat, field)
    os.system('cp -f tof_save.txt tof_test_9e.txt')

def test9f():
    os.system('rm -rf tof_save.txt')
    field = mytof.Field(1., 1., 2.)
    pat = mytof.Particle(-13, 1., 0.1056583745, -.1, 0., .1, 0., 0.33, -1.5)
    tof = mytof.Tof(.5, 2., 3., 25.e-9, .999)
    print tof.get_time_v2(pat, field)
    os.system('cp -f tof_save.txt tof_test_9f.txt')

def test9g():
    os.system('rm -rf tof_save.txt')
    field = mytof.Field(1., 1., 2.)
    pat = mytof.Particle(13, -1., 0.1056583745, -.1, 0., .1, 0., 0.33, -1.5)
    tof = mytof.Tof(.5, 2., 3., 25.e-9, .999)
    print tof.get_time_v2(pat, field)
    os.system('cp -f tof_save.txt tof_test_9g.txt')

def test9j():
    os.system('rm -rf tof_save.txt')
    field = mytof.Field(1., 1., 2.)
    pat = mytof.Particle(-13, 1., 0.1056583745, -.33, 0.01, .1, 0., 0.33, -1.5)
    tof = mytof.Tof(.5, 2., 3., 25.e-9, .999)
    print tof.get_time_v2(pat, field)
    os.system('cp -f tof_save.txt tof_test_9j.txt')

def test9k():
    os.system('rm -rf tof_save.txt')
    field = mytof.Field(1., 1., 2.)
    pat = mytof.Particle(13, -1., 0.1056583745, -.1, 0., .02, 0., 0.33, -1.5)
    tof = mytof.Tof(.5, 2., 3., 25.e-9, .999)
    print tof.get_time_v2(pat, field)
    os.system('cp -f tof_save.txt tof_test_9k.txt')


def main():
    test9f()
    # return 0
    test6c()
    test3()
    test9j()
    test9g()
    test9d()
    test6b()
    test1()
    test1b()
    test2()
    test4()
    test5()
    test6()
    test7()
    test8()
    test9()
    test9b()
    test9c()
    test9e()
    test9k()


if __name__ == '__main__':
    main()
