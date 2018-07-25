#! /usr/bin/env python

import ROOT as r
# import os
import sys
# import glob
# import math
import numpy as np
import array
import tof as mytof


def beforeFill_pid(n0, pdg0, q0, en0, px0, py0, pz0, tof0, mu0, tof, mu, tofdet, mu0det, field):
    for i in xrange(n0[0]):
        ptcl = mytof.Particle(pdg0[i], q0[i], mytof.MASS[pdg0[i]], px0[i], py0[i], pz0[i])
        tof0[i] = tofdet.get_time( ptcl, field )
        mu0[i] = mu0det.get_time( ptcl, field )


def bbeforeFill(evn, tin, n0, pdg0, q0, en0, px0, py0, pz0):
    tin.GetEntry(evn-1)
    n0[0] = 0
    for i, v in enumerate(tin.allGenParticles):
        pdg0[i] = v.core.pdgId
        q0[i] = v.core.charge
        px0[i] = v.core.p4.px
        py0[i] = v.core.p4.py
        pz0[i] = v.core.p4.pz
        en0[i] = ( v.core.p4.mass**2 + px0[i]**2 + py0[i]**2 + pz0[i]**2 )**.5
        n0[0] += 1


def beforeFill(evn, tin, genver, nv0, vx0, vy0, vz0, n, recgen, px, py, pz, n0, px0, py0, pz0):
    tin.GetEntry(evn-1)
    vert = getattr(tin, 'allGenParticles#0')
    for i, v in enumerate(vert):
        genver[i] = v.index
    nv0[0] = 0
    for i, v in enumerate(tin.allGenVertices):
        vx0[i] = v.position.x
        vy0[i] = v.position.y
        vz0[i] = v.position.z
        nv0[0] += 1
    chi2 = []
    for i in xrange(n[0]):
        chi2.append([])
        for j in xrange(n0[0]):
            chi2[-1].append( (px[i]-px0[j])**2 + (py[i]-py0[j])**2 + (pz[i]-pz0[j])**2 )
        recgen[i] = chi2[-1].index( min(chi2[-1]) )


def main():
    #
    # PID systems
    # create ToF and mu as a ToF.
    field = mytof.Field(1., 1.575, 3.780/2)
    tofdet = mytof.Tof(.250, 1.090, 2.550/2, .1, .99)
    mu0det = mytof.Tof(.950, 1.930, 3.900/2, .1, .99)
    #
    Fin = r.TFile('/home/vvorob/public/tuples/fccedm/dkspipi.root')
    tin = Fin.Get('events')
    #
    finname = 'output/tmp_out.txt'
    Fout = r.TFile('output/txt2root.root', 'recreate')
    t = r.TTree('t', 'ctau papas tree :)')
    #
    n = array.array('i', [0])
    recgen = array.array('i', [0]*100)
    pdg = array.array('i', [0]*100)
    q = array.array('f', [0.]*100)
    en = array.array('f', [0.]*100)
    pt = array.array('f', [0.]*100)
    px = array.array('f', [0.]*100)
    py = array.array('f', [0.]*100)
    pz = array.array('f', [0.]*100)
    th = array.array('f', [0.]*100)
    phi = array.array('f', [0.]*100)
    tof = array.array('f', [0.]*100)
    tofmu = array.array('f', [0.]*100)
    #
    t.Branch('n', n, 'n/I')
    t.Branch('recgen', recgen, 'recgen[n]/I')
    t.Branch('pdg', pdg, 'pdg[n]/I')
    t.Branch('q', q, 'q[n]/F')
    t.Branch('en', en, 'en[n]/F')
    t.Branch('pt', pt, 'pt[n]/F')
    t.Branch('px', px, 'px[n]/F')
    t.Branch('py', py, 'py[n]/F')
    t.Branch('pz', pz, 'pz[n]/F')
    t.Branch('th', th, 'th[n]/F')
    t.Branch('phi', phi, 'phi[n]/F')
    t.Branch('tof', tof, 'tof[n]/F')
    t.Branch('tofmu', tofmu, 'tofmu[n]/F')
    #
    n0 = array.array('i', [0])
    pdg0 = array.array('i', [0]*100)
    genver = array.array('i', [0]*100)
    q0 = array.array('f', [0.]*100)
    en0 = array.array('f', [0.]*100)
    pt0 = array.array('f', [0.]*100)
    px0 = array.array('f', [0.]*100)
    py0 = array.array('f', [0.]*100)
    pz0 = array.array('f', [0.]*100)
    th0 = array.array('f', [0.]*100)
    phi0 = array.array('f', [0.]*100)
    tof0 = array.array('f', [0.]*100)
    tofmu0 = array.array('f', [0.]*100)
    #
    t.Branch('n0', n0, 'n0/I')
    t.Branch('genver', genver, 'genver[n0]/I')
    t.Branch('pdg0', pdg0, 'pdg0[n0]/I')
    t.Branch('q0', q0, 'q0[n0]/F')
    t.Branch('en0', en0, 'en0[n0]/F')
    # t.Branch('pt0', pt0, 'pt0[n0]/F')
    t.Branch('px0', px0, 'px0[n0]/F')
    t.Branch('py0', py0, 'py0[n0]/F')
    t.Branch('pz0', pz0, 'pz0[n0]/F')
    # t.Branch('th0', th0, 'th0[n0]/F')
    # t.Branch('phi0', phi0, 'phi0[n0]/F')
    t.Branch('tof0', tof0, 'tof0[n]/F')
    t.Branch('tofmu0', tofmu0, 'tofmu0[n]/F')
    #
    nv0 = array.array('i', [0])
    vx0 = array.array('f', [0.]*100)
    vy0 = array.array('f', [0.]*100)
    vz0 = array.array('f', [0.]*100)
    #
    t.Branch('nv0', nv0, 'nv0/I')
    t.Branch('vx0', vx0, 'vx0[nv0]/F')
    t.Branch('vy0', vy0, 'vy0[nv0]/F')
    t.Branch('vz0', vz0, 'vz0[nv0]/F')
    #
    # type evn prt pdgid en pt px py pz th eta phi  m  q
    #        0   1     2  3  4  5  6  7  8   9 10  11 12
    fin = open(finname)
    evn = 1
    n[0] = 0
    n0[0] = 0
    for line in fin:
        if line[0] == '#':
            continue
        # print line.rstrip('\n')
        lin = line.split()
        li = [float(x) for x in lin[4:]]
        li = [int(x) for x in lin[1:4]] + li
        if lin[0] == 'gen':
            isGen = True
            isRec = False
        elif lin[0] == 'rec':
            isRec = True
            isGen = False
        if isRec:
            pdg[n[0]] = li[2]
            en[n[0]] = li[3]
            pt[n[0]] = li[4]
            px[n[0]] = li[5]
            py[n[0]] = li[6]
            pz[n[0]] = li[7]
            th[n[0]] = li[8]
            phi[n[0]] = li[10]
            q[n[0]] = li[12]
            n[0] += 1
        ########
        if isGen:
            if evn != li[0] and evn != -1:
                ########
                bbeforeFill(evn, tin, n0, pdg0, q0, en0, px0, py0, pz0)
                beforeFill(evn, tin, genver, nv0, vx0, vy0, vz0, n, recgen, px, py, pz, n0, px0, py0, pz0)
                beforeFill_pid(n0, pdg0, q0, en0, px0, py0, pz0, tof0, tofmu0, tof, tofmu, tofdet, mu0det, field)
                ########
                t.Fill()
                evn = li[0]
                n[0] = 0
                n0[0] = 0
        ########
    bbeforeFill(evn, tin, n0, pdg0, q0, en0, px0, py0, pz0)
    beforeFill(evn, tin, genver, nv0, vx0, vy0, vz0, n, recgen, px, py, pz, n0, px0, py0, pz0)
    beforeFill_pid(n0, pdg0, q0, en0, px0, py0, pz0, tof0, tofmu0, tof, tofmu, tofdet, mu0det, field)
    t.Fill()
    t.Write()
    Fin.Close()
    Fout.Close()

if __name__ == '__main__':
    main()
