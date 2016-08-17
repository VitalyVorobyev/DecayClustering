# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 17:47:11 2016

@author: Vitaly Vorobyev
"""

import ROOT             as r
#import decay_tree_tools as dtt

cuts = {
 "de"      : { "s"  : (0     , 0.025)},
 "mbc"     : { "s"  : (5.2795, 0.006)},
 "md0"     : { "s"  : (1.865 , 0.012)},
 "r_hp_d0" : { "s"  : (0     , 0.04 )},
 "r_hm_d0" : { "s"  : (0     , 0.04 )},
 "z_hp_d0" : { "s"  : (0     , 1.50 )},
 "z_hm_d0" : { "s"  : (0     , 1.50 )},
 "mdpip"   : { "si" : (4.0412, 0.015),  # D* veto
               "a"  : (4.    , 27   )}, # only converged mass fits
 "mdpim"   : { "si" : (4.0412, 0.015),  # D* veto
               "a"  : (4.    , 27   )},
}

def asymcutstr(par,cuts,inv = False):
    g1, g2, g3, bra, ket = ("<", ">", " || ", "(", ")") if inv else (">", "<", " && ", "", "")
    return bra + par + g1 + str(cuts[0]) + g3 + par + g2 + str(cuts[1]) + ket

def symcutstr(par,cuts,inv = False):
    g1 = ">" if inv else "<"
    if cuts[0] == 0:
        return "abs(" + par + ")" + g1 + str(cuts[1])
    else:
        return "abs(" + par + "-" + str(cuts[0]) + ")" + g1 + str(cuts[1])

def cutstr(par,cutd):
    cuts = []
    for c in cutd:
        if c[0] == "a":
            cuts.append(asymcutstr(par,cutd[c],c == "ai"))
        elif c[0] == "s":
            cuts.append(symcutstr( par,cutd[c],c == "si"))
    return " && ".join(cuts)

selections = " && ".join([cutstr(par,cuts[par]) for par in cuts])
print(selections)

#filename = "/home/vitaly/B0toD0pipi/gen_mixed_s0_kpi.root"
filename = "/home/vitaly/B0toD0pipi/Analysis/data/gen_mixed_s0_kpi.root"
##filename = "/home/vitaly/B0toD0pipi/sigmc_kpi_fil.root"
rfile  = r.TFile(filename)
intree = rfile.Get('TEvent')

ofile = r.TFile("test_tree.root","RECREATE")
outtree = intree.CopyTree(selections)
print(outtree.GetEntries())
outtree.Write()
ofile.Close()
rfile.Close()
#
#selections = "abs(mbc-5.2795)<0.006 && abs(de)<0.02"
#
#b_gv_decays = dtt.b_meson_gv_decay_trees(intree,5555)
#for gv in b_gv_decays:
#    name = "test-output/decaytree.gv" + str(b_gv_decays.index(gv)+1)
#    gv.render(name, view=True)
#
#decay_tree = dtt.full_gv_decay_tree(intree,5555)
#decay_tree.render('test-output/decaytree.gv', view=True)
#
#rfile.Close()
