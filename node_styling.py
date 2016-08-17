# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 08:21:14 2016

@author: Vitaly Vorobyev
"""
import evtpdlparser as pdl
"""                [ EvtGen ID    : ( Node Shape,      Font Color,       Fill Color) ]  """
node_format = {
 pdl.pname_id_dict['Upsilon(5S)'] : ( 'tripleoctagon', 'black'         , 'gold')          , # 9000553
 pdl.pname_id_dict['Upsilon(4S)'] : ( 'doubleoctagon', '#006699'       , 'gold')          , # 300553
 pdl.pname_id_dict['Upsilon(3S)'] : ( 'octagon'      , 'black'         , 'gold')          , # 200553
 pdl.pname_id_dict['Upsilon(2S)'] : ( 'octagon'      , 'black'         , 'gold')          , # 100553
 pdl.pname_id_dict['Upsilon']     : ( 'octagon'      , 'black'         , 'gold')          , # 553
 pdl.pname_id_dict['B0']          : ( 'polygon'      , 'black'         , 'antiquewhite')  , # 511
 pdl.pname_id_dict['B+']          : ( 'polygon'      , 'black'         , 'antiquewhite')  , # 521
 pdl.pname_id_dict['B_s0']        : ( 'polygon'      , 'black'         , 'antiquewhite')  , # 531
 pdl.pname_id_dict['D+']          : ( 'polygon'      , 'black'         , 'azure')         , # 411
 pdl.pname_id_dict['D0']          : ( 'polygon'      , 'black'         , 'azure')         , # 421
 pdl.pname_id_dict['D_s+']        : ( 'polygon'      , 'black'         , 'beige')         , # 431
 pdl.pname_id_dict['J/psi']       : ( 'octagon'      , 'black'         , 'gold')          , # 443
 pdl.pname_id_dict['psi(2S)']     : ( 'octagon'      , 'black'         , 'gold')          , # 100443
 pdl.pname_id_dict['D*+']         : ( 'invtriangle'  , 'black'         , 'cyan')          , # 413
 pdl.pname_id_dict['D*0']         : ( 'invtriangle'  , 'black'         , 'cyan')          , # 423
 pdl.pname_id_dict['D_s*+']       : ( 'invtriangle'  , 'black'         , 'cyan')          , # 433
 pdl.pname_id_dict['psi(3770)']   : ( 'octagon'      , 'black'         , 'gold')          , # 30443
 pdl.pname_id_dict['K_S0']        : ( 'doublecircle' , 'black'         , 'cornflowerblue'), # 310
 pdl.pname_id_dict['K_L0']        : ( 'doublecircle' , 'black'         , 'cornflowerblue'), # 130
 pdl.pname_id_dict['K+']          : ( 'oval'         , 'black'         , 'deeppink')      , # 321
 pdl.pname_id_dict['pi+']         : ( 'oval'         , 'black'         , 'deepskyblue')   , # 211
 pdl.pname_id_dict['pi0']         : ( 'oval'         , 'black'         , 'coral')         , # 111
 pdl.pname_id_dict['p+']          : ( 'pentagon'     , 'black'         , 'darkseagreen1') , # 2212
 pdl.pname_id_dict['n0']          : ( 'pentagon'     , 'black'         , 'darkseagreen1') , # 2112
 pdl.pname_id_dict['e+']          : ( 'circle'       , 'black'         , 'floralwhite')   , # 11
 pdl.pname_id_dict['mu+']         : ( 'circle'       , 'black'         , 'floralwhite')   , # 13
 pdl.pname_id_dict['tau+']        : ( 'circle'       , 'black'         , 'floralwhite')   , # 15
 pdl.pname_id_dict['gamma']       : ( 'oval'         , 'black'         , 'chartreuse')      # 22
}

def node_style(idhep):
    nst = {'label': pdl.particle_name(idhep)}
    idh = abs(idhep)
    if idh in node_format:
        nst['style']     = 'filled'
        nst['shape']     = node_format[idh][0]
        nst['fontcolor'] = node_format[idh][1]
        nst['fillcolor'] = node_format[idh][2]
    return nst