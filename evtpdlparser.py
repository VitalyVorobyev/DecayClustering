# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 18:36:38 2016
@author: Vitaly Vorobyev
"""

import pandas as pd
import os.path

evtpdl        = pd.DataFrame()
numpdl        = 0
pid_tuple     = ()
pname_tuple   = ()
pid_name_dict = {}
pname_id_dict = {}

def update_evtpdl(infile = "evtpdl.csv"):
    global evtpdl
    global numpdl
    global pid_tuple
    global pname_tuple
    global pid_name_dict
    global pname_id_dict
    evtpdl        = pd.read_csv(infile)
    numpdl        = evtpdl.shape[0]
    pid_tuple     = tuple(evtpdl['id'])
    pname_tuple   = tuple(evtpdl['name'])
    pid_name_dict = dict(zip(pid_tuple,pname_tuple))
    pname_id_dict = dict(zip(pname_tuple,pid_tuple))

if os.path.isfile("evtpdl.csv"):
    update_evtpdl()
    print("evtpdl is initialized with evtpdl.csv, table size is %d" % numpdl)
else:
    print("!!! Cannot initialize evtpdl table !!!")

def particle_name(idhep):
    if idhep in pid_tuple:
        return pid_name_dict[idhep]
    else:
        return "noname"
        
def particle_id(name):
    if name in pname_tuple:
        return pname_id_dict[name]
    else:
        return -1

def generate_evtpdl():
    """ Parsing the evl.pdl file and save the table in evtpdl.csv
    """
    flag = True
    while flag:
        print("Give me evt.pdl, please...")
        evtpdl = raw_input().strip()
        if evtpdl == "":
            evtpdl = "/home/vitaly/EvtGen/EvtGen/R01-06-00/evt.pdl"
            if os.path.isfile(evtpdl):
                print("Thanks.")
                inputfile = open(evtpdl)
                break
            else:
                print("The file does not exist. Try again...")
        
        labels = tuple(['name',
                        'id',
                        'mass/GeV',
                        'width/GeV',
                        'max_Dm/GeV',
                        '3*charge',
                        '2*spin',
                        'lifetime*c/mm',
                        'PythiaId'])

    df = pd.DataFrame([dict(zip(labels,x.split()[3:])) for x in inputfile if len(x.split()) == 12])
    inputfile.close()
    print(df.shape)
    df.to_csv("evtpdl.csv")
    update_evtpdl()
    