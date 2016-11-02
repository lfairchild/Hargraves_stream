#!/usr/bin/env python
import os
import sys
from time import time
from time import sleep
if sys.version_info[0] < 3:
    from StringIO import StringIO
else:
    from io import StringIO
import pmagpy.ipmag as ipmag
import pmagpy.pmag as pmag
import pmagpy.pmagplotlib as pmagplotlib
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pmagpy.demag_gui_utilities import *
import mpld3
from mpld3 import plugins

def magic_write(ofile,Recs,file_type):
    """
    called by magic_write(outputfile,records_list,magic_file_type)
    writes out a magic format list of dictionaries to ofile

    """
    if len(Recs)<1:
        return False, 'No records to write to file {}'.format(ofile)
    #else:
        #print len(Recs),' records written to file ',ofile
    pmag_out=''
    outstring="tab \t"+file_type+"\n"
    pmag_out+=outstring
    keystring=""
    keylist=[]
    for key in Recs[0].keys():
        keylist.append(key)
    keylist.sort()
    for key in keylist:
        keystring=keystring+'\t'+key.strip()
    keystring=keystring + '\n'
    pmag_out+=keystring[1:]
    for Rec in Recs:
        outstring=""
        for key in keylist:
           try:
              outstring=outstring+'\t'+str(Rec[key]).strip()
           except KeyError:
              if 'er_specimen_name' in Rec.keys():
                  print Rec['er_specimen_name']
              elif 'er_specimen_names' in Rec.keys():
                  print Rec['er_specimen_names']
              print("No data for %s"%key)
              #raw_input()
        outstring=outstring+'\n'
        pmag_out+=outstring[1:]
    return str(pmag_out)


def cit_magic(command_line=True, **kwargs):
    """
    NAME
        cit_magic.py

    DESCRIPTION
        converts CIT and .sam  format files to magic_measurements format files

    SYNTAX
        cit_magic.py [command line options]

    OPTIONS
        -h: prints the help message and quits.
        -usr USER:   identify user, default is ""
        -f FILE: specify .sam format input file, required
        -fsi SITEFILE : specify file with site names and locations [tab delimited magic file]
        -F FILE: specify output  measurements file, default is magic_measurements.txt
        -Fsp FILE: specify output er_specimens.txt file, default is er_specimens.txt
        -Fsi FILE: specify output er_sites.txt file, default is er_sites.txt
        -Fsa FILE: specify output er_samples.txt file, default is er_samples.txt  # LORI
        -n [gm,kg,cc,m3]: specify normalization
        -A: don't average replicate measurements
        -spc NUM : specify number of characters to designate a  specimen, default = 0
        -ncn NCON: specify naming convention
        -loc LOCNAME : specify location/study name, must have either LOCNAME or SITEFILE or be a synthetic
        -mcd [FS-FD:SO-MAG,.....] colon delimited list for method codes applied to all specimens in .sam file
        -dc B PHI THETA: dc lab field (in micro tesla) and phi,theta, default is none
              NB: use PHI, THETA = -1 -1 to signal that it changes, i.e. in anisotropy experiment
        -ac B : peak AF field (in mT) for ARM acquisition, default is none

    INPUT
        Best to put separate experiments (all AF, thermal, thellier, trm aquisition, Shaw, etc.)

    NOTES:
         Sample naming convention:
        [1] XXXXY: where XXXX is an arbitrary length site designation and Y
            is the single character sample designation.  e.g., TG001a is the
            first sample from site TG001.    [default]
        [2] XXXX-YY: YY sample from site XXXX (XXX, YY of arbitary length)
        [3] XXXX.YY: YY sample from site XXXX (XXX, YY of arbitary length)
        [4-Z] XXXX[YYY]:  YYY is sample designation with Z characters from site XXX
        [5] site name = sample name
        [6] site name entered in site_name column in the orient.txt format input file  -- NOT CURRENTLY SUPPORTED
        [7-Z] [XXX]YYY:  XXX is site designation with Z characters from samples  XXXYYY
        NB: all others you will have to either customize your
            self or e-mail ltauxe@ucsd.edu for help.
    """
    #
    #initialize variables
    norm='cc'
    samp_con,Z='3',1
    meas_file='magic_measurements.txt'
    spec_file='er_specimens.txt'
    samp_file='er_samples.txt'
    site_file='er_sites.txt'
    ErSpecs,ErSamps,ErSites,ErLocs,ErCits=[],[],[],[],[]
    MeasRecs=[]
    specnum,units,locname=0,"1","unknown"
    citation="This study"
    dir_path='.'
    args=sys.argv
    if command_line:
        if '-WD' in args:
            ind=args.index("-WD")
            dir_path=args[ind+1]
        if "-h" in args:
            print main.__doc__
            return False
        if "-usr" in args:
            ind=args.index("-usr")
            user=args[ind+1]
        if '-F' in args:
            ind=args.index("-F")
            meas_file=args[ind+1]
        if '-Fsp' in args:
            ind=args.index("-Fsp")
            spec_file=args[ind+1]
        if '-Fsa' in args:
            ind=args.index("-Fsa")
            samp_file=args[ind+1]
        if '-Fsi' in args:   # LORI addition
            ind=args.index("-Fsi")
            site_file=args[ind+1]
        if '-loc' in args:
            ind=args.index("-loc")
            locname=args[ind+1]
        if '-mcd' in args:
            ind=args.index("-mcd")
            methods=args[ind+1]
        else:
            methods='SO-MAG'
        if '-spc' in args:
            ind=args.index("-spc")
            specnum=-int(args[ind+1])
        if '-n' in args:
            ind=args.index("-n")
            norm=args[ind+1]
        if "-A" in args:
            avg=1
        else:
            avg=0
        if "-ncn" in args:
            ind=args.index("-ncn")
            samp_con=sys.argv[ind+1]
            if "4" in samp_con:
                if "-" not in samp_con:
                    print "option [4] must be in form 4-Z where Z is an integer"
                    return False, "naming convention option [4] must be in form 4-Z where Z is an integer"
                else:
                    Z=samp_con.split("-")[1]
                    samp_con="4"
        if '-f' in args:
            ind=args.index("-f")
            magfile=args[ind+1]
        if '-ID' in args:
            ind = args.index('-ID')
            input_dir_path = args[ind+1]
        else:
            input_dir_path = dir_path
        output_dir_path = dir_path
        # LJ

    # if you are running as a module:
    elif not command_line:
        dir_path = kwargs.get('dir_path', '.')
        user = kwargs.get('user', '')
        meas_file = kwargs.get('meas_file', 'magic_measurements.txt') # outfile
        spec_file = kwargs.get('spec_file', 'er_specimens.txt') # specimen outfile
        samp_file = kwargs.get('samp_file', 'er_samples.txt') # sample outfile
        site_file = kwargs.get('site_file', 'er_sites.txt') # site outfile
        locname = kwargs.get('locname', '')
        methods = kwargs.get('methods', ['SO-MAG'])
        specnum = -int(kwargs.get('specnum', 0))
        norm = kwargs.get('norm', 'cc')
        avg = kwargs.get('avg', 0)  # 0 means do average, 1 means don't
        samp_con = kwargs.get('samp_con', '3')
        magfile = kwargs.get('magfile', '')
        input_dir_path = kwargs.get('input_dir_path', dir_path)
        output_dir_path = dir_path
    # done with module-specific stuff

    # formatting and checking variables
    if "4" in samp_con:
        if "-" not in samp_con:
            print "option [4] must be in form 4-Z where Z is an integer"
            return False, "naming convention option [4] must be in form 4-Z where Z is an integer"
        else:
            Z=samp_con.split("-")[1]
            samp_con="4"

    magfile = os.path.join(input_dir_path, magfile)
    spec_file = os.path.join(output_dir_path, spec_file)
    samp_file = os.path.join(output_dir_path, samp_file)
    site_file = os.path.join(output_dir_path, site_file)
    meas_file= os.path.join(output_dir_path, meas_file)
    try:
        file_input=open(magfile,'r')
    except Exception as ex:
        print "bad sam file name: ", magfile
        return False, "bad sam file name"
    File = file_input.readlines()
    if len(File) == 1: File = File[0].split('\r'); File = map(lambda x: x+"\r\n", File)
    sids,ln,format=[],0,'CIT'
    formats=['CIT','2G','APP','JRA']
    if File[ln].strip()=='CIT': ln+=1
    ErLocRec={}
    ErLocRec["er_location_name"]=locname
    ErLocRec["er_citation_names"]=citation
    comment=File[ln]
    if comment=='CIT':
       format=comment
       ln+=1
    comment=File[ln]
    #print comment
    ln+=1
    specimens,samples,sites=[],[],[]
    if format=='CIT':
        line=File[ln].split()
        site_lat=line[0]
        site_lon=line[1]
        ErLocRec["location_begin_lat"]=site_lat
        ErLocRec["location_begin_lon"]=site_lon
        ErLocRec["location_end_lat"]=site_lat
        ErLocRec["location_end_lon"]=site_lon
        ErLocs.append(ErLocRec)
        try: Cdec=float(line[2])
        except ValueError: pdb.set_trace()
        for k in range(ln+1,len(File)):
            line=File[k]
            rec=line.split()
            if rec == []: continue
            specimen=rec[0]
            specimens.append(specimen)
    for specimen in specimens:
        ErSpecRec,ErSampRec,ErSiteRec={},{},{}
        if specnum!=0:
            sample=specimen[:specnum]
        else: sample=specimen
        site=pmag.parse_site(sample,samp_con,Z)
        ErSpecRec['er_specimen_name']=specimen
        ErSpecRec['er_sample_name']=sample
        ErSpecRec['er_site_name']=site
        ErSpecRec['er_location_name']=locname
        ErSpecRec['er_citation_names']=citation
        ErSampRec['er_sample_name']=sample
        ErSampRec['er_site_name']=site
        ErSampRec['er_location_name']=locname
        ErSampRec['er_citation_names']=citation
        ErSampRec['magic_method_codes']=methods
        ErSampRec['sample_declination_correction']='%7.1f'%(Cdec)
        ErSiteRec['er_site_name']=site
        ErSiteRec['er_location_name']=locname
        ErSiteRec['er_citation_names']=citation
        ErSiteRec['site_lat']=site_lat
        ErSiteRec['site_lon']=site_lon
        f=open(input_dir_path+'/'+specimen,'rU')
        Lines=f.readlines()
        comment=""
        line=Lines[0].split()
        if len(line)>2:
            comment=line[2]
        info=Lines[1].split()
        vol=float(info[-1])
        if vol!=1.0:
            if norm=='cc':units="1"
            if norm=='m3':units="2"
            ErSpecRec['specimen_weight']=""
            if units=="1" or "":
                ErSpecRec['specimen_volume']='%10.3e'%(vol*1e-6)
            else:
                ErSpecRec['specimen_volume']='%10.3e'%(vol)
        else:
            if norm=='cc':units="1"
            if norm=='m3':units="2"
            ErSpecRec['specimen_volume']=""
            if units=="1" or "":
                ErSpecRec['specimen_weight']='%10.3e'%(vol*1e-3)
            else:
                ErSpecRec['specimen_weight']='%10.3e'%(vol)
        dip=float(info[-2])
        dip_direction=float(info[-3])+Cdec+90.
        sample_dip=-float(info[-4])
        sample_azimuth=float(info[-5])+Cdec-90.
        if len(info)>5:
            ErSampRec['sample_height']=info[-6]
        else:
            ErSampRec['sample_height']='0'
        ErSampRec['sample_azimuth']='%7.1f'%(sample_azimuth)
        ErSampRec['sample_dip']='%7.1f'%(sample_dip)
        ErSampRec['sample_bed_dip']='%7.1f'%(dip)
        ErSampRec['sample_bed_dip_direction']='%7.1f'%(dip_direction)
        ErSampRec['sample_class']=''
        ErSampRec['sample_type']=''
        ErSampRec['sample_lithology']=''
        if Cdec!=0 or Cdec!="":
            ErSampRec['magic_method_codes']='SO-CMD-NORTH'
        else:
            ErSampRec['magic_method_codes']='SO-MAG'
        for line in Lines[2:len(Lines)]:
            #print 'line:', line
            MeasRec=ErSpecRec.copy()

#           Remove specimen_volume and specimen_weight as they do not exits in the magic_measurement table
            del MeasRec["specimen_volume"]
            del MeasRec["specimen_weight"]

            treat_type=line[0:3]
            treat=line[2:6]
            try: float(treat)
            except ValueError: treat = line[3:6]
            if treat_type.startswith('NRM'):
                MeasRec['magic_method_codes']='LT-NO'
                MeasRec['measurement_temp']='273'
                MeasRec['treatment_temp']='273'
                MeasRec['treatment_dc_field']='0'
                MeasRec['treatment_ac_field']='0'
            elif treat_type.startswith('AF'):
                MeasRec['magic_method_codes']='LT-AF-Z'
                MeasRec['measurement_temp']='273'
                MeasRec['treatment_temp']='273'
                MeasRec['treatment_dc_field']='0'
                if treat.strip() == '':
                    MeasRec['treatment_ac_field']='0'
                else:
                    MeasRec['treatment_ac_field']='%10.3e'%(float(treat)*1e-3)
            elif treat_type.startswith('ARM'):
                MeasRec['magic_method_codes']="LP-ARM"
                MeasRec['measurement_temp']='273'
                MeasRec['treatment_temp']='273'
                MeasRec['treatment_dc_field']='0'
                if treat.strip() == '':
                    MeasRec['treatment_ac_field']='0'
                else:
                    MeasRec['magic_method_codes']="LP-ARM-AFD"
                    MeasRec['treatment_ac_field']='%10.3e'%(float(treat)*1e-3)
            elif treat_type.startswith('TT'):
                MeasRec['magic_method_codes']='LT-T-Z'
                MeasRec['measurement_temp']='273'
                if treat.strip() == '':
                    MeasRec['treatment_temp']='273'
                else:
                    MeasRec['treatment_temp']='%7.1f'%(float(treat)+273)
                MeasRec['treatment_dc_field']='0'
                MeasRec['treatment_ac_field']='0'
            elif treat_type.startswith('LT') or treat_type.startswith('LN2'):
                MeasRec['magic_method_codes']='LT-LT-Z'
                MeasRec['measurement_temp']='273'
                MeasRec['treatment_temp']='77'
                MeasRec['treatment_dc_field']='0'
                MeasRec['treatment_ac_field']='0'
            else:
                print "trouble with your treatment steps"
            MeasRec['measurement_dec']=line[46:51]
            MeasRec['measurement_inc']=line[52:58]
            M='%8.2e'%(float(line[31:39])*vol*1e-3) # convert to Am2
            MeasRec['measurement_magn_moment']=M
            MeasRec['measurement_csd']='%7.1f'%(eval(line[41:46]))
            MeasRec["measurement_positions"]='1'
            MeasRec['measurement_standard']='u'
            if len(line)>60:
                MeasRec['magic_instrument_codes']=line[85:]
                MeasRec['measurement_sd_x']='%8.2e'%(float(line[58:67])*1e-8) #(convert e-5emu to Am2)
                MeasRec['measurement_sd_y']='%8.2e'%(float(line[67:76])*1e-8)
                MeasRec['measurement_sd_z']='%8.2e'%(float(line[76:85])*1e-8)
            MeasRecs.append(MeasRec)
        ErSpecs.append(ErSpecRec)
        if sample not in samples:
            samples.append(sample)
            ErSamps.append(ErSampRec)
        site=pmag.parse_site(sample,samp_con,Z)
        if site not in sites:
            sites.append(site)
            ErSites.append(ErSiteRec)
    er_specs = magic_write(spec_file,ErSpecs,'er_specimens')
    #print 'specimens stored in "er_specs"'
    er_samps = magic_write(samp_file,ErSamps,'er_samples')
    #print 'samples stored in "er_samps"'
    er_sites = magic_write(site_file,ErSites,'er_sites')
    #print 'sites stored in "er_sites"'
    Fixed=pmag.measurements_methods(MeasRecs,avg)
    magic_meas = magic_write(meas_file,Fixed,'magic_measurements')
    #print 'data stored in "magic_meas"'
    return True, meas_file, str(er_specs), str(er_samps), str(er_sites), str(magic_meas)

def plotZ(datablock,angle,s,norm):
    """
    function to make Zijderveld diagrams
    """
    amin,amax=0.,-100.
    fact=1./(datablock[0][3])   # normalize to NRM=1
    if norm==0:fact=1.
    x,y,z=[],[],[]
    xb,yb,zb=[],[],[]
    forVDS=[]
    labels = []
# convert to cartesian
    recnum,delta=0,""
    for plotrec in datablock:
        forVDS.append([plotrec[1],plotrec[2],plotrec[3]/datablock[0][3]])
        rec= pmag.dir2cart([(plotrec[1]-angle),plotrec[2],plotrec[3]*fact])
        if len(plotrec)==4:plotrec.append('0') # fake the ZI,IZ step for old data
        if len(plotrec)==5:plotrec.append('g') # assume good measurement if not specified
        if plotrec[5]=='g':
          #  z.append(-rec[2])
            z.append(rec[2])
            x.append(rec[0])
          #  y.append(-rec[1])
            y.append(rec[1])
            if x[-1]>amax:amax=x[-1]
            if y[-1]>amax:amax=y[-1]
            if z[-1]>amax:amax=z[-1]
            if x[-1]<amin:amin=x[-1]
            if y[-1]<amin:amin=y[-1]
            if z[-1]<amin:amin=z[-1]
            if delta=="":delta=.02*x[-1]
            if recnum%2==0 and len(x)>0:
                # labels.append(' '*3 +str(int(datablock[recnum][0] - 273))+'$\degree$C')
                plt.text(x[-1]-delta,z[-1]+delta,' '*3 +str(int(datablock[recnum][0] - 273))+'$\degree$C',fontsize=9)
            recnum+=1
        elif len(plotrec)>=6 and plotrec[5]=='b':
          #  zb.append(-rec[2])
            zb.append(rec[2])
            xb.append(rec[0])
          #  yb.append(-rec[1])
            yb.append(rec[1])
            if xb[-1]>amax:amax=xb[-1]
            if yb[-1]>amax:amax=yb[-1]
            if zb[-1]>amax:amax=zb[-1]
            if xb[-1]<amin:amin=xb[-1]
            if yb[-1]<amin:amin=yb[-1]
            if zb[-1]<amin:amin=zb[-1]
            if delta=="":delta=.02*xb[-1]
            # labels.append(' '*3 +str(int(datablock[recnum][0] - 273))+'$\degree$C')
            plt.text(xb[-1]-delta,zb[-1]+delta,' '*3 + str(int(datablock[recnum][0] - 273))+'$\degree$C',fontsize=9)
            recnum+=1
# plotting stuff
    if angle !=0:tempstr= "\n Declination rotated by: "+str(angle)+'\n'
#         globals.text.insert(globals.END,tempstr)
    Zlist  = x
    Zlisty = y
    Zlistz = z
    if len(xb)>0:
        plt.scatter(xb,yb,marker='d',c='w',s=30)
        plt.scatter(xb,zb,marker='d',c='w',s=30)
    plt.plot(x,y,'r')
    plt.plot(x,z,'b')
    dec_points = plt.scatter(x,y,marker='o',c='r')
    inc_points = plt.scatter(x,z,marker='s',c='w')
    xline=[amin,amax]
   # yline=[-amax,-amin]
    yline=[amax,amin]
    zline=[0,0]
    plt.plot(xline,zline, c='k')
    plt.plot(zline,xline, c='k')
    # if angle!=0:xlab="X: rotated to Dec = "+'%7.1f'%(angle)
    # if angle==0:xlab="X: rotated to Dec = "+'%7.1f'%(angle)
    # plt.xlabel(xlab)
    # plt.ylabel("Circles: Y; Squares: Z")
    tstring=s+': NRM = '+'%9.2e'%(datablock[0][3])
    plt.axis([amin,amax,amax,amin])
    plt.axis("equal")
    plt.axis('off')
    plt.tight_layout()
    plt.title(tstring)
    # tooltip_dec = plugins.PointHTMLTooltip(dec_points, labels, voffset=10, hoffset=10)
    # tooltip_inc = plugins.PointHTMLTooltip(inc_points, labels, voffset=10, hoffset=10)
    # plugins.connect(plt.gcf(), tooltip_dec)
    # plugins.connect(plt.gcf(), tooltip_inc)
    # mpld3.enable_notebook()


def plotMT(datablock,s,num,units,norm):
    Ints=[]
    for plotrec in datablock:
        Ints.append(plotrec[3])
    Ints.sort()
    T,M,Tv,recnum=[],[],[],0
    Mex,Tex,Vdif=[],[],[]
    recbak=[]
    for rec in datablock:
        if rec[5]=='g':
            if units=="T":
                T.append(rec[0]*1e3)
                Tv.append(rec[0]*1e3)
                if recnum>0:Tv.append(rec[0]*1e3)
            elif units=="U":
                T.append(rec[0])
                Tv.append(rec[0])
                if recnum>0:Tv.append(rec[0])
            elif units=="K":
                T.append(rec[0]-273)
                Tv.append(rec[0]-273)
                if recnum>0:Tv.append(rec[0]-273)
            elif "T" in units and "K" in units:
                if rec[0]<1.:
                    T.append(rec[0]*1e3)
                    Tv.append(rec[0]*1e3)
                else:
                    T.append(rec[0]-273)
                    Tv.append(rec[0]-273)
                    if recnum>0:Tv.append(rec[0]-273)
            else:
                T.append(rec[0])
                Tv.append(rec[0])
                if recnum>0:Tv.append(rec[0])
            if norm==1:
                M.append(rec[3]/Ints[-1])
            else:
                M.append(rec[3])
            if recnum>0 and len(rec)>0 and len(recbak)>0:
                v=[]
                if recbak[0]!=rec[0]:
                    V0=pmag.dir2cart([recbak[1],recbak[2],recbak[3]])
                    V1=pmag.dir2cart([rec[1],rec[2],rec[3]])
                    for el in range(3):v.append(abs(V1[el]-V0[el]))
                    vdir=pmag.cart2dir(v)
                    Vdif.append(vdir[2]/Ints[-1]) # append vector difference
                    Vdif.append(vdir[2]/Ints[-1]) #
            recbak=[]
            for el in rec: recbak.append(el)
            delta=.02*M[0]
            if num==1:
                if recnum%2==0: plt.text(T[-1]+delta,M[-1],' '*3 + str(int(datablock[recnum][0] - 273))+'$\degree$C',fontsize=9)
            recnum+=1
        else:
            if rec[0]<200:Tex.append(rec[0]*1e3)
            if rec[0]>=200:Tex.append(rec[0]-273)
            Mex.append(rec[3]/Ints[-1])
            recnum+=1
        MTlist =T
        MTlisty=M
    if len(Mex)>0 and len(Tex)>0:
        plt.scatter(Tex,Mex,marker='d',color='k')
    if len(Vdif)>0:
        Vdif.append(vdir[2]/Ints[-1]) #
        Vdif.append(0)
    Tv.append(Tv[-1])
    plt.plot(T,M)
    plt.plot(T,M,'ro')
    if len(Tv)==len(Vdif) and norm==1:plt.plot(Tv,Vdif,'g-')
    if units=="T":
        plt.xlabel("Step (mT)")
    elif units=="K":
        plt.xlabel("Step (C)")
    elif units=="J":
        plt.xlabel("Step (J)")
    else:
        plt.xlabel("Step [mT,C]")
    if norm==1:plt.ylabel("Fractional Magnetization")
    if norm==0:plt.ylabel("Magnetization")
    plt.axvline(0,color='k')
    plt.axhline(0,color='k')
    tstring= s
    plt.title(tstring)

def generate_inp_file(full_sam_path):
    """
    DESCRIPTION
        Uses sample and site DataFrames from mk_sam_file.main function to generate

        @param: od - output directory
        @param: df - sample Dataframe
        @param: hdf - site DataFrame

    OUTPUT
        .inp file

    """
    sam_file = open(full_sam_path,'r')
    sam_doc = sam_file.read()
    sam_file.close()
    if sam_doc.find('\r\n') != -1:
        sam_doc = sam_doc.replace('\r\n','\n')
    else:
        sam_doc = sam_doc.replace('\r','\n')

# #     od = os.path.join(os.getcwd(), od)
#     if od != '' and not od.endswith('/'):
#         od += '/'

    inps = ""

    inps += "CIT\n"
    inps += "sam_path\tfield_magic_codes\tlocation\tnaming_convention\tnum_terminal_char\tdont_average_replicate_measurements\tpeak_AF\ttime_stamp\n"
    inps += (full_sam_path + '\t')
#     if all(df.T['comment'] == 'sun compass orientation'): inps += 'SO-SUN\t'
#     elif all(df.T['comment'] == 'mag compass orientation (IGRF corrected)'): inps += 'SO-MAG\t'
#     else: inps += 'SO-SM\t'
    inps += 'SO-SM\t'
    site_name = full_sam_path.split('/')[-2]
    inps += (site_name + '\t')

    first_sample_id = sam_doc.split('\n')[3]

    if '-' in first_sample_id:
        inps += '2\t'
        if '-' in site_name:
            first_sample_id = first_sample_id.replace(site_name,'')
        else:
            first_sample_id = first_sample_id.replace(site_name+'-','')
    elif '.' in first_sample_id:
        inps += '3\t'
        if '.' in site_name:
            first_sample_id = first_sample_id.replace(site_name,'')
        else:
            first_sample_id = first_sample_id.replace(site_name+'.','')
    elif site_name == first_sample_id:
        inps += '5\t'
    else: inps += '4\t'

    inps += str(len(first_sample_id)-1) + '\t'
    inps += "True\t"
    inps += "None\t"
    inps += '0.0\n'

    return inps

#     print('Writing file - ' + od + hdf['site_info']['site_id'] + '.inp')
#     if od != '' and not os.path.exists(od):
#         os.makedirs(od)
#     inpf = open(od + hdf['site_info']['site_id'] + '.inp', 'w+')
#     inpf.write(inps)
#     inpf.close()

def read_inp(inp_file_name, full_sam_path='.'):
    magic_files = []
    # try:
    #     inp_file = open(inp_file_name, "r")
    #     lines = inp_file.read().split("\n")
    # except:
    #     inp_file = generate_inp_file(full_sam_path)
    #     lines = inp_file.split("\n")
    inp_file = generate_inp_file(full_sam_path)
    lines = inp_file.split("\n")
    new_inp_file = ""

    if len(lines) < 3: print(".inp file improperly formated"); return
    new_inp_file = lines[0] + "\n" + lines[1] + "\n"
    [lines.remove('') for i in range(lines.count(''))]
    format = lines[0].strip()
    header = lines[1].split('\t')
    update_files = lines[2:]
    update_data = False
    for i,update_file in enumerate(update_files):
        update_lines = update_file.split('\t')
        if not os.path.isfile(update_lines[0]):
            #print("%s not found searching for location of file"%(update_lines[0]))
            sam_file_name = os.path.split(update_lines[0])[-1]
            new_file_path = find_file(sam_file_name, os.path.curdir)
            if new_file_path == None or not os.path.isfile(new_file_path):
                print("%s does not exist in any subdirectory of %s and will be skipped"%(update_lines[0], self.WD))
                new_inp_file += update_file+"\n"
                continue
            else:
                #print("new location for file found at %s"%(new_file_path))
                update_lines[0] = new_file_path
        d = reduce(lambda x,y: x+"/"+y, update_lines[0].split("/")[:-1])+"/"
        f = update_lines[0].split("/")[-1].split(".")[0] + ".magic"
        if (d+f) in magic_files:
            new_inp_file += update_file+"\n"
            continue
        if float(update_lines[-1]) >= os.path.getctime(update_lines[0]):
            if os.path.isfile(d+f):
                magic_files.append(d+f)
                new_inp_file += update_file+"\n"
                continue
        if len(header) != len(update_lines):
            print("length of header and length of enteries for the file %s are different and will be skipped"%(update_lines[0]))
            new_inp_file += update_file+"\n"
            continue
        update_dict = {}
        for head,entry in zip(header,update_lines):
            update_dict[head] = entry
        if format == "CIT":
            CIT_kwargs = {}
            CIT_name = update_dict["sam_path"].split("/")[-1].split(".|-")[0]
            # print CIT_name
#             os.mkdir(os.path.abspath(os.path.dirname(inp_file_name)) + '/.tmp/')
            CIT_kwargs["dir_path"] = os.path.abspath(os.path.dirname(full_sam_path))# + '/.tmp/'
            CIT_kwargs["user"] = ""
            CIT_kwargs["meas_file"] = CIT_name + ".magic"
            CIT_kwargs["spec_file"] = CIT_name + "_er_specimens.txt"
            CIT_kwargs["samp_file"] = CIT_name + "_er_samples.txt"
            CIT_kwargs["site_file"] = CIT_name + "_er_sites.txt"
            CIT_kwargs["locname"] = update_dict["location"]
            CIT_kwargs["methods"] = update_dict["field_magic_codes"]
            CIT_kwargs["specnum"] = update_dict["num_terminal_char"]
            CIT_kwargs["avg"] = update_dict["dont_average_replicate_measurements"]
            CIT_kwargs["samp_con"] = update_dict["naming_convention"]
            CIT_kwargs["peak_AF"] = update_dict["peak_AF"]
            CIT_kwargs["magfile"] = update_dict["sam_path"].split("/")[-1]
            CIT_kwargs["input_dir_path"] = reduce(lambda x,y: x+"/"+y, update_dict["sam_path"].split("/")[:-1])

            program_ran, error_message, er_specs_str, er_samps_str, er_sites_str, magic_measurements_str = cit_magic(command_line=False, **CIT_kwargs)

            if program_ran:
                update_data = True
                update_lines[-1] = time()
                new_inp_file += reduce(lambda x,y: str(x)+"\t"+str(y), update_lines)+"\n"
                magic_files.append(CIT_kwargs["dir_path"]+CIT_kwargs["meas_file"])
            else:
                new_inp_file += update_file
                if os.path.isfile(CIT_kwargs["dir_path"]+CIT_kwargs["meas_file"]):
                    magic_files.append(CIT_kwargs["dir_path"]+CIT_kwargs["meas_file"])
    try:
        inp_file.close()
    except:
        pass

    er_specimens = pd.read_csv(StringIO(er_specs_str), sep='\t', skiprows=1)
    er_samples = pd.read_csv(StringIO(er_samps_str), sep='\t', skiprows=1)
    er_sites = pd.read_csv(StringIO(er_sites_str), sep='\t', skiprows=1)
    magic_measurements = pd.read_csv(StringIO(magic_measurements_str), sep='\t', skiprows=1)

    return er_specimens, er_sites, er_samples, magic_measurements


#     out_file = open(inp_file_name, "w")
#     out_file.write(new_inp_file)
    # return update_data, er_specimens, er_sites, er_samples, magic_measurements

def plot_net():
    Dcirc=np.arange(0,361.)
    Icirc=np.zeros(361,'f')
    Xcirc,Ycirc=[],[]
    for k in range(361):
        XY= pmag.dimap(Dcirc[k],Icirc[k])
        Xcirc.append(XY[0])
        Ycirc.append(XY[1])
    plt.plot(Xcirc,Ycirc,'k')

    # put on the tick marks
    Xsym,Ysym=[],[]
    for I in range(10,100,10):
        XY=pmag.dimap(0.,I)
        Xsym.append(XY[0])
        Ysym.append(XY[1])
    plt.plot(Xsym,Ysym,'k+')
    Xsym,Ysym=[],[]
    for I in range(10,90,10):
        XY=pmag.dimap(90.,I)
        Xsym.append(XY[0])
        Ysym.append(XY[1])
    plt.plot(Xsym,Ysym,'k+')
    Xsym,Ysym=[],[]
    for I in range(10,90,10):
        XY=pmag.dimap(180.,I)
        Xsym.append(XY[0])
        Ysym.append(XY[1])
    plt.plot(Xsym,Ysym,'k+')
    Xsym,Ysym=[],[]
    for I in range(10,90,10):
        XY=pmag.dimap(270.,I)
        Xsym.append(XY[0])
        Ysym.append(XY[1])
    plt.plot(Xsym,Ysym,'k+')
    for D in range(0,360,10):
        Xtick,Ytick=[],[]
        for I in range(4):
            XY=pmag.dimap(D,I)
            Xtick.append(XY[0])
            Ytick.append(XY[1])
        plt.plot(Xtick,Ytick,'k')
    plt.axis("equal")
    plt.axis((-1.05,1.05,-1.05,1.05))
    plt.tight_layout()

def do_help():
    return main.__doc__
