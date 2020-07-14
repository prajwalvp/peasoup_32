import glob
import xml.etree.ElementTree as ET
import subprocess
import optparse
from optparse import OptionParser
import numpy as np
import sys
import time
import logging
import commands
import re


log = logging.getLogger('manual_presto_fold')
FORMAT = "[%(levelname)s - %(asctime)s - %(filename)s:%(lineno)s] %(message)s"
logging.basicConfig(format=FORMAT)

LIGHT_SPEED = 2.99792458e8                 # Speed of Light in SI

def period_modified(p0,pdot,no_of_samples,tsamp,fft_size):
    if (fft_size==0.0):
        return p0 - pdot*float(1<<(no_of_samples.bit_length()-1))*tsamp/2
    else:
        return p0 - pdot*float(fft_size)*tsamp/2

def period_modified_new(p0,pdot,pdd,no_of_samples,tsamp,fft_size):
    if (fft_size==0.0):
        T_by_2 = float(1<<(no_of_samples.bit_length()-1))*tsamp/2
        return p0 - pdot*T_by_2 - pdd*0.5*T_by_2*T_by_2
    else:
        T_by_2 = float(fft_size)*tsamp/2 
        return p0 - pdot*T_by_2 - pdd*0.5*T_by_2*T_by_2


def a_to_pdot(P_s, acc_ms2):
    return P_s * acc_ms2 /LIGHT_SPEED

def j_to_pdd(P_s, jerk_ms3):
    return 0.33*P_s * jerk_ms3 /LIGHT_SPEED

def middle_epoch(epoch_start, no_of_samples, tsamp):
     return epoch_start +0.5*no_of_samples*tsamp 

def extract_and_fold(opts):

    xml={}

    xml_file = opts.path+'/overview.xml'
    tree = ET.parse(xml_file)
    root = tree.getroot()

    #initiate empty arrays
    mod_period=[]
    period=[]
    acc=[]
    jerk=[]
    pdot=[]
    pdd=[]
    dm=[]
    snr=[]

    #datetime parameters
    xml['datetime'] = root.find('misc_info/utc_datetime').text.replace(":","-")
    #Header Parameters
    xml['ra'] = root.find('header_parameters/src_raj').text
    xml['dec'] = root.find('header_parameters/src_dej').text
    source_name = root.find('header_parameters/source_name').text
    xml['source_name'] = source_name.replace(" ","").replace(":","").replace(",","")
    #raw_data_filename=root.find('header_parameters/rawdatafile').text 
    xml['epoch_start'] = float(root.find("header_parameters/tstart").text)
    xml['tsamp'] = float(root.find("header_parameters/tsamp").text)
    xml['no_of_samples'] = int(root.find("header_parameters/nsamples").text)

    #Search Parameters
    xml['infile_name'] = root.find("search_parameters/infilename").text
    xml['fft_size'] = float(root.find('search_parameters/size').text)


    for P in root.findall("candidates/candidate/period"):
        period.append(float(P.text))
    for A in root.findall("candidates/candidate/acc"):
        acc.append(float(A.text))
    for J in root.findall("candidates/candidate/jerk"):
        jerk.append(float(J.text))
    for D in root.findall("candidates/candidate/dm"):
        dm.append(float(D.text))
    for s in root.findall("candidates/candidate/snr"):
        snr.append(float(s.text))
    
    for i in range(len(period)):
        Pdot = a_to_pdot(period[i],acc[i])
        Pdd = j_to_pdd(period[i],jerk[i])
        #mod_period.append(period_modified(period[i],Pdot,xml['no_of_samples'],xml['tsamp'],xml['fft_size']))
        mod_period.append(period_modified_new(period[i],Pdot,Pdd,xml['no_of_samples'],xml['tsamp'],xml['fft_size']))
        pdot.append(Pdot)
        pdd.append(Pdd)

   


     
    # Set all paths
    
    input_file = xml['infile_name']
    output_path = opts.outp
 
    #input_name = opts.fil_path+"/"+input_file
    input_name = opts.fil_path
    source_name = xml['source_name']    
    #mask_path = opts.mask



    # Run in batches
    folding_packet={}
    folding_packet['period'] = mod_period[opts.cand_no]
    folding_packet['acc'] = acc[opts.cand_no]
    folding_packet['jerk'] = jerk[opts.cand_no]
    folding_packet['pdot'] = pdot[opts.cand_no] 
    folding_packet['pdd'] = pdd[opts.cand_no] 
    folding_packet['dm'] = dm[opts.cand_no] 
    output_name="refold_%d_dm_%.2f_acc_%.2f_jerk_%.2f"%(opts.cand_no,folding_packet['dm'],folding_packet['acc'],folding_packet['jerk'])
    #output_name2="refold_%d_dm0_acc_%.2f"%(opts.cand_no,folding_packet['acc'])
    print output_name
    #subprocess.check_call("prepfold  -ncpus 12 %s -noxwin -nosearch -nodmsearch -topo -p %s -pd %s -dm 0.0 -o %s %s"%(opts.po,str(folding_packet['period']),str(folding_packet['pdot']),output_name2,input_name),shell=True,cwd=opts.outp)
    subprocess.check_call("prepfold  -ncpus 12 -n 128  -noxwin  -nodmsearch -topo -p %s -pd %s -pdd %s -dm %s -o %s %s"%(str(folding_packet['period']),str(folding_packet['pdot']),str(folding_packet['pdd']),str(folding_packet['dm']),output_name,input_name),shell=True,cwd=opts.outp)
 

                
if __name__=='__main__':


    # Update all input arguments
    parser = optparse.OptionParser()
    parser.add_option("-P",type=str,help="Path to XML file",dest="path")
    parser.add_option("-F",type=str,help="Path to filterbank file to fold",dest="fil_path")
    parser.add_option("-C",type=int,help="Candidate number",dest="cand_no")
    parser.add_option("--prepfold_options",type=str,help="Extra prepfold options",dest="po")
    parser.add_option("-O",type=str,help="Path to output files",dest="outp")
    #parser.add_option("-M",type=str,help="full path of mask file",dest="mask")
    parser.add_option("--log_level",dest='log_level',type=str,help="Logging level for PREPFOLD fold logger",default="INFO")
    #parser.add_option("--batch_number",dest='batch',type=int,help="number of PREPFOLD runs in one batch",default=24)
    opts,args = parser.parse_args()


    # Setup logging config
    log_type=opts.log_level
    log.setLevel(log_type.upper())

    #Extract and Fold
    extract_and_fold(opts)


