#!/usr/bin/env python

import struct
import os
import numpy as np
import itertools
import pylab as plt
import sys
import tkMessageBox
from Tkinter import Tk
from tkFileDialog import askopenfilename
import tkSimpleDialog
from Tkinter import *
import io
import platform

'''
Tool for converting Siemens TWIX format MRS data into dangerplot format for use with the TARQUIN analysis software (http://tarquin.sourceforge.net/).

Most of this code is a loose Python translation from the Philipp Ehses' MATLAB code mapVBVD.m found in the GANNET 2.0 source tree. Some inspiration was also taken from the twix_parser.py code by Brian Soher and Philip Semanchuk on behalf of the Vespa project. The Vespa project TWIX parser code also acknowledges Assaf Tal, Uzay Emir, and Jonathan Polimeni as providing code to inform thiers. Jamie Near also contributed MATLAB code used for reading parameters and raw data reconstruction.

Martin Wilson
April 2015
'''

def get_fd(td):
    return np.real(np.fft.fftshift(np.fft.fft(td)))

def bits(f):
    bytes = (ord(b) for b in f.read())
    for b in bytes:
        for i in xrange(8):
            yield (b >> i) & 1

def bin(a, pad):
    s=''
    t={'0':'000','1':'001','2':'010','3':'011',
       '4':'100','5':'101','6':'110','7':'111'}
    for c in oct(a)[1:]:
            s+=t[c]

    s = (pad-len(s))*'0' + s
    return s

def write_dpt(outputfile, N, fs, ft, phi0, phi1, ref, te, data):
    StrList = []
    StrList.append("Dangerplot_version\t2.0\n")
    StrList.append("Number_of_points\t%d\n" % N)
    StrList.append("Sampling_frequency\t%8.8e\n" % fs)
    StrList.append("Transmitter_frequency\t%8.8e\n" % ft)
    StrList.append("Phi0\t%8.8e\n" % phi0)
    StrList.append("Phi1\t%8.8e\n" % phi1)
    StrList.append("PPM_reference\t%8.8e\n" % ref)
    StrList.append("Echo_time\t%8.8e\n" % te)
    StrList.append("Real_FID\tImag_FID\t\n")
    for x in data:
        StrList.append("%8.8e %8.8e\n" % (x.real, x.imag))

    output = open(outputfile, 'w')
    output.writelines(StrList)
    output.close

def read_uint16(data):
    return struct.unpack('<H', data.read(2))[0]

def read_uint32(data):
    return struct.unpack('<I', data.read(4))[0]

def read_uint64(data):
    return struct.unpack('<Q', data.read(8))[0]

def read_ushort(data):
    return struct.unpack('<H', data.read(2))[0]

def read_byte(data):
    return struct.unpack('<B', data.read(1))[0]

def read_int(data):
    return struct.unpack('<i', data.read(4))[0]

def read_uint(data):
    return struct.unpack('<I', data.read(4))[0]

def read_float(data):
    return struct.unpack('<f', data.read(4))[0]

def convert_mrs_data(fname, block_size = 8):
    twix_data = open(fname, "rb")
    firstInt = read_uint32(twix_data)
    secondInt = read_uint32(twix_data)

    # check for VB or VD foramt
    if ( (firstInt < 10000) and (secondInt <= 64) ):
        print 'File is VD format.'
        version = 'vd'
        NScans = secondInt
        measID = read_uint32(twix_data)
        fileID = read_uint32(twix_data)
        measOffset = read_uint64(twix_data) # usually 10240 bytes
        measLength = read_uint64(twix_data)
        twix_data.seek(measOffset, 0)
        hdrLength = read_uint32(twix_data)
        dataStart = measOffset + hdrLength
    else:
        print 'File is VB format.'
        version = 'vb'
        dataStart = firstInt
        NScans = 1

    twix_data.close()
    #print dataStart

    if platform.system() == 'Windows':
        twix_data_txt = io.open(fname, "r",newline='')
    else:
        twix_data_txt = open(fname, "r")

    array = []
    for line in twix_data_txt:
        if line.startswith("ulVersion") and (version == "vb"):
            break

        if line.startswith("ulVersion") and (version == "vd"):
            tSequenceFilename = twix_data_txt.next()
            tProtocolName = twix_data_txt.next()
            if tProtocolName != 'tProtocolName\t = \t"AdjCoilSens"\n' and tProtocolName != 'tProtocolName\t = \t"CBU_MPRAGE_32chn"\n':
                print tProtocolName
                break

    for line in twix_data_txt:
        if "### ASCCONV END ###" in line:
            break

        elif line.startswith("lAverages"):
            averages = int(line.split("=")[1])

        elif line.startswith("sRXSPEC.alDwellTime[0]"):
            fs = 1.0e9/(float(line.split("=")[1]))

        elif line.startswith("sTXSPEC.asNucleusInfo[0].lFrequency"):
            ft = float(line.split("=")[1])

        elif line.startswith("alTE[0]"):
            te = float(line.split("=")[1])/1.0e6

        elif line.startswith("sSpecPara.lVectorSize"):
            N = int(line.split("=")[1])

    twix_data_txt.close()

    twix_data = open(fname, "rb")
    twix_data.seek(0, os.SEEK_END)
    file_size = twix_data.tell()
    data_bytes = file_size - dataStart
    twix_data.seek(dataStart)

    last_scan = False
    scans = 0
    ima_echoes = 0
    fid_list = []
    cPos = dataStart

    for scans in range(NScans):
        n = 0
        while last_scan == False:
            twix_data.seek(cPos, 0)
            print str(round(float(cPos) / float(file_size) * 100.0,1)) + "% complete."
            #print "Scan : " + str(scans)
            n = n + 1
            #print "N = " + str(n)
            ulDMALength_bin = bin(read_uint32(twix_data),32)
            ulDMALength = int(ulDMALength_bin[-25:],2)
            #print "ulDMALength = " + str(ulDMALength)
            #print "cPos = " + str(cPos)
            twix_data.seek(cPos, 0)

            if version == "vb":
                twix_data.seek(20, 1) # move ahead 20 bytes
            else:
                twix_data.seek(40, 1) # move ahead 40 bytes

            eval_info_mask     = [ bin(read_byte(twix_data),8)[::-1] for i in range(8) ]
            info_bits = eval_info_mask[0]+eval_info_mask[1]+eval_info_mask[2]+eval_info_mask[3]+eval_info_mask[4]+eval_info_mask[5]+eval_info_mask[6]+eval_info_mask[7]
            #print info_bits
            samples_in_scan    = read_ushort(twix_data)
            used_channels      = read_ushort(twix_data)
            #print "Samples  : " + str(samples_in_scan)
            #print "Channels : " + str(used_channels)
            twix_data.seek(8, 1)
            Necho                 = read_ushort(twix_data)
            twix_data.seek(22, 1)
            kspace_center_column  = read_ushort(twix_data)

            MDH_ACQEND            = info_bits[0] == "1"
            MDH_RTFEEDBACK        = info_bits[1] == "1"
            MDH_HPFEEDBACK        = info_bits[2] == "1"
            MDH_SYNCDATA          = info_bits[5] == "1"
            MDH_RAWDATACORRECTION = info_bits[10] == "1"
            MDH_REFPHASESTABSCAN  = info_bits[14] == "1"
            MDH_PHASESTABSCAN     = info_bits[15] == "1"
            MDH_SIGNREV           = info_bits[17] == "1"
            MDH_PHASCOR           = info_bits[21] == "1"
            MDH_PATREFSCAN        = info_bits[22] == "1"
            MDH_PATREFANDIMASCAN  = info_bits[23] == "1"
            MDH_REFLECT           = info_bits[24] == "1"
            MDH_NOISEADJSCAN      = info_bits[25] == "1"
            MDH_IMASCAN           = True

            if ( MDH_ACQEND or MDH_RTFEEDBACK or MDH_HPFEEDBACK or MDH_PHASCOR  or MDH_NOISEADJSCAN or MDH_SYNCDATA ):
                MDH_IMASCAN = False

            if ( MDH_PHASESTABSCAN or MDH_REFPHASESTABSCAN ):
                MDH_PATREFSCAN = False;
                MDH_PATREFANDIMASCAN = False;
                MDH_IMASCAN = False; 
            
            if ( MDH_PATREFSCAN and not(MDH_PATREFANDIMASCAN) ):
                MDH_IMASCAN = False

            if version == "vb":
                if ( not(MDH_SYNCDATA) and not(MDH_ACQEND) ):
                    ulDMALength = (2*4*samples_in_scan + 128) * used_channels;

            elif version == "vd":
                if ( not(MDH_SYNCDATA) and not(MDH_ACQEND) and ulDMALength != 0 ):
                    ulDMALength = 192 + (2*4*samples_in_scan + 32) * used_channels;
            
            if MDH_IMASCAN:
                #print "Found IMA"
                ima_coils = used_channels
                ima_samples = samples_in_scan
                if Necho > ima_echoes:
                    ima_echoes = Necho
                ima_kspace_center_column = kspace_center_column 
                # Lets read this average in
                for x in range(used_channels):
                    if version == "vb":
                        twix_data.seek(128 + cPos + x * (2*4*samples_in_scan+128), 0)
                        raw = np.fromfile(twix_data, dtype=np.dtype('<f'), count = samples_in_scan*2)
                        data_iter = iter(raw)
                        raw_cplx = [complex(r, i) for r, i in itertools.izip(data_iter, data_iter)]
                        fid_list.append(np.conj(np.complex64(raw_cplx)))

                    elif version == "vd":
                        twix_data.seek(192 + 32 + cPos + x * (2*4*samples_in_scan+32), 0)
                        raw = np.fromfile(twix_data, dtype=np.dtype('<f'), count = samples_in_scan*2)
                        data_iter = iter(raw)
                        raw_cplx = [complex(r, i) for r, i in itertools.izip(data_iter, data_iter)]
                        fid_list.append(np.conj(np.complex64(raw_cplx)))

            if ( MDH_ACQEND or ulDMALength == 0 ): # break out to the next scan
                if (scans < NScans-1):
                    cPos = cPos + ulDMALength
                    cPos = cPos + 512 - (cPos % 512)
                    twix_data.seek(cPos)
                    hdrLength  = read_uint32(twix_data);
                    cPos = cPos + hdrLength
                    twix_data.seek(cPos, 0)

                twix_data.seek(cPos, 0)
                break

            if MDH_SYNCDATA:
                cPos = cPos + ulDMALength
                twix_data.seek(cPos, 0)
                continue

            cPos = cPos + ulDMALength
            twix_data.seek(cPos, 0)

            if twix_data.tell() >= file_size:
                print "I don't belong here."
                break
            
    data_avs = len(fid_list)/ima_coils

    print str(data_avs) + " averages."
    print str(ima_coils) + " coils."
    print str(ima_kspace_center_column/2+1) + " points to FID offset."
    print str(len(fid_list[0])) + " complex data points."

    fid_np = np.reshape(fid_list,[data_avs,ima_coils,len(fid_list[0])])
    del fid_list

    if (ima_kspace_center_column/2 > 0):
        fid_np = fid_np[:,:,ima_kspace_center_column/2+1:]

    '''
    for x in range(32):
        plt.plot(get_fd(fid_np_rep_av[x]))

    plt.show()
    '''

    # sum across the first point in the repeats for each coil
    #av_first_point_coil = np.angle(fid_np[::2,:,0].sum(0))
    av_first_point_coil = fid_np[:,:,0].sum(0)

    # take the phase and amplitude
    phase = np.angle(av_first_point_coil)
    amp = np.abs(av_first_point_coil)

    fid_np_corr = fid_np.copy()
    # phase and scale data from each coil
    for n in range(len(phase)):
        fid_np_corr[:,n,:] = fid_np_corr[:,n,:] * np.exp(-1.0j*phase[n]) * amp[n]

    # sum across the coils
    fid_np_coil_av = fid_np_corr.sum(1)

    # get odd scans
    fid_np_rep_coil_av_odd = fid_np_coil_av[::2].sum(0)

    # get even scans
    fid_np_rep_coil_av_even = fid_np_coil_av[1::2].sum(0)

    # sum across the coils and repeats
    fid_np_rep_coil_av = fid_np_coil_av.sum(0)

    '''
    for x in range(2):
        plt.plot(get_fd(fid_np_coil_av[x]))

    plt.show()
    '''


    blocks = block_size
    if (ima_echoes == 1):
        n = np.shape(fid_np_coil_av[::2])[0]
        new_n = n/blocks
        odd_blocks = np.zeros([new_n, np.shape(fid_np_coil_av[::2])[1]],np.complex64)
        even_blocks = np.zeros([new_n, np.shape(fid_np_coil_av[::2])[1]],np.complex64)
        inter_blocks = np.zeros([(new_n*2), np.shape(fid_np_coil_av[::2])[1]],np.complex64)
        for x in range(new_n):
            ind1 = blocks*x
            ind2 = blocks*x+blocks
            odd_blocks[x,:] = np.average(fid_np_coil_av[::2][ind1:ind2,:],0)
            even_blocks[x,:] = np.average(fid_np_coil_av[1::2][ind1:ind2,:],0)
            inter_blocks[(2*x),:] = np.average(fid_np_coil_av[::2][ind1:ind2,:],0)
            inter_blocks[(2*x+1),:] = np.average(fid_np_coil_av[1::2][ind1:ind2,:],0)
    else:
        n = np.shape(fid_np_coil_av)[0]
        new_n = n/blocks
        inter_blocks = np.zeros([new_n, np.shape(fid_np_coil_av)[1]],np.complex64)
        for x in range(new_n):
            ind1 = blocks*x
            ind2 = blocks*x+blocks
            inter_blocks[x,:] = np.average(fid_np_coil_av[ind1:ind2,:],0)

    #plt.plot(get_fd(fid_np_rep_coil_av_odd))
    #plt.plot(get_fd(fid_np_rep_coil_av_even))
    #plt.show()

    #plt.plot(get_fd(fid_np_rep_coil_av_odd-fid_np_rep_coil_av_even))
    #plt.show()

    #spec = get_fd(fid_np_rep_coil_av)
    #plt.plot(spec/np.max(spec))
    #plt.show()
    
    # get the full path to data
    full_path = os.path.abspath(fname)
    dpt_path = full_path[:-3] + "dpt"
    print "Writing to : " + dpt_path

    # write to file
    write_dpt(dpt_path, np.shape(fid_np)[2], fs, ft, 0, 0, 4.65, te, inter_blocks.flatten())
    print "Finished"


class MyDialog(tkSimpleDialog.Dialog):
    def body(self, master):
        Label(master, text="Block size:").grid(row=0)
        #Label(master, text="Block size:").grid(row=1)
        self.e1 = Entry(master)
        self.e1.insert(0,8)
        #self.e2 = Entry(master)
        #self.e2.insert(0,8)
        self.e1.grid(row=0, column=1)
        #self.e2.grid(row=1, column=1)
        #return self.e1 # initial focus

    def apply(self):
        first = int(self.e1.get())
        #second = int(self.e2.get())
        self.result = first #, second
        #print first, second # or something


Tk().withdraw() 
filename = askopenfilename(title="Choose MRS TWIX file.")

if filename == "":
    print "Exiting"
    sys.exit()

print(filename)
root = Tk()
root.withdraw() 
d = MyDialog(root,title="Data options")

if d.result == None:
    print "Exiting"
    sys.exit()

convert_mrs_data(filename, d.result)
tkMessageBox.showinfo(title="Update", message="Conversion complete")
