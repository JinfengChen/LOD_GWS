#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import re
import os
import argparse
import matplotlib as mp
mp.use('pdf')
import matplotlib.pyplot as plt
import pandas as pd
from urllib2 import urlopen
import numpy as np
import pylab as pl
import matplotlib.colors as colors
import random
from matplotlib import gridspec

def usage():
    test="name"
    message='''
python LOD_GWAS.py --input QTL.regions.GWAS.list

    
'''
    print message

def set_ticks_X(ax, ylen, ylab):
    fig = plt.gcf()
    #fig.set_size_inches(8, 8)

    # turn off the frame
    ax.set_frame_on(False)

    # put the major ticks at the middle of each cell
    #ax.set_yticks(np.arange(qtl.shape[0]) + 0.5, minor=False)
    ax.set_xticks(np.arange(ylen) + 0.5, minor=False)

    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    # Set the labels
    ax.set_xticklabels(ylab, minor=False, fontsize=8)
    ax.set_yticklabels([])

    # rotate the
    plt.xticks(rotation=40)

    ax.grid(False)

    # Turn off all the ticks
    ax = plt.gca()

    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False

    return ax

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def set_ticks_XY(ax, xlen, ylen, xlab, ylab):
    fig = plt.gcf()
    #fig.set_size_inches(8, 11)

    # turn off the frame
    ax.set_frame_on(False)

    # put the major ticks at the middle of each cell
    ax.set_yticks(np.arange(xlen) + 0.5, minor=False)
    ax.set_xticks(np.arange(ylen) + 0.5, minor=False)

    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    # Set the labels
    ax.set_xticklabels(xlab, minor=False)
    ax.set_yticklabels(ylab, minor=False)

    # rotate the
    plt.xticks(rotation=90)

    ax.grid(False)

    # Turn off all the ticks
    ax = plt.gca()

    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False

    return ax

def set_ticks_XY_Right(ax, xlen, ylen, xlab, ylab):
    fig = plt.gcf()
    #fig.set_size_inches(8, 11)

    # turn off the frame
    ax.set_frame_on(False)
    ax.yaxis.set_ticks_position('right')

    # put the major ticks at the middle of each cell
    ax.set_yticks(np.arange(xlen) + 0.5, minor=False)
    ax.set_xticks(np.arange(ylen) + 0.5, minor=False)

    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    # Set the labels
    ax.set_xticklabels(xlab, minor=False, fontsize=8)
    ax.set_yticklabels(ylab, minor=False)

    # rotate the
    plt.xticks(rotation=40)

    ax.grid(False)

    # Turn off all the ticks
    ax = plt.gca()

    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False

    return ax



def colorbar(heatmap,  cbaxes, label, ticks, ticke):
    # Add a colorbar below the heatmap
    # add_axes refer to [left, bottom, width, height], where the coordinates are just fractions that go from 0 to 1 of the plotting area.
    #cbaxes = fig.add_axes([0.15, 0.15, 0.2, 0.02]) 
    cb = plt.colorbar(heatmap, cax = cbaxes, orientation='horizontal', ticks=np.linspace(ticks, ticke, 5)) 
    #cb = pl.colorbar(heatmap, orientation='horizontal', shrink=0.5825,
    #                 fraction=0.02, pad=-0.035,ticks=np.linspace(-1, 1, 5),
    #                 use_gridspec=True)
    cb.set_label('$\mathrm{%s}\ $' %(label))
    #return fig

def heatmap(gene_list, rnaseq1_list, rnaseq2_list, rnaseq3_list, annotation, project):
    # Create a figure.
    figsize=(8,11)
    fig = pl.figure(figsize=figsize)
    pl.subplots_adjust(bottom=0.2)
    gs = gridspec.GridSpec(1, 4, width_ratios=[1, 3, 3, 3])

    # Read Dataframe from file
    qtl = pd.read_table(gene_list, index_col=0)
    rnaseq1 = pd.read_table(rnaseq1_list, index_col=0)
    rnaseq2 = pd.read_table(rnaseq2_list, index_col=0)
    rnaseq3 = pd.read_table(rnaseq3_list, index_col=0)
    anno = pd.read_table(annotation, sep='\t')

    #use truncated color table from 0.2 to 0.8
    cmap = plt.cm.get_cmap("Blues")
    new_cmap = truncate_colormap(cmap, 0.2, 0.8)
 

    # Draw heatmap for variations
    #ax = fig.add_subplot(131, frame_on=False)
    #ax = fig.add_subplot(133, frame_on=True)
    ax0 = plt.subplot(gs[0])
    heatmap = ax0.pcolor(qtl, cmap=new_cmap, alpha=0.8, vmin= 0, vmax= 4)
    # Set names for x,y_ticks
    ax0 = set_ticks_XY(ax0, qtl.shape[0], qtl.shape[1], qtl.columns, qtl.index)
    
    ticks0 = 0
    ticke0 = 4 
    cbaxes0 = fig.add_axes([0.08, 0.18, 0.2, 0.02])
    colorbar(heatmap, cbaxes0, 'Variation Effect', ticks0, ticke0)
    
    # Draw heatmap for gene expression: cold
    #ax = fig.add_subplot(132, frame_on=False) 
    ax1 = plt.subplot(gs[1])
    heatmap = ax1.pcolor(rnaseq1, cmap=plt.cm.cool, vmin= -2, vmax= 2, alpha=0.8)  
    plt.title('Cold', y=1.065) 
    # Set names for x_ticks only
    ax1 = set_ticks_X(ax1, rnaseq1.shape[1], rnaseq1.columns)
    #ax = fig.add_subplot(133, frame_on=True)
    
    ticks1 = -2
    ticke1 = 2
    cbaxes1 = fig.add_axes([0.47, 0.18, 0.2, 0.02])
    colorbar(heatmap, cbaxes1, 'Log2(HEG4/NB)', ticks1, ticke1)

   # Draw heatmap for gene expression: drought
    #ax = fig.add_subplot(132, frame_on=False) 
    ax2 = plt.subplot(gs[2])
    heatmap = ax2.pcolor(rnaseq2, cmap=plt.cm.cool, vmin= -2, vmax= 2, alpha=0.8)
    plt.title('Drought', y=1.065) 
    # Set names for x_ticks only
    ax2 = set_ticks_X(ax2, rnaseq2.shape[1],  rnaseq2.columns)

    # Draw heatmap for gene expression: salt
    #ax = fig.add_subplot(132, frame_on=False) 
    ax3 = plt.subplot(gs[3])
    heatmap = ax3.pcolor(rnaseq3, cmap=plt.cm.cool, vmin= -2, vmax= 2, alpha=0.8)
    plt.title('Salt', y=1.065) 
    # Set names for x_ticks only
    ylabs = anno['Annotation']
    ax3 = set_ticks_XY_Right(ax3, rnaseq3.shape[0], rnaseq3.shape[1], rnaseq3.columns, ylabs)
    
 
    # Draw LOD curve
    #ax2 = plt.subplot(gs[2])
    #x=random.sample(range(1,20),6)
    #y=np.arange(6)
    #plt.xlim=(10,20)
    #ax2.plot(x,y)
    # Set names for x_ticks 
    #ax2 = set_ticks_X(ax2)

    # Save file
    fig.savefig('%s.pdf' %(project), bbox_inches='tight')

'''

input:
Chr	Start	End	Name	Length	QTL_TRAIT	GWAS:TRAIT
Chr3	1255241	1479666	HeadingDate1	224425	HeadingDays	Naive:FloweringtimeatArberdeen,Naive:FloweringtimeatArkansas

lodfile:
Position        Bin_start       Bin_end chr     pos     HeadingDays     Plant.Height..cm..in.Field      Biomass Number.of.Tillers       Single.Plant..Grain.Yield
111023  0       222046  Chr1       0       0.00726849104738392     1.75434974017193        0.706883828535311       0.745160416472634       1.50828253471261        1.6934514845984
361453  222046  500860  Chr1       0.38023 0.00249508279722122     2.48456985459131        0.606550483403985       0.662605759329817       1.29921364448124        1.99225214187723

GWAS:
CHR     BP      SNP     P
Chr1    14147   id1000001       0.312923632962147
Chr1    74192   id1000003       0.288067975202483
Chr1    75969   id1000005       0.116683519464606
Chr1    76852   id1000007       0.665625290523823

'''
def subfiles(infile, lodfile, gwas_prefix):
    flank = 1000000
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if line.startswith('Chr'): 
                unit = re.split(r'\t',line)
                chrs = unit[0]
                start= int(unit[1]) - flank if int(unit[1]) - flank > 0 else 0
                end  = int(unit[2]) + flank
                qtl_trait = unit[5]
                gwas_trait = unit[6]
                prefix = '%s' %(unit[3])
                #LODfile 
                lodfile, start, end = sublodfile(lodfile, chrs, start, end, qtl_trait, prefix)
                #GWASfile
                gwasfiles= subgwasfile(gwas_prefix, chrs, start, end, gwas_trait, prefix)

def sublodfile(lodfile, chrs, start, end, qtl_trait, prefix):
    outfile = prefix + '.LOD.subfile'
    ofile = open(outfile, 'w')
    print >> ofile, 'Position\tBin_start\tBin_end\tChr\tLOD'
    lod_data = pd.read_table(lodfile)
    position = []
    for i in range(len(lod_data['Position'])):
        #print lod_data['chr'][i], lod_data['Position'][i], chrs, start, end
        if lod_data['chr'][i] == chrs and int(lod_data['Position'][i]) > start and int(lod_data['Position'][i]) < end:
            line = '%s\t%s\t%s\t%s\t%s' %(lod_data['Position'][i], lod_data['Bin_start'][i], lod_data['Bin_end'][i], lod_data['chr'][i], lod_data[qtl_trait][i])
            position.append(lod_data['Bin_start'][i])
            position.append(lod_data['Bin_end'][i])
            print >> ofile, line
    ofile.close()
    position = sorted(position, key=int)
    start = position[0]
    end   = position[-1]
    return outfile, start, end

def subgwasfile(gwas_prefix, chrs, start, end, gwas_trait, prefix):
    gwasfiles = []
    traits = re.split(r',', gwas_trait)
    for t in traits:
        pop, trait = re.split(r':', t)
        #print pop, trait
        gwasfile = '%s%s/%s.gwas' %(gwas_prefix, pop, trait)
        outfile  = '%s.%s.%s.GWAS.subfile' %(prefix, pop, trait)
        gwasfiles.append(outfile)
        ofile = open(outfile, 'w')
        print >> ofile, 'CHR\tBP\tSNP\tP'
        gwas_data = pd.read_table(gwasfile)
        for i in range(len(gwas_data['CHR'])):
            #print lod_data['chr'][i], lod_data['Position'][i], chrs, start, end
            if gwas_data['CHR'][i] == chrs and int(gwas_data['BP'][i]) > start and int(gwas_data['BP'][i]) < end:
                line = '%s\t%s\t%s\t%s' %(gwas_data['CHR'][i], gwas_data['BP'][i], gwas_data['SNP'][i], gwas_data['P'][i])
                print >> ofile, line
        ofile.close()
    return gwasfiles


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    lodfile = '../input/MPR.cross.uniq.QTL.mr.table.new'
    gwas_prefix = '../input/Rice.44k_' 
    subfiles(args.input, lodfile, gwas_prefix)
    

if __name__ == '__main__':
    main()

