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
from matplotlib.patches import Rectangle
import math
import matplotlib.lines as mlines

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


def set_ticks_XY_Left_LOD(ax, xmin, xmax, ymin, ymax):
    fig = plt.gcf()

    # turn off the frame
    ax.set_frame_on(False)
    ax.yaxis.set_ticks_position('left')

    # Set the labels
    ax.set_xticklabels([], minor=False)

    ax.grid(False)

    # Turn off all the ticks
    ax = plt.gca()

    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False

    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    
    #reset y tick location and labels 
    yloc =  map(int, ax.yaxis.get_majorticklocs())
    ax.set_yticks(yloc, minor=False)
    ax.set_yticklabels(yloc, minor=False)
   
    ymax1 = max(yloc)

    #draw x and y axix line  
    ly = mlines.Line2D([xmin, xmin], [ymin, ymax1], color='black')
    #lx = mlines.Line2D([xmin, xmax], [ymin, ymin], color='black')
    ax.add_line(ly)
    #ax.add_line(lx)

    ax.set_ylabel('LOD')
    ax.yaxis.set_label_coords(-0.03, 0.5)

    return ax


def set_ticks_XY_Left_Bottom(ax, xmin, xmax, ymin, ymax, xstart, chrs):
    fig = plt.gcf()

    # turn off the frame
    ax.set_frame_on(False)
    ax.yaxis.set_ticks_position('left')

    # Set the labels
    #ax.set_xticklabels([], minor=False)
    ax.xaxis.set_ticks_position('bottom')
    #xlabeltext = [item.get_text() for item in ax.get_xticklabels()]
    #print xlabeltext[1]
    #ax.set_xticklabels(xlabels, minor=False)
    xlabels = []
    xpos    = []
    for loc in ax.xaxis.get_majorticklocs():
        xpos.append(loc)
        xlab = (int(loc) + int(xstart))/1000
        xlabels.append(xlab)
    print xpos
    print xlabels
    ax.set_xticks(xpos, minor=False)
    ax.set_xticklabels(xlabels, minor=False)
    
    ax.set_xlabel('Chromosome %s (kb)' %(chrs))
    ax.xaxis.set_label_coords(0.5, -0.3) 

    ax.grid(False)

    # Turn off all the ticks
    ax = plt.gca()

    #for t in ax.xaxis.get_major_ticks():
        #t.tick1On = False
        #t.tick2On = False

    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)

    #reset x tick location and labels
    #xloc =  map(int, ax.xaxis.get_majorticklocs())
    #ax.set_xticks(xloc, minor=False)
    #ax.set_xticklabels(xloc, minor=False)
   
 
    #reset y tick location and labels 
    yloc =  map(int, ax.yaxis.get_majorticklocs())
    ax.set_yticks(yloc, minor=False)
    ax.set_yticklabels(yloc, minor=False)
   
    ymax1 = max(yloc)

    #draw x and y axis line  
    ly = mlines.Line2D([xmin, xmin], [ymin, ymax1], color='black')
    lx = mlines.Line2D([xmin, xmax], [ymin, ymin], color='black')
    ax.add_line(ly)
    ax.add_line(lx)

    #set y axis annotation
    ax.set_ylabel(r'$-log_{10}(P)$')
    ax.yaxis.set_label_coords(-0.03, 0.5)

    #reset x tick location and labels
    #xloc =  map(int, ax.xaxis.get_majorticklocs())
    #ax.set_xticks(xloc, minor=False)
    #ax.set_xticklabels(xloc, minor=False)
    
    #ax.set_xticklabels([], minor=False)

    return ax


def set_ticks_XY_Left(ax, xmin, xmax, ymin, ymax):
    fig = plt.gcf()

    # turn off the frame
    ax.set_frame_on(False)
    ax.yaxis.set_ticks_position('left')

    # Set the labels
    ax.set_xticklabels([], minor=False)

    ax.grid(False)

    # Turn off all the ticks
    ax = plt.gca()

    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False

    #plt.xlim(xmin, xmax)
    #plt.ylim(ymin, ymax)
    
    #reset y tick location and labels 
    yloc =  map(int, ax.yaxis.get_majorticklocs())
    ax.set_yticks(yloc, minor=False)
    ax.set_yticklabels(yloc, minor=False)
   
    ymax1 = max(yloc)

    #draw x and y axis line  
    ly = mlines.Line2D([xmin, xmin], [ymin, ymax1], color='black')
    lx = mlines.Line2D([xmin, xmax], [ymin, ymin], color='black')
    ax.add_line(ly)
    ax.add_line(lx)

    #set y axis annotation
    ax.set_ylabel(r'$-log_{10}(P)$')
    ax.yaxis.set_label_coords(-0.03, 0.5)
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
def subfiles(infile, lodfile, genefile, gwas_prefix):
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
                lodfiles, start, end = sublodfile(lodfile, chrs, start, end, qtl_trait, prefix)
                #genefile
                genefiles = subgenefile(genefile, chrs, start, end, prefix)
                #GWASfile
                gwasfiles = subgwasfile(gwas_prefix, chrs, start, end, gwas_trait, prefix)
                #Plot
                plot_region(lodfiles, gwasfiles, genefiles, start, end, prefix)

def plot_region(lodfiles, gwasfiles, genefiles, start, end, prefix):
    nrow = len(lodfiles) + len(gwasfiles)
    xmin = 0
    xmax = end - start +1 
    # Create a figure.
    figsize=(12,6)
    fig = plt.figure(figsize=figsize)
    plt.subplots_adjust(bottom=0.2)
  
    # split plot into 4 rows X 1 cols. gs[1, 0] is 2th row and 1th column
    # gs[0:, 0]  is all rows of 1th column
    print 'number of row: %s' %(nrow)
    gs = gridspec.GridSpec(int(nrow), 1, hspace=0.2)

    #plot qtl lod curve
    #Position        Bin_start       Bin_end Chr     LOD
    #343751  290088  397414  Chr3    3.09188403328
    rank = 0
    axs = []
    for i in range(len(lodfiles)):
        print i
        lod_data = pd.read_table(lodfiles[i])
        positions= map(lambda n:n-start, lod_data['Position'])
        lod_score= lod_data['LOD']
        ax0 = fig.add_subplot(gs[i])
        ax0.plot(positions, lod_score, marker='.', color='lightblue', lw=1)
        rank += 1
        ymin = int(min(lod_score)*0.95)
        ymax = max(lod_score)*1.02
        plt.ylim(ymin, ymax)
        plt.xlim(xmin, xmax)
        # set axis
        ax0 = set_ticks_XY_Left_LOD(ax0, xmin, xmax, ymin, ymax)

        for j in range(len(lod_data['Bin_start'])):   
            bin_s = lod_data['Bin_start'][j] - start
            bin_e = lod_data['Bin_end'][j] - start
            width= bin_e - bin_s + 1
            height= (ymax-ymin)*0.1
            rx = bin_s
            ry = ymin
            #print j, rx, ry, width, height
            color = 'grey' if j % 2 == 0 else 'black'
            ax0.add_patch(Rectangle((rx, ry), width, height, color=color))
        axs.append(ax0)
    #plot gwas curve
    #CHR     BP      SNP     P
    #Chr3    307962  id3000025       0.391076770445
    #Chr3    309403  id3000030       0.131840540423
    for i in range(len(gwasfiles)):
        gwas_data = pd.read_table(gwasfiles[i])
        positions = map(lambda n:n-start, gwas_data['BP'])
        #pvalue    = map(math.log10, gwas_data['P'])
        pvalue    = map(lambda n:-n, map(math.log10, gwas_data['P']))
        ymax = max(pvalue) * 1.1
        print '%s' %(i+rank) 
        ax0 = fig.add_subplot(gs[i+rank])
        #ax0.plot(positions, pvalue, 'bo')
        ax0.scatter(positions, pvalue, c='lightskyblue', edgecolors='none')
        plt.xlim(xmin, xmax)
        plt.ylim(0, ymax)
        chrs = gwas_data['CHR'][1][3:]
        if i < len(gwasfiles)-1:
            ax0 = set_ticks_XY_Left(ax0, xmin, xmax, 0, ymax)
        else:
            ax0 = set_ticks_XY_Left_Bottom(ax0, xmin, xmax, 0, ymax, start, chrs) 

        axs.append(ax0)

    #plot gene line
    #Chr3    MSU7    gene    717447  720837  .       +       .       ID=Os03g0112700;MSU=LOC_Os03g02160;Name=Ehd4;Gene=early heading date 4;
    #Chr3    MSU7    gene    1269856 1300384 .       -       .       ID=Os03g0122600;MSU=LOC_Os03g03070,LOC_Os03g03100;Name=OsMADS50;Gene=osmads50;
    rx = re.compile(r'Name=(.*?);') 
    for i in range(len(genefiles)):
        gene_data = pd.read_table(genefiles[i], header=None)
        for r in range(len(gene_data[0])):
            gstart = int(gene_data[3][r]) - start
            gend   = int(gene_data[4][r]) - start
            gmidpoint = gstart + (gend - gstart)/2
            plt.axvline(x=gmidpoint, ymin=0, ymax=3.5, clip_on=False, linestyle='dashed', color='black')
            m = rx.search(gene_data[8][r])
            name = m.groups(0)[0] if m else 'None'
            plt.text(gmidpoint*0.98, 70, name, style='italic', horizontalalignment='left', verticalalignment='bottom', rotation=50)
            #for n in range(len(axs)):
                #axs[n].axvline(x=gmidpoint, ymin=-2, ymax=2, clip_on=False, color='black') 
    #save file
    fig.savefig('%s.pdf' %(prefix), bbox_inches='tight')

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
    return [outfile], start, end

def subgwasfile(gwas_prefix, chrs, start, end, gwas_trait, prefix):
    gwasfiles = []
    traits = re.split(r',', gwas_trait)
    for t in traits:
        pop, trait = re.split(r':', t)
        #print pop, trait
        gwasfile = '%s%s/%s.gwas' %(gwas_prefix, pop, trait)
        outfile  = '%s.%s.%s.GWAS.subfile' %(prefix, pop, trait)
        gwasfiles.append(outfile)
        #we will not rewrite file if file exists
        if os.path.isfile(outfile):
            continue
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

'''
Chr1    MSU7    gene    134291  135685  .       +       .       ID=Os01g0102300;MSU=LOC_Os01g01280;Name=OsTLP27;Gene=Thylakoid lumen
'''
def subgenefile(genefile, chrs, start, end, prefix):
    outfile = prefix + '.Gene.subfile'
    if os.path.isfile(outfile):
        return [outfile]
    ofile = open(outfile, 'w')
    with open (genefile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if unit[0] == chrs and int(unit[3]) > start and int(unit[3]) < end:
                    print >> ofile, line
    ofile.close()
    return [outfile]


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
    genefile= '../input/QTARO.rice.gene.gff'
    gwas_prefix = '../input/Rice.44k_' 
    subfiles(args.input, lodfile, genefile, gwas_prefix)
    

if __name__ == '__main__':
    main()

