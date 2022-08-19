#!/usr/bin/env python
# -*- coding: utf-8 -*-
# By Elena

from numpy import *
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.patches as mpatches
import matplotlib.colors as mc
import seaborn as sns
import colorsys
import datatable as dt
import math

# For VSCode + cluster: export PATH=/home/ebernab3/anaconda3/bin:$PATH

# Aesthetic stuff
##########################################

font = "Galvji"
path = '/Cluster_Filespace/Marioni_Group/Elena/other/fonts/Galvji/Galvji-01.ttf'
fe = matplotlib.font_manager.FontEntry(
    fname= path,
    name=font)
matplotlib.font_manager.fontManager.ttflist.insert(0, fe)
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.family'] = font


# CpG trajectory plot function
##########################################

def lighten_color(color, amount=1.2):  
    # --------------------- SOURCE: @IanHincks ---------------------
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])

def cpg_trajectory(df, cpg, gene_name, output_dir, cpg_col = "CpG", age_col = "Age", color = "#07B1D8", wave = "w3"):
    # Plot - all together
    fig, axes = plt.subplots(1, 2, figsize = (6.5, 3.5), sharey=True, gridspec_kw=dict(width_ratios=[4,1]))
    sns.regplot(ax = axes[0], data = df, x = age_col, y = cpg_col, color = color, marker = '.', scatter_kws={'alpha':1, 'rasterized':True}, line_kws={'color':lighten_color(color)}, lowess = True)
    sns.kdeplot(ax = axes[1], data = df, y = cpg_col, shade = True, color = color, linewidth = 2)
    if "beta" not in wave:
        axes[0].set_ylabel("CpG M-Value")
    else:
        axes[0].set_ylabel("CpG Beta-Value")
        axes[0].set_ylim([0,1])
    axes[0].set_xlabel("Age")
    plt.suptitle("%s | %s" % (cpg, gene_name))
    plt.tight_layout()
    plt.savefig(output_dir + "%s_%s_trajectory_%s.pdf" % (cpg, gene_name, wave), dpi = 300)
    plt.close(fig)

def genomic_pos(pos, col_chr, col_bp):
    l = pos.groupby(col_chr)[col_bp].max()
    offset = l.cumsum() - l
    return pos[col_bp] + list(offset.loc[pos[col_chr]])

def chromosome_colors(N, name='custom2'):
    colors = []
    if name == "default":
        colors = sns.color_palette("Set2", 6)
    elif name == "custom1":
        colors = ["#DD1155", "#7D82B8", "#FCBA04", "#63C7B2", "#45364B"] #Hotpink, lilacky-blue, teal, yellow, dark purple
    elif name == "custom2":
        colors = ["#47C3F4", "#FF9429", "#63C7B2", "#FFCE00", "#EF553C", "#2C819B"] #Lightblue, orange, teal, yellow, red, dark blue, midtone blue, very light blue, 
    colors = (colors*int(math.ceil(N/len(colors))))[:N]
    return colors

def manhattan(values,meta,col_chr=None,col_bp=None,ylabel=None,ax=None,ylim=None,colors=None,annotate=None):
    meta = meta.loc[values.index].copy()
    meta['genomic_position'] = genomic_pos(meta,col_chr,col_bp)
    data = pd.concat([meta,values.to_frame('value')],axis=1)
    Nchr = len(data[col_chr].unique())

    if colors is None:
        colors = chromosome_colors(Nchr)
    elif type(colors) is str:
        colors = chromosome_colors(Nchr, name=colors)
    
    ticks = {}
    for c,g in data.groupby(col_chr):
        bounds = [g.genomic_position.min(),g.genomic_position.max()]
        ticks[c] = bounds[0] + (bounds[1]-bounds[0])/2
    print(ticks)
    
    if ax is None:
        fig, ax = subplots(figsize=(15,8))
        
    which = data.index if ylim is None else (ylim[0] <= data.value) & (data.value <= ylim[1])
    i = 0
    for c,g in data.loc[which].groupby(col_chr):
        if c != "X":
            i += 1
            color = colors[i-1]
        else:
            color = "slategrey"
        sns.scatterplot(g.genomic_position,g.value, color=color, edgecolor=None, marker ='.', rasterized = True)
    
    ax.set_xticks(list(ticks.values()))
    ax.set_xticklabels(list(ticks.keys()))
    ax.tick_params()
    ax.set_xlabel('Genomic Position')
    if not ylabel is None:
        ax.set_ylabel(ylabel)


