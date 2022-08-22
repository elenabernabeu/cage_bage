#!/usr/bin/env python
# -*- coding: utf-8 -*-
# By Elena

import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.patches as mpatches
import matplotlib.colors as mc
import seaborn as sns
import random
import math
import colorsys

######################################################
# Separate plates into folds

# Font config
font = "Galvji"
path = '/Cluster_Filespace/Marioni_Group/Elena/other/fonts/Galvji/Galvji-01.ttf'
fe = matplotlib.font_manager.FontEntry(
    fname= path,
    name=font)
matplotlib.font_manager.fontManager.ttflist.insert(0, fe)
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.family'] = font

def flatten(t):
    return [item for sublist in t for item in sublist]

def lighten_color(color, amount=0.5):  
    # --------------------- SOURCE: @IanHincks ---------------------
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])

def set_colors(N, name = "default"):
    colors = []
    if name == "default":
        colors = sns.color_palette("Set2", 6)
    elif name == "custom1":
        colors = ["#DD1155", "#7D82B8", "#FCBA04", "#63C7B2", "#45364B"] #Hotpink, lilacky-blue, teal, yellow, dark purple
    elif name == "custom2":
        colors = ["#47C3F4", "#FF9429", "#63C7B2", "#FFCE00", "#EF553C", "#2C819B", "#07B1D8", "#C5E3EA", "#326B64", "#FFD15F", "#EA4C15"] #Lightblue, orange, teal, yellow, red, dark blue, midtone blue, very light blue, 
    colors = (colors*int(math.ceil(N/len(colors))))[:N]
    return colors

def get_folds(df, plate_dictionary, platesperfold_dictionary, output_dir, name, waves, plate_col = "Batch", index_col = "Sample_Sentrix_ID", wave_col = "Set"):
    df["Fold"] = 0
    t = flatten([platesperfold_dictionary[i] for i in waves]) # Flatten all folds into a single list

    z = 1 # Fold counter
    for wave in waves:
        print("############ Wave %s." % wave)
        df_w = df[df[wave_col] == "W%s" % wave]
        print("A total of %s indivs." % len(df_w.index))
        plates_w = sorted(list(set(df_w[plate_col])))
        print("A total of %s plates." % len(plates_w))
        
        # Export number of samples per plate
        n_plates = df_w[plate_col].value_counts()
        n_plates.to_csv(output_dir + "plate_counts_wave%s_%s.tsv" % (wave, name), sep = "\t", header = False)

        j = 0 # Plate counter
        base = ["W%s_%s" % (wave, i) for i in plate_dictionary[wave]]
        for i in platesperfold_dictionary[wave]:
            plates4fold = base[j:(j+i)]
            df.loc[df[plate_col].isin(plates4fold), "Fold"] = z # Add fold to df
            print("Plates in fold %s: %s" % (z, ", ".join(plates4fold)))
            j += i
            z += 1
    
    # Export table with total number of individuals and plates per fold
    df_ppf = pd.DataFrame(index = range(1, len(t)+1), columns = ["nPlates", "nSamples"])
    for f in range(1, len(t)+1):
        df_ppf.loc[f, "nPlates"] = t[f-1]
        df_ppf["nSamples"] = df["Fold"].value_counts()
    df_ppf.to_csv(output_dir + "fold_counts_%s.tsv" % name, sep = "\t", index_label = "Fold")
    return df["Fold"]


def explore_folds(df, name, output_dir, fold_col = "Fold", age_col = "age", sex_col = "sex", wave_col = "Set", plate_col = "Batch", color = "custom2", external = "F"):
    # Aesthetic
    palette = set_colors(6, color)
    sns.set_style("white", {'legend.frameon':True})
    sns.set_palette(sns.color_palette(palette))
    matplotlib.rcParams['font.family'] = font

    # Keep essentials
    df = df[[fold_col, age_col, sex_col]]

    # Make fold column integers
    df[fold_col] = pd.to_numeric(df[fold_col], downcast='integer')
    
    # Initialize plot
    fig, axes = plt.subplots(2, 1, figsize = (7, 5.5))
    plt.suptitle("Elastic Net CV Folds")

    # Boxplot
    col = palette[0]
    sns.boxplot(ax = axes[0], data = df, x = fold_col, y = age_col, linewidth = 1, saturation = 1, color = col)  
    axes[0].set_ylabel("Age")
    axes[0].set_xlabel("")
    axes[0].set_xticklabels([])
    for i,artist in enumerate(axes[0].artists):
        # Set the linecolor on the artist to the facecolor, and set the facecolor to None
        col = lighten_color(artist.get_facecolor(), 1.5)
        artist.set_edgecolor(col)    

        # Each box has 6 associated Line2D objects (to make the whiskers, fliers, etc.)
        # Loop over them here, and use the same colour as above
        for j in range(i*6,i*6+6):
            line = axes[0].lines[j]
            line.set_color(col)
            line.set_mfc(col)
            line.set_mec(col)

    # Barplot number of samples and sex
    count = df.groupby(fold_col)[sex_col].value_counts()
    count = count.rename('Count').reset_index()  
    sns.barplot(ax = axes[1], x = fold_col, y = "Count", data = count, hue = sex_col, saturation = 1, linewidth = 0)
    axes[1].set_xlabel("Fold")
    axes[1].get_legend().remove()
    axes[1].set_ylim(bottom = 0, top = 800)

    # Final bits and legend
    elements = [mpatches.Patch(edgecolor = palette[0], facecolor = palette[0], label = "F"), mpatches.Patch(edgecolor = palette[1], facecolor = palette[1], label = "M")]
    fig.legend(handles=elements, loc = 3, ncol = 2, fancybox = False)
    plt.tight_layout(rect=[0, 0.03, 1, 0.94])
    fig.subplots_adjust(bottom=0.13)
    
    # Export
    if external == "F":
        plt.savefig(output_dir + "fold_age.sex_%s.pdf" % name)
    else:
        plt.savefig(output_dir + "fold_age.sex_%s_withexternal.pdf" % name)
    plt.close(fig)




    

