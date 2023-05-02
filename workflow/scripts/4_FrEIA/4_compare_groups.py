#!/usr/bin/env python3
# Author: Norbert Moldován
import os
import sys
import gc
from multiprocessing import Pool
from functools import partial
# import time
import argparse
import math
from itertools import cycle
import matplotlib
matplotlib.use('Agg')  # Matplotlib can't use interactive backend on HPC.
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import scipy.spatial as sp
import scipy.cluster.hierarchy as hc
from scipy.stats import (mannwhitneyu,
                         kruskal,
                         linregress)
#  from statsmodels.stats.multitest import multipletests
from scipy import fft
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from skbio.stats.composition import ilr
import logomaker
from statannot import add_stat_annotation
from FrEIA_tools import (RegroupSamples,
                         DetectControl,
                         GeneratePalette,
                         CastDataTypes)

import pandas as pd



# Parsing the arguments + some nice help text.We
def ParsingArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",
                        dest="ResultsFolder",
                        type=str,
                        required=True,
                        help="The path "
                        "to the folder produced in previous steps.")
    parser.add_argument("-st", "--sampTable",
                        dest="SampleTable",
                        type=str,
                        required=True,
                        help="The path "
                        "to the original sample file.")
    parser.add_argument("-pid", "--pairID",
                        dest="pairID",
                        type=str,
                        required=False,
                        help="The variable "
                        "used to pair the datasets.")
    parser.add_argument("-t", "--threads",
                        dest="threads",
                        type=int,
                        required=True,
                        help="Number of threads.")
    parser.add_argument("-rgr", "--regroup",
                        dest="Regroup",
                        type=str,
                        required=False,
                        default="na",
                        help="The path "
                        "to a sample file which re-groups the samples."
                        "'-cg' and '-sg' are still mandatory.")
    parser.add_argument("-sN", "--subsetN",
                        dest="SubsetN",
                        type=int,
                        required=False,
                        nargs="?",
                        const=0,
                        help="Subset longer results, "
                        "like trinucleotide content to the most abundant N "
                        "motifs. This is useful when plotting large datasets. "
                        "0 results in all elements being plotted. This will "
                        "only constrain the resulting plots and not the "
                        "dataset. Thus the whole dataset is submitted for "
                        "PCA!")
    parser.add_argument("--plot",
                        action="store_true",
                        required=False,
                        help="Output plots. This inkreases runtime!")
    # Option -f may be implemented in the future.
    # parser.add_argument("-f", "--fullChromosomeList",
    #                     dest="fullChrom",
    #                     action="store_true",
    #                     help ="This option will plot every contig "
    #                     "to which reads mapped "
    #                     "when chromosomes are plotted.By default only the 22"
    #                     "autosomes, X and Y are plotted.")

    return parser.parse_args()


# Constant color palettes.
colBase = {'A': 'green',
           'T': 'red',
           'G': 'blue',
           'C': 'orange'}
# List of chromosomes.
ChrList = ["chr1", "chr2", "chr3", "chr4", "chr5",
           "chr6", "chr7", "chr8", "chr9", "chr10",
           "chr11", "chr12", "chr13", "chr14", "chr15",
           "chr16", "chr17", "chr18", "chr19", "chr20",
           "chr21", "chr22", "chrM"]


# Fourier transform.
def PlotFourier(data, minLen, maxLen, amplitude, Groups, Palette):
    FourierDf = pd.DataFrame(columns=["Group", "Power", "Freq"])

    for gr in Groups:
        hFourierDf = pd.DataFrame(columns=["Group", "Power", "Freq"])

        hData = data[(data["Base"] == "G") & (data["WhichGroup"] == gr)]
        # Subsetting the dataframe.
        SubsetDf = hData[(hData["Length"] >= minLen) &
                         (hData["Length"] <= maxLen)]

        # Center signal around zero.
        signal = SubsetDf[amplitude]
        signal = np.array(signal-signal.mean().tolist())

        # Calculate power in the frequency domain
        fft_output = fft.fft(signal)
        power = np.abs(fft_output)
        freq = fft.fftfreq(len(signal))
        mask = freq >= 0
        freq = freq[mask]
        power = power[mask]

        hFourierDf["Power"] = power
        hFourierDf["Freq"] = freq*((maxLen-minLen))
        hFourierDf["Group"] = gr

        FourierDf = FourierDf.append([hFourierDf])
    # print(FourierDf.sort_values(by=["Power"], ascending=False))
    PlotOut = sns.lineplot(data=FourierDf,
                           x="Freq",
                           y="Power",
                           hue="Group",
                           palette=Palette)
    PlotOut.set(xlabel='Frequency (1/nt)', ylabel='Amplitude')
    plt.show()

    return PlotOut


# Calculate the log2 foldchange of control/affected.
def Log2FC(data, ControlGroup, affected, groupBy, value):
    AbDf = data.groupby(["WhichGroup", "Base", groupBy],
                        observed=True)[value].mean().reset_index()
    zipCol = list(zip(AbDf["Base"], AbDf[groupBy]))
    AbDf["ZipCol"] = zipCol
    pivAb = AbDf.pivot(index="ZipCol", columns="WhichGroup", values=value)
    pivAb.dropna(axis=0, inplace=True)
    pivAb.columns = pd.Index(list(pivAb.columns))
    pivAb["s/c"] = (pivAb[affected] / pivAb[ControlGroup])
    pivAb.drop(columns=[affected, ControlGroup], inplace=True)

    hOutDf = pd.DataFrame()
    hOutDf["Base"], hOutDf[groupBy] = zip(*pivAb.index)
    hOutDf["s/c"] = pivAb["s/c"].reset_index(drop=True)

    Log2FCOut = hOutDf.pivot(index="Base", columns=groupBy, values="s/c")
    Log2FCOut = np.log2(Log2FCOut.mask(Log2FCOut <= 0))

    return Log2FCOut


# Calculate motif diversity score per sample.
def MotifDiversityScore(data, groupBy, value, element):
    data = data[data["WhichGroup"] == element]

    # The function for normalized Shannon entropy
    # (a.k.a. motif diversity score).
    motCount = len(data["Base"].unique())
    premdsF = lambda x: -x*np.log2(x) / np.log2(motCount)
    simpsonF = lambda x: x**2
    hMDSDf = data.groupby(["WhichSample", "WhichEnd", "WhichGroup",
                           "Base", groupBy],
                          observed=True)[value].mean().reset_index()

    hMDSDf["MDS"] = hMDSDf[value].apply(premdsF)
    hMDSDf["Simpson"] = hMDSDf[value].apply(simpsonF)

    MDSDf = hMDSDf.groupby(["WhichSample", "WhichEnd", "WhichGroup", groupBy],
                           observed=True)["MDS"].sum().reset_index()
    SimpsonDf = hMDSDf.groupby(["WhichSample", "WhichEnd", "WhichGroup",
                               groupBy],
                               observed=True)["Simpson"].sum().reset_index()
    MDSDf["Gini"] = 1-SimpsonDf["Simpson"]

    return MDSDf


# Hypothesis testing.
def HypoTest(data, group1, Groups, x, y):
    if x != "WhichGroup":
        StatDf = data.groupby([x, "WhichGroup", "WhichSample"],
                              observed=True)[y].mean().reset_index()
    else:
        StatDf = data.groupby([x, "WhichSample"],
                              observed=True)[y].mean().reset_index()
    # For more than 20 samples use the Mann-Whitney U test.
    # For less use the Kruskal-Wallis test.
    GrElement = data["WhichGroup"].iloc[0]
    if len(set(data[data["WhichGroup"] == GrElement]["WhichSample"])) >= 20:
        StatFunc = lambda x, y: mannwhitneyu(x, y, alternative="two-sided")
        testType = "Mann-Whitney_rank_test"
    else:
        StatFunc = lambda x, y: kruskal(x, y)
        testType = "Kruskal-Wallis_H_test"

    # Significance annotation.
    AnnotFunc = lambda x: "ns" if x > 0.05 else("***" if x <= 0.001 else
                                                ("**" if x <= 0.01 else
                                                 ("*")))
    OutDf = pd.DataFrame()

    for xAx in set(data[x]):
        for group2 in Groups:
            if group1 != group2:
                """This was commented out, becuase MDS errored out on the else part.
                if x != "WhichGroup":
                    Gr1=StatDf[(StatDf[x]==xAx)&(StatDf["WhichGroup"]==group1)][y]
                    Gr2=StatDf[(StatDf[x]==xAx)&(StatDf["WhichGroup"]==group2)][y]
                    #print(len(Gr1),"-",Gr1,"\t",len(Gr2),"-",Gr2)
                else:
                    Gr1=StatDf[(StatDf[x]==xAx)][y]
                    Gr2=StatDf[(StatDf[x]==xAx)][y]
                """
                Gr1 = StatDf[(StatDf[x] == xAx) &
                             (StatDf["WhichGroup"] == group1)][y]
                Gr2 = StatDf[(StatDf[x] == xAx) &
                             (StatDf["WhichGroup"] == group2)][y]

                OutDic = {"Test": testType,
                          x: xAx,
                          "Groups": [group1, group2],
                          "P-val": StatFunc(Gr1, Gr2).pvalue,
                          "Stats": StatFunc(Gr1, Gr2).statistic,
                          "Annot": AnnotFunc(StatFunc(Gr1, Gr2).pvalue)}
                OutDf = OutDf.append(OutDic, ignore_index=True)
        #! !!!!!!!!!!Do we need to correct for multiple hypothesis testing?
        # Multiple hypotesis testing correction with the Bonferroni method.
        # if len(Groups)>=20:
        # hOutDf["P-corrected"] = multipletests(OutDf["P-val"],
        #                                        method = 'bonferroni')[1]

    # Cast values to categorical type.
    for c in ["Annot", x, "Test"]:
        OutDf[c] = OutDf[c].astype("category")

    return OutDf


def StatAnnotation(ax, data, x, y, hue, Groups, ControlGroup, Order, HypoDf):
    # Subset to contain only those lines that were plotted.
    SubHypoDf = HypoDf[HypoDf[x].isin(Order)].reset_index(drop=True)
    SubData = data[data[x].isin(Order)].reset_index(drop=True)
    # Plot annotation only if significant.
    SubHypoDf = SubHypoDf[SubHypoDf["Annot"] != "ns"].reset_index(drop=True)
    # Create motif-group pairs. This is expected by add_stat_annotation.
    MotifGroupPairsL = list()
    for pair in zip(SubHypoDf[x], SubHypoDf["Groups"]):
        MotifGroupPairsL.append(tuple(zip(cycle([pair[0]]), pair[1])))
    # Extract the axes element of the plot.
    # ax=plot.axes
    # print(SubHypoDf)
    # print(MotifGroupPairsL)
    # print(len(MotifGroupPairsL))
    PlotOut = add_stat_annotation(ax, data=SubData,
                                  x=x, y=y, hue=hue,
                                  box_pairs=MotifGroupPairsL,
                                  perform_stat_test=False,
                                  pvalues=SubHypoDf["P-val"],
                                  text_annot_custom=SubHypoDf["Annot"],
                                  loc='inside',
                                  linewidth=0.3,
                                  order=Order,
                                  verbose=0)

    return PlotOut


# Plotting scripts.
def CatPlot(data, kind, x, y, hue, order, orientation, row, col, palette,
            height, aspect, legend, legend_out, sharex, sharey,
            xlab, ylab, xrotation, title, doStrip, Groups, ControlGroup):
    sns.set_style("ticks")
    CatFig = sns.FacetGrid(data=data,
                           row=row,
                           col=col,
                           sharex=sharex,
                           sharey=sharey,
                           height=height,
                           aspect=aspect)

    CatFig.map_dataframe(sns.boxplot,
                         x, y,
                         hue=hue,
                         order=order,
                         palette=palette,
                         linewidth=0.5,
                         showfliers=False).add_legend()

    # Create empty boxplots with colored lines.
    ax = CatFig.axes
    # Iterate over plots.
    for r in range(0, len(ax)):
        for i, artist in enumerate(ax[r, 0].artists):
            # Set the linecolor on the artist to the facecolor,
            # and set the facecolor to None.
            col = artist.get_facecolor()
            artist.set_edgecolor(col)
            artist.set_facecolor('None')
            # Each box has 5 associated Line2D objects (to make the whiskers,
            # fliers, etc.).
            # Loop over them here, and use the same colour as above.
            for j in range(i*5, i*5+5):
                line = ax[r, 0].lines[j]
                line.set_color(col)
                line.set_mfc(col)

        HypoDf = pd.DataFrame()
        if (x == "Unit") & ((y == "RelAbSamp") | (y == "RelAbGroup")):
            # Hypothesis testing.
            for u in set(data["Unit"]):
                HypoDfH = HypoTest(data[data["Unit"] == u], ControlGroup,
                                   Groups, "Base", y)
                HypoDfH["Unit"] = u
                HypoDf = HypoDf.append([HypoDfH]).reset_index(drop=True)
        else:
            # Hypothesis testing.
            for e in set(data["WhichEnd"]):
                HypoDfH = HypoTest(data[data["WhichEnd"] == e],
                                   ControlGroup, Groups, x, y)
                HypoDfH["WhichEnd"] = e
                HypoDf = HypoDf.append([HypoDfH]).reset_index(drop=True)

        # Statistical significance plotting.
        if row is not None:
            axRow = ax[r, 0].get_title().split("= ")[-1]
            try:
                StatAnnotation(ax[r, 0], data[data[row] == axRow],
                               x, y, hue, Groups, ControlGroup,
                               order, HypoDf[HypoDf[row] == axRow])
            except:
                pass
        else:
            try:
                StatAnnotation(ax[r, 0], data,
                               x, y, hue, Groups, ControlGroup,
                               order, HypoDf)
            except:
                pass

    # Also fix the legend.
    for legpatch in CatFig.legend.get_patches():
        col = legpatch.get_facecolor()
        legpatch.set_linewidth(1)
        legpatch.set_edgecolor(col)
        legpatch.set_facecolor('None')

    # Overplotting stripplot.
    CatFig.map_dataframe(sns.stripplot,
                         x, y,
                         hue=hue,
                         dodge=True,
                         order=order,
                         palette=palette,
                         size=2,
                         alpha=0.4)

    ax[0, 0].set_title("")

    CatFig.set(xlabel=xlab, ylabel=ylab)
    CatFig.set_xticklabels(rotation=xrotation)
    CatFig.fig.subplots_adjust(top=0.8)
    CatFig.fig.suptitle(title)

    return {"CatFig": CatFig, "HypoDf": HypoDf}


def RelPlot(data, kind, x, y, hue, marker, col, col_wrap, order, palette,
            height, aspect, legend, xlab, ylab, title):
    sns.set_style("ticks")

    RelFig = sns.relplot(data=data,
                         x=x, y=y,
                         hue=hue,
                         marker=marker,
                         kind=kind,
                         col=col,
                         col_wrap=col_wrap,
                         col_order=order,
                         palette=palette,
                         legend=legend,
                         height=height,
                         aspect=aspect)
    RelFig.despine()

    if len(data[x].unique()) <= 15:
        RelFig.set(xticks=data[x].unique())

    RelFig.set(xlabel=xlab, ylabel=ylab)
    RelFig.fig.subplots_adjust(top=0.8)
    RelFig.fig.suptitle(title)
    return RelFig


def RelPlotSubchr(data, kind, isMDS, x, y, hue, col, col_wrap, col_order,
                  palette, height, aspect, legend, xlab, ylab, title):
    def cleanPlot(**kwargs):
        data = kwargs.pop("data")
        plt.margins(x=0)

    def plot_hline(y, **kwargs):
        data = kwargs.pop("data")
        plt.axhline(y=y, c='black', zorder=-1)
    # These will be the width of subplots in inches.
    # Prevent plotting the mitochondria as it is usually missing.
    ChrList.remove("chrM")
    ChrDatapoints = data.groupby(["Chr"],
                                 observed=True
                                 )["Unit"].count().reindex(ChrList)
    AxWidth = list(ChrDatapoints.fillna(0)/100)

    sns.set_style("white")

    if not isMDS:
        RelFig = sns.relplot(data=data,
                             x=x, y=y,
                             hue=hue,
                             kind=kind,
                             col=col,
                             col_wrap=col_wrap,
                             palette=palette,
                             alpha=0.3,
                             s=5,
                             legend=legend,
                             height=height,
                             aspect=aspect,
                             col_order=col_order,
                             facet_kws={'sharex': False,
                                        'sharey': True,
                                        'gridspec_kws': {'width_ratios':
                                                         AxWidth}})
        if (data[y].min() > -1) and (data[y].max() < 1):
            RelFig.set(ylim=(-1, 1))
        RelFig.map_dataframe(plot_hline, y=0)
    else:
        RelFig = sns.relplot(data=data,
                             x=x, y=y,
                             hue=hue,
                             kind=kind,
                             col=col,
                             col_wrap=col_wrap,
                             palette=palette,
                             height=height,
                             aspect=aspect,
                             col_order=col_order,
                             facet_kws={'sharex': False,
                                        'sharey': True,
                                        'gridspec_kws': {'width_ratios':
                                                         AxWidth}})

    RelFig.despine(right=True, top=True, left=True, bottom=True)
    RelFig.set(xlabel=xlab, ylabel=ylab, xticklabels=[])
    RelFig.map_dataframe(cleanPlot)
    RelFig.fig.subplots_adjust(top=0.8, wspace=0.2)
    RelFig.set_titles("{col_name}")
    RelFig.fig.suptitle(title)

    # Add the mitochondria back in the list for other functions.
    ChrList.append("chrM")
    return RelFig


def PCA_analysis(data, components, labels, palette, PlotTitle):
    hDataDf = data.groupby(["WhichGroup", "WhichSample", "Base"],
                           observed=True).RelAbSamp.mean().reset_index()
    hDataDf.sort_values(["WhichSample"], inplace=True)

    PCADf = hDataDf.pivot("WhichSample", "Base", "RelAbSamp")
    PCADf.dropna(axis="columns", inplace=True)
    features = set(PCADf.columns)

    # Separating out the features
    x = PCADf.loc[:, features].values

    xILR = ilr(x)  # Izometric logratio transform because compositional data.

    # Separating out the target
    GroupL = hDataDf.drop_duplicates(subset=["WhichSample"]
                                     )["WhichGroup"].reset_index(drop=True)

    # Standardizing the features
    xILR = StandardScaler().fit_transform(xILR)

    pca = PCA(n_components=len(components[:3]))
    principalComponents = pca.fit_transform(xILR)
    principalDf = pd.DataFrame(data=principalComponents,
                               columns=components[:3])
    PCA_out = pd.concat([principalDf, GroupL], axis=1)
    sns.set_style("ticks")
    plt.axvline(x=0,
                ymin=0,
                ymax=1,
                color="gray",
                linestyle="--")
    plt.axhline(y=0,
                xmin=0,
                xmax=1,
                color="gray",
                linestyle="--")
    PcaFig = sns.scatterplot(data=PCA_out,
                             x="PC0",
                             y="PC1",
                             hue="WhichGroup",
                             alpha=0.5,
                             palette=palette)

    PcaFig.set(title=PlotTitle)

    # Plot screeplot for PCA.
    SPDf = pd.DataFrame({"PC": components[:3],
                         "EigVal": pca.explained_variance_ratio_})
    ScreeFig = SPDf.plot(kind="bar",
                         x="PC",
                         y="EigVal",
                         color="black")

    return {"PCA": PcaFig.get_figure(),
            "Scree": ScreeFig.get_figure()}


# Plot heatmap.
def Heatmap(data, Order, ColisClustered, RowisClustered, ylab, widthK,
            SampColors, colPal, zScore, title):
    # Subsampling dataframe if rows exceed 1000,
    # becaus clustermap can't handle it.
    if data.shape[0] > 1000:
        data.columns = pd.Index(list(data.columns))
        data["rowMeans"] = data.mean(axis=1)
        data.sort_values("rowMeans", inplace=True, ascending=False)
        data.drop(columns=["rowMeans"], inplace=True)
        data = data[:1000]
    # Setting the font size.
    # sns.set(font_scale=1.6)
    # Order motifs
    # data=data.reindex(Order)
    # Make sure that colors and data columns are in the same order.
    # try:
    #    SampColors=list(SampColors.reindex(list(data.columns)))
    # except:
    #    SampColors=None

    width = widthK*len(data.columns)

    # Precalculate distances to avoid error caused by NaN values.
    row_dism = 1 - data.fillna(0).T.corr()
    row_linkage = hc.linkage(sp.distance.squareform(row_dism),
                             method='ward')
    if ColisClustered:
        col_dism = 1 - data.corr()
        col_linkage = hc.linkage(sp.distance.squareform(col_dism),
                                 method='ward')
    else:
        col_linkage = None

    # Plot clustered heatmap.
    HeatFig = sns.clustermap(data=data,
                             mask=data.isin([np.nan, np.inf, -np.inf]),
                             row_cluster=True,
                             row_linkage=row_linkage,
                             col_cluster=ColisClustered,
                             col_linkage=col_linkage,
                             dendrogram_ratio=0.05,
                             center=0,
                             col_colors=SampColors,
                             cmap=colPal,
                             z_score=zScore,
                             figsize=(width, 20),
                             rasterized=True)

    plt.yticks(rotation=0)
    ax = HeatFig.ax_heatmap
    ax.set_xlabel("")
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    ax.set_ylabel("")
    HeatFig.fig.subplots_adjust(top=0.9, bottom=0.5)
    HeatFig.fig.suptitle(title)

    return HeatFig


def SequenceLogo(data, ControlGroup, affected, title):
    for e in set(data["WhichEnd"]):
        L2fcDf = Log2FC(data[data["WhichEnd"] == e],
                        ControlGroup,
                        affected,
                        "Pos",
                        "RelAbGr")
        # Plot logo.
        LogoFig = logomaker.Logo(L2fcDf.transpose())
        LogoFig.style_spines(visible=False)
        LogoFig.style_spines(spines=['left', 'bottom'], visible=True)
        LogoFig.ax.set_xticks(range(0, 10))
        ymin = math.ceil(L2fcDf.min().min())-1
        ymax = math.ceil(L2fcDf.max().max())+1
        yticks = range(ymin, ymax)
        LogoFig.ax.set_ylim(ymin, ymax)
        LogoFig.ax.set_yticks(yticks)
        LogoFig.ax.set_ylabel("Log2(affected/control)", fontsize=12)
        LogoFig.fig.suptitle(t=title, fontsize="medium")

    return LogoFig


def CreateFolders(WhichLvl, args):
    if args.plot:
        if not os.path.exists("".join((args.ResultsFolder,
                                       "4_FrEIA/4_Compare/Plots"))):
            os.makedirs("".join((args.ResultsFolder, "4_FrEIA/4_Compare/Plots")
                                ))
    if not os.path.exists("".join((args.ResultsFolder,
                                   "4_FrEIA/4_Compare/Data"))):
        os.makedirs("".join((args.ResultsFolder, "4_FrEIA/4_Compare/Data")))


def ReadFile(args, WhichGroup, prefix, gr, lvl, filename):
    DataDf = pd.read_parquet("".join((args.ResultsFolder,
                                      "4_FrEIA/3_Abundances/",
                                      gr, "/",
                                      lvl, "/",
                                      filename)))
    DataDf["WhichSample"] = filename.strip(".pq").split("__")[1]

    # Regroup or subset samples.
    if args.Regroup.endswith(".csv"):
        DataDf = RegroupSamples(args.Regroup, DataDf)
    else:
        # If no regrouping than set group as folder name.
        DataDf["WhichGroup"] = gr

    if len(DataDf) > 0:
        DataDf = CastDataTypes(DataDf)
    gc.collect()  # freeing memory for following style_spines
    return DataDf


def ReadData(args, sampTDf, WhichGroup, prefix, lvl):
    DataDf = pd.DataFrame()
    pool = Pool(processes=args.threads)

    for gr in WhichGroup:
        filenameL = []
        if os.path.exists("".join((args.ResultsFolder,
                                   "4_FrEIA/3_Abundances/",
                                   gr, "/",
                                   lvl, "/"))):
            for files in os.walk("".join((args.ResultsFolder,
                                          "4_FrEIA/3_Abundances/",
                                          gr, "/",
                                          lvl, "/"))):
                for filename in files[2]:
                    # print("1***", prefix, filename)
                    if filename.find(prefix) == 0:
                        # print("2***", filename.find(prefix))
                        filenameL.append(filename)
        else:
            sys.exit("File not found error! Input path does not exist. "
                     "Check if folders and files are located under "
                     "[input]/4_FrEIA/3_Abundances/")

        DataDf = DataDf.append(pool.map(partial(ReadFile,
                                                args,
                                                WhichGroup,
                                                prefix,
                                                gr,
                                                lvl),
                                        filenameL),
                               ignore_index=True)
        if len(DataDf) > 0:
            DataDf = CastDataTypes(DataDf)
        gc.collect()  # freeing memory for following steps.
        # print(DataDf.info(memory_usage="deep"))
    pool.close()
    pool.join()

    # Add sample table metadata or regroup table metadata to the dataframe.
    metaDf = sampTDf.drop(["group"], axis=1)
    metaDf = metaDf.rename(columns={"sample_name": "WhichSample"})
    DataDf = DataDf.merge(metaDf, on="WhichSample", how="left")

    # Cast data types for memory and run efficiency.
    DataDf = CastDataTypes(DataDf)

    return DataDf


def SavePlot(plot, path):
    try:
        plot.savefig(path)
    except:
        plt.savefig(path)
    plt.close("all")


# Calculate the fractions per group and per group per read length.
def CalcRelAbGroup(data, prefix):
    # print(data)
    if prefix == "P__":
        data["RelAbGr"] = (data.groupby(["WhichGroup",
                                         "WhichEnd",
                                         "Unit",
                                         "Base",
                                         "Pos"],
                                        observed=True)["CountSamp"
                                                       ].transform("sum")
                             / data.groupby(["WhichGroup",
                                             "WhichEnd",
                                             "Unit",
                                             "Pos"],
                                            observed=True)["CountSamp"
                                                           ].transform("sum"))
    if prefix != "P__":
        data["RelAbGr"] = (data.groupby(["WhichGroup",
                                         "WhichEnd",
                                         "Unit",
                                         "Base"],
                                        observed=True)["CountSamp"
                                                       ].transform("sum")
                            / data.groupby(["WhichGroup",
                                            "WhichEnd",
                                            "Unit"],
                                           observed=True)["CountSamp"
                                                          ].transform("sum"))
        data["RelAbLen"] = (data.groupby(["WhichGroup",
                                          "WhichSample",
                                          "WhichEnd",
                                          "Unit",
                                          "Base",
                                          "Length"],
                                         observed=True)["CountLen"
                                                        ].transform("sum")
                    / data.groupby(["WhichGroup",
                                    "WhichSample",
                                    "WhichEnd",
                                    "Unit",
                                    "Length"],
                                   observed=True)["CountLen"
                                                  ].transform("sum"))
    return data


def MotifOrder(data, orderBy, args):
    # Order motifs according to their fraction.
    Order = data.sort_values(by=[orderBy],
                             ascending=False)["Base"].drop_duplicates()
    Order = Order.reset_index(drop=True)
    # If data is subset for better visualization keel the first
    # args.SubsetN of the series.
    if args.SubsetN != 0:
        Order = Order[:args.SubsetN]
    return Order


def ColorSamples(groupSampleDf, Palette):
    # Pari samples to group colors.
    ColoredSamples = groupSampleDf.drop_duplicates()
    ColoredSamples.set_index(["WhichSample"], inplace=True)
    # P = {k: list(v) for k, v in Palette.items()}
    ColoredSamples = ColoredSamples["WhichGroup"].map(Palette)
    # print(ColoredSamples["WhichGroup"].map(pd.DataFrame.from_dict(Palette)))
    return ColoredSamples


"""
def CreateFigLegends(nbase, WhichPlot, ControlGroup, Groups, hyptest):
    #SignAnnot="The results of the "hyptest" test is only shown for significant differences. (p<0.001: ***; p<0.01: **; p<0.05: *;)"

    "Fig_GP1.svg"       "The base composition of the first " nbase " bases."
    "Fig_GP_",n,".svg"  "The log2 proportion of nucleotides of the first " nbase " bases of the " affected " and " ControlGroup " samples."
    "Fig_GM1.svg"       "The base composition of cohorts based on the first base of the fragments. "SignAnnot
    "Fig_GM2.svg"       "The base composition of cohorts by fragment lentgh."
    "Fig_GM3.svg"       "Principal component analysis of the samples. Variables are the fractions of mononucleotides."
    "Fig_GM4.svg"       "Eigenvalues of the principal components calculated from the fractions of mononucleotides."
    "Fig_GT1.svg"       "The 3-mer composition of cohorts based on the first 3 base of the fragments. "SignAnnot
    "Fig_GT2.svg"       "The 3-mer frequences of cohorts based on the first 3 base of the fragments. The "args.SubsetN" most abundant 3-mers are shown. "
                        "Z-scores were calculated per 3-mer. "
    "Fig_GT3.svg"       "Principal component analysis of the samples. Variables are the fractions of 3-mers."
    "Fig_GT4.svg"       "Eigenvalues of the principal components calculated from the fractions of 3-mers."
    "Fig_GT5.svg"       "Motif diversity scores per length of the 3-mers."
    "Fig_GJM1.svg"      "The joined mononucleotide fractions of cohorts based on the first and last bases of the fragments. "
                        "The "args.SubsetN" most abundant joined ending are shown. "
                        "Z-scores were calculated per joined ending. "
    "Fig_GJM2.svg"      "Principal component analysis of the samples. Variables are the fractions of joined mononucleotide."
    "Fig_GJM3.svg"      "Eigenvalues of the principal components calculated from the fractions of joined mononucleotide."
    "Fig_GJM4.svg"      "Motif diversity scores per length of the joined mononucleotides."
    "Fig_GJT1.svg"      "The joined 3-mer fractions of cohorts based on the first and last 3-mers of the fragments. "
                        "The "args.SubsetN" most abundant joined endings are shown. "
                        "Z-scores were calculated per joined ending. "
    "Fig_GJT2.svg"      "Principal component analysis of the samples. Variables are the fractions of joined 3-mers."
    "Fig_GJT3.svg"      "Eigenvalues of the principal components calculated from the fractions of joined 3-mers."
    "Fig_GJT4.svg"      "Motif diversity scores per length of the joined 3-mers."
    "Fig_CM1.svg"       "The base composition of cohorts based on the first base of the fragments per chromosome. "SignAnnot
    "Fig_CT1.svg"       "The base composition of cohorts based on the first base of the fragments per chromosome. "SignAnnot

pass
"""


def ThreadMDS(args, DataDf, groupBy, value):
    GroupS = set(DataDf["WhichGroup"])
    MDSDf = pd.DataFrame()

    data = DataDf[["WhichSample",
                   "WhichEnd",
                   "WhichGroup",
                   "Base",
                   groupBy,
                   value]]  # Subsetting necesary columns.

    pool = Pool(processes=args.threads)

    try:
        MDSDf = MDSDf.append(pool.map(partial(MotifDiversityScore,
                                              data,
                                              groupBy,
                                              value),
                                      GroupS), ignore_index=True)
        # Creates a list of dataframes.
    except MemoryError:
        pass
    pool.close()
    pool.join()

    return MDSDf


def GenomeLevelAnalysis(DataDf, sampTDf, Groups, ControlGroup,
                        Palette, Order, prefix, lvl, args, MotifAbbr):
    PlotOutPath = "".join((args.ResultsFolder, "4_FrEIA/4_Compare/Plots/"))

    if prefix == "P__":   # START OF POLYNUCLEOTIDE ANALYSIS
        # Save data file.
        DataOutPath = "".join((args.ResultsFolder, "4_FrEIA/4_Compare/Data/"))
        DatOutDf = DataDf.groupby(["Base",
                                   "WhichSample",
                                   "WhichGroup",
                                   "Pos"],
                                  observed=True)["RelAbSamp"
                                                 ].mean().reset_index()
        DatOutDf = DatOutDf.pivot(index=["WhichSample",
                                         "WhichGroup",
                                         "Pos"],
                                  columns="Base",
                                  values="RelAbSamp")
        DatOutDf = DatOutDf[Order]
        DatOutDf.to_csv("".join((DataOutPath, "Dat_G", MotifAbbr, ".csv")))

        if args.plot:
            # Fraction of the first n base.
            Plot = RelPlot(DataDf,
                           "line",
                           "Pos",
                           "RelAbSamp",
                           "WhichGroup",
                           None,
                           "Base",
                           2,
                           Order,
                           Palette,
                           5,
                           1.5,
                           "auto",
                           "Position within read",
                           "Fraction of reads",
                           "")

            SavePlot(Plot, "".join((PlotOutPath, "Fig_G", MotifAbbr, "1.svg")))
            if ControlGroup is not None:
                Groups.remove(ControlGroup)
                n = 2
                for affected in Groups:
                    # Sequence logo of the first n nucleotides
                    # based on the log2 foldchange
                    # between affected/control.
                    Plot = SequenceLogo(DataDf[(DataDf["WhichGroup"] == affected)
                                        | (DataDf["WhichGroup"] == ControlGroup)],
                                        ControlGroup, affected, "")
                    SavePlot(Plot, "".join((PlotOutPath, "Fig_G",
                                            MotifAbbr, str(n), ".svg")))
                    n += 1
                Groups.append(ControlGroup)

    elif prefix == "M__":  # START OF MUNUNUCLEOTIDE ANALYSIS

        # Fourier transform of the signal.
        # PlotFourier(DataDf, 50, 150, "RelAbLen", Groups, Palette)

        # Save data files grouped by samples.
        DataOutPath = "".join((args.ResultsFolder, "4_FrEIA/4_Compare/Data/"))
        DatOutSampDf = DataDf.groupby(["Base",
                                       "WhichSample",
                                       "WhichGroup"],
                                      observed=True)["RelAbSamp"
                                                     ].mean().reset_index()
        DatOutSampDf = DatOutSampDf.pivot(index=["WhichSample",
                                                 "WhichGroup"],
                                          columns="Base",
                                          values="RelAbSamp")
        DatOutSampDf = DatOutSampDf[Order]
        DatOutSampDf.to_csv("".join((DataOutPath, "Dat_G",
                                     MotifAbbr, "sample.csv")))

        # Save data files grouped by length and sample.
        DatOutLenDf = DataDf.groupby(["Base",
                                      "WhichSample",
                                      "WhichGroup",
                                      "Length"],
                                     observed=True)["RelAbLen"
                                                    ].mean().reset_index()
        DatOutLenDf = DatOutLenDf.pivot(index=["WhichSample",
                                               "WhichGroup",
                                               "Length"],
                                        columns="Base",
                                        values="RelAbLen")
        DatOutLenDf = DatOutLenDf[Order]
        DatOutLenDf.to_csv("".join((DataOutPath, "Dat_G",
                                    MotifAbbr, "length.csv")))

        if args.plot:
            # Genome level abundance of the first base.
            BoxDf = DataDf.groupby(["Base",
                                    "WhichSample",
                                    "WhichGroup",
                                    "WhichEnd"],
                                   observed=True)["RelAbSamp"].mean().reset_index()
            Plot = CatPlot(BoxDf,
                           "box",
                           "Base",
                           "RelAbSamp",
                           "WhichGroup",
                           Order,
                           "v",
                           "WhichEnd",
                           None,
                           Palette,
                           3,
                           1.5,
                           True,
                           True,
                           True,
                           True,
                           None,
                           "Fraction of reads",
                           0,
                           "",
                           True,
                           Groups,
                           ControlGroup)
            SavePlot(Plot.get("CatFig"), "".join((PlotOutPath, "Fig_G",
                                                  MotifAbbr, "1.svg")))
            Plot.get("HypoDf").to_csv("".join((DataOutPath, "Dat_G",
                                               MotifAbbr, "hyptest.csv")))
            # Length distribution of the first base.
            Plot = RelPlot(DataDf,
                           "line",
                           "Length",
                           "RelAbLen",
                           "WhichGroup",
                           None,
                           "Base",
                           2,
                           Order,
                           Palette,
                           5,
                           1.5,
                           "auto",
                           "Fragment length",
                           "Fraction of reads",
                           "")
            SavePlot(Plot, "".join((PlotOutPath, "Fig_G", MotifAbbr, "2.svg")))

            # PCA and scree plots.
            pcaL = ["PC0", "PC1", "PC2", "PC3"]
            Plot = PCA_analysis(DataDf, pcaL, None, Palette, "")
            SavePlot(Plot.get("PCA"), "".join((PlotOutPath, "Fig_G",
                                               MotifAbbr, "3.svg")))
            SavePlot(Plot.get("Scree"), "".join((PlotOutPath, "Fig_G",
                                                 MotifAbbr, "4.svg")))

            # Set figure nr.
            n = 5
            if args.pairID:
                # Trajectory plots.
                # Subsample Df for samples with more than one 'pairID' occurrence.
                TimePointSamples = sampTDf.groupby(args.pairID
                                                   ).filter(lambda x: len(x) >= 2)
                traDf = DataDf[DataDf[args.pairID
                                      ].isin(TimePointSamples[args.pairID])
                               ].reset_index()
                Plot = RelPlot(traDf,
                               "line",
                               "timepoint",
                               "RelAbSamp",
                               args.pairID,
                               "o",
                               "Base",
                               2,
                               Order,
                               None,
                               5,
                               1.5,
                               "auto",
                               "",
                               "Fraction of reads",
                               "")
                SavePlot(Plot, "".join((PlotOutPath, "Fig_G", MotifAbbr,
                                        str(n), ".svg")))
                # Set figure nr.
                n += 1

        # Correlations between numerical variables
        # and fragment end base proportions.
        numColL = sampTDf.select_dtypes(np.number).columns
        LRDf = pd.DataFrame()
        for numCol in numColL:
            for gr in Groups:
                for base in Order:
                    maskNA = ~np.isnan(DataDf[(DataDf["WhichGroup"] == gr)
                                              & (DataDf["Base"] == base
                                                 )][numCol]
                                       )
                    LinReg = linregress(DataDf[(DataDf["WhichGroup"] == gr)
                                               & (DataDf["Base"] == base)
                                               ][numCol][maskNA],
                                        DataDf[(DataDf["WhichGroup"] == gr)
                                        & (DataDf["Base"] == base)]["RelAbSamp"
                                                                    ][maskNA])
                    LRDf = LRDf.append({"Group": gr,
                                        "Base": base,
                                        "Metric": numCol,
                                        "Slope": LinReg.slope,
                                        "Intercept": LinReg.intercept,
                                        "R-val": LinReg.rvalue,
                                        "P-val": LinReg.pvalue},
                                       ignore_index=True)
            if args.plot:
                Plot = sns.lmplot(data=DataDf,
                                  x=numCol,
                                  y="RelAbSamp",
                                  hue="WhichGroup",
                                  col="Base",
                                  col_wrap=2,
                                  col_order=Order,
                                  palette=Palette)
                SavePlot(Plot, "".join((PlotOutPath, "Fig_G",
                                        MotifAbbr, str(numCol), str(n), ".png")
                                       ))
        # LRDf.to_csv("".join((DataOutPath, "Dat_G", MotifAbbr, "linreg.csv")))

    else:  # START OF TRI-, JM- & JT-NUNUCLEOTIDE ANALYSIS

        # Save data files grouped by sample.
        DataOutPath = "".join((args.ResultsFolder, "4_FrEIA/4_Compare/Data/"))
        DatOutSampDf = DataDf.groupby(["Base",
                                       "WhichSample",
                                       "WhichGroup"],
                                      observed=True)["RelAbSamp"
                                                     ].mean().reset_index()
        DatOutSampDf = DatOutSampDf.pivot(index=["WhichSample", "WhichGroup"],
                                          columns="Base",
                                          values="RelAbSamp")
        DatOutSampDf = DatOutSampDf[Order]
        DatOutSampDf.to_csv("".join((DataOutPath, "Dat_G",
                                     MotifAbbr, "sample.csv")))

        # Save data files grouped by length and sample.
        DatOutLenDf = DataDf.groupby(["Base",
                                      "WhichSample",
                                      "WhichGroup",
                                      "Length"],
                                     observed=True)["RelAbLen"
                                                    ].mean().reset_index()
        DatOutLenDf = DatOutLenDf.pivot(index=["WhichSample",
                                               "WhichGroup",
                                               "Length"],
                                        columns="Base",
                                        values="RelAbLen")
        DatOutLenDf = DatOutLenDf[Order]
        DatOutLenDf.to_csv("".join((DataOutPath, "Dat_G",
                                    MotifAbbr, "length.csv")))

        # Genome level abundance of the first trinucleotide.
        del DatOutSampDf
        del DatOutLenDf
        gc.collect()  # freeing memory for following steps.

        if args.plot:
            # Set figure nr.
            n = 1

            if prefix == "T__":   # START OF TRINUNUCLEOTIDE ANALYSIS
                BoxDf = DataDf.groupby(["Base",
                                        "WhichSample",
                                        "WhichGroup",
                                        "WhichEnd"],
                                       observed=True)["RelAbSamp"
                                                      ].mean().reset_index()
                Plot = CatPlot(BoxDf,
                               "box",
                               "Base",
                               "RelAbSamp",
                               "WhichGroup",
                               Order,
                               "v",
                               "WhichEnd",
                               None,
                               Palette,
                               5,
                               1.5,
                               True,
                               True,
                               True,
                               True,
                               None,
                               "Fraction of fragments",
                               90,
                               "",
                               True,
                               Groups,
                               ControlGroup)

                SavePlot(Plot.get("CatFig"), "".join((PlotOutPath, "Fig_G",
                                                      MotifAbbr, "1.svg")))
                Plot.get("HypoDf").to_csv("".join((DataOutPath, "Dat_G",
                                                   MotifAbbr, "hyptest.csv")))

                del BoxDf
                gc.collect()  # freeing memory for following steps.

                if args.pairID:
                    # Trajectory plots.
                    # Subsample Df for samples with
                    # more than one 'pairID' occurrence.
                    TimePointSamples = sampTDf.groupby(args.pairID
                                                       ).filter(lambda x:
                                                                len(x) >= 2)
                    # Trajectory by abundance.
                    traDf = DataDf[DataDf[args.pairID
                                          ].isin(TimePointSamples[args.pairID])
                                   ].reset_index()
                    Plot = RelPlot(traDf,
                                   "line",
                                   "timepoint",
                                   "RelAbSamp",
                                   args.pairID,
                                   "o",
                                   "Base",
                                   2,
                                   Order,
                                   None,
                                   5,
                                   1.5,
                                   "auto",
                                   "",
                                   "Fraction of reads",
                                   "")
                    SavePlot(Plot, "".join((PlotOutPath, "Fig_G",
                                            MotifAbbr, "2.svg")))

                gc.collect()  # freeing memory for following steps.
                # Set figure nr.
                n = 3

                # The normalized abundances per sample.
                # Z-scores are calculated per sequence motif (rows).
                PivDf = DataDf.groupby(["Base", "WhichSample"],
                                       observed=True)["RelAbSamp"
                                                      ].mean().reset_index()
                PivDf = PivDf.pivot(index="Base",
                                    columns="WhichSample",
                                    values="RelAbSamp")
                Plot = Heatmap(PivDf,
                               Order,
                               True,
                               True,
                               1,
                               0.5,
                               ColorSamples(DataDf[["WhichGroup", "WhichSample"
                                                    ]],
                                            Palette),
                               "vlag",
                               0,
                               "Fraction of fragments per sample.")
                SavePlot(Plot, "".join((PlotOutPath, "Fig_G",
                                        MotifAbbr, str(n), ".svg")))

                del PivDf
                gc.collect()  # freeing memory for following steps.

                n += 1  # Set figure nr.

            if ControlGroup is not None:
                Groups.remove(ControlGroup)
                for affected in Groups:
                    # The log2fold change by length.
                    Log2FCDf = Log2FC(DataDf[(DataDf["WhichGroup"] == affected)
                                             | (DataDf["WhichGroup"]
                                             == ControlGroup)],
                                      ControlGroup,
                                      affected,
                                      "Length",
                                      "RelAbLen")
                    Plot = Heatmap(Log2FCDf,
                                   Order,
                                   False,
                                   True,
                                   10,
                                   0.1,
                                   None,
                                   "vlag",
                                   None,
                                   "Log2Ratio "+affected+"/"+ControlGroup)
                    SavePlot(Plot, "".join((PlotOutPath, "Fig_G",
                                            MotifAbbr, str(n), ".svg")))

                    n += 1  # Set figure nr.
                Groups.append(ControlGroup)

            # PCA and scree plots of abundances.
            pcaL = ["PC0", "PC1", "PC2", "PC3"]
            Plot = PCA_analysis(DataDf, pcaL, None, Palette, "")
            SavePlot(Plot.get("PCA"), "".join((PlotOutPath, "Fig_G",
                                               MotifAbbr, str(n), ".svg")))
            # Set figure nr.
            n += 1
            SavePlot(Plot.get("Scree"), "".join((PlotOutPath, "Fig_G",
                                                 MotifAbbr, str(n), ".svg")))
            # Set figure nr.
            n += 1

        # Correlations between numerical variables and
        # fragment end motif proportions.
        numColL = sampTDf.select_dtypes(np.number).columns
        LRDf = pd.DataFrame()
        for numCol in numColL:
            for gr in Groups:
                for base in Order:
                    maskNA = ~np.isnan(DataDf[(DataDf["WhichGroup"] == gr)
                                              & (DataDf["Base"] == base)
                                              ][numCol]
                                       )
                    LinReg = linregress(DataDf[(DataDf["WhichGroup"] == gr)
                                               & (DataDf["Base"] == base)
                                               ][numCol][maskNA],
                                        DataDf[(DataDf["WhichGroup"] == gr)
                                               & (DataDf["Base"] == base)
                                               ]["RelAbSamp"][maskNA])
                    LRDf = LRDf.append({"Group": gr,
                                        "Base": base,
                                        "Metric": numCol,
                                        "Slope": LinReg.slope,
                                        "Intercept": LinReg.intercept,
                                        "R-val": LinReg.rvalue,
                                        "P-val": LinReg.pvalue},
                                       ignore_index=True)
            if args.plot:
                Plot = sns.lmplot(data=DataDf,
                                  x=numCol,
                                  y="RelAbSamp",
                                  hue="WhichGroup",
                                  col="Base",
                                  col_wrap=2,
                                  col_order=Order,
                                  palette=Palette)
                SavePlot(Plot, "".join((PlotOutPath, "Fig_G",
                                        MotifAbbr, str(numCol), str(n), ".png")
                                        ))
        #LRDf.to_csv("".join((DataOutPath, "Dat_G", MotifAbbr, "linreg.csv")))

        del LRDf
        gc.collect()  # freeing memory for following steps.

        # Motif diversity score per length.
        # Iterate through Groups separately as MDSDf calculation
        # would take forever otherwise, because of the .apply() function.
        MDSDf = ThreadMDS(args,
                          DataDf,
                          "Length",
                          "RelAbLen")
        if args.plot:
            # Set figure nr.
            n += 1
            Plot = RelPlot(MDSDf,
                           "line",
                           "Length",
                           "MDS",
                           "WhichGroup",
                           None,
                           "WhichEnd",
                           None,
                           None,
                           Palette,
                           5,
                           1.5,
                           "auto",
                           "Fragment length",
                           "Motif diversity score",
                           "")
            SavePlot(Plot, "".join((PlotOutPath, "Fig_G",
                                    MotifAbbr, str(n), ".svg")))
            # Set figure nr.
            n += 1

            # Gini index per length.
            Plot = RelPlot(MDSDf,
                           "line",
                           "Length",
                           "Gini",
                           "WhichGroup",
                           None,
                           "WhichEnd",
                           None,
                           None,
                           Palette,
                           5,
                           1.5,
                           "auto",
                           "Fragment length",
                           "Gini index",
                           "Gini index by fragment lenght")
            SavePlot(Plot, "".join((PlotOutPath, "Fig_G",
                                    MotifAbbr, str(n), ".svg")))
            # Set figure nr.
            n += 1

        # Save data file.
        DataOutPath = "".join((args.ResultsFolder, "4_FrEIA/4_Compare/Data/"))
        MDSDf.rename(columns={"WhichSample": "sample_name",
                              "WhichGroup": "group"},
                     inplace=True)
        MDSDf.to_csv("".join((DataOutPath, "Dat_G",
                              MotifAbbr, "MDS_length.csv")),
                     columns=["sample_name",
                              "group",
                              "Length",
                              "MDS",
                              "Gini"],
                     index=False)

        # Calculating MDS per sample.
        MDSDf_unit = ThreadMDS(args,
                               DataDf,
                               "Unit",
                               "RelAbSamp")

        metaDf = sampTDf.drop(["group"], axis=1)
        metaDf = metaDf.rename(columns={"sample_name": "WhichSample"})
        MDSDf_unit = MDSDf_unit.merge(metaDf, on="WhichSample", how="left")

        if args.plot:
            # Motif diversity of the groups.
            BoxDf = MDSDf_unit.groupby(["WhichSample",
                                        "WhichGroup",
                                        "WhichEnd"],
                                       observed=True)["MDS"].mean().reset_index()
            Plot = CatPlot(BoxDf,
                           "box",
                           "WhichGroup",
                           "MDS",
                           None,
                           None,
                           "v",
                           None,
                           None,
                           Palette,
                           3,
                           1.5,
                           True,
                           True,
                           True,
                           True,
                           None,
                           "MDS",
                           0,
                           "",
                           True,
                           Groups,
                           ControlGroup)
            SavePlot(Plot.get("CatFig"), "".join((PlotOutPath, "Fig_G",
                                                  MotifAbbr, str(n), ".svg")))
            n += 1

            # Gini index of the groups.
            BoxDf = MDSDf_unit.groupby(["WhichSample", "WhichGroup", "WhichEnd"],
                                       observed=True)["Gini"].mean().reset_index()
            Plot = CatPlot(BoxDf,
                           "box",
                           "WhichGroup",
                           "Gini",
                           None,
                           None,
                           "v",
                           None,
                           None,
                           Palette,
                           3,
                           1.5,
                           True,
                           True,
                           True,
                           True,
                           None,
                           "Gini",
                           0,
                           "",
                           True,
                           Groups,
                           ControlGroup)
            SavePlot(Plot.get("CatFig"), "".join((PlotOutPath, "Fig_G",
                                                  MotifAbbr, str(n), ".svg")))
            n += 1

        # Save data file.
        DataOutPath = "".join((args.ResultsFolder, "4_FrEIA/4_Compare/Data/"))
        MDSDat_unit = MDSDf_unit.rename(columns={"WhichSample": "sample_name",
                                                 "WhichGroup": "group"})
        MDSDat_unit.to_csv("".join((DataOutPath, "Dat_G",
                                    MotifAbbr, "MDS_sample.csv")),
                           columns=["sample_name",
                                    "group",
                                    "MDS",
                                    "Gini"],
                           index=False)

        # Correlations between numerical variables and
        # fragment end motif proportions.
        numColL = sampTDf.select_dtypes(np.number).columns
        LRDf = pd.DataFrame()
        for numCol in numColL:
            for gr in Groups:
                maskNA = ~np.isnan(MDSDf_unit[(MDSDf_unit["WhichGroup"] == gr)
                                              ][numCol])
                LinReg = linregress(MDSDf_unit[(MDSDf_unit["WhichGroup"] == gr)
                                               ][numCol][maskNA],
                                    MDSDf_unit[(MDSDf_unit["WhichGroup"] == gr)
                                               ]["MDS"][maskNA])
                LRDf = LRDf.append({"Group": gr,
                                    "Metric": numCol,
                                    "Slope": LinReg.slope,
                                    "Intercept": LinReg.intercept,
                                    "R-val": LinReg.rvalue,
                                    "P-val": LinReg.pvalue}, ignore_index=True)
            if args.plot:
                Plot = sns.lmplot(data=MDSDf_unit,
                                  x=numCol,
                                  y="MDS",
                                  hue="WhichGroup",
                                  palette=Palette)
                SavePlot(Plot, "".join((PlotOutPath, "Fig_G", MotifAbbr,
                                    str(numCol), str(n), ".png")))
        # LRDf.to_csv("".join((DataOutPath, "Dat_G",
        #                     MotifAbbr, "MDS_linreg.csv")))


def ChrLevelAnalysis(DataDf, sampTDf, Groups, ControlGroup, Palette,
                     Order, prefix, lvl, args, MotifAbbr):
    PlotOutPath = "".join((args.ResultsFolder, "4_FrEIA/4_Compare/Plots/"))
    DataOutPath = "".join((args.ResultsFolder, "4_FrEIA/4_Compare/Data/"))

    if prefix == "P__":
        pass
    elif prefix == "M__":
        # Save data file.
        DatOutDf = DataDf.groupby(["Unit",
                                   "Base",
                                   "WhichSample",
                                   "WhichGroup"],
                                  observed=True)["RelAbSamp"
                                                 ].mean().reset_index()
        DatOutDf = DatOutDf.pivot(index=["Unit",
                                         "WhichSample",
                                         "WhichGroup"],
                                  columns="Base",
                                  values="RelAbSamp")
        DatOutDf = DatOutDf[Order]
        DatOutDf.to_csv("".join((DataOutPath, "Dat_C",
                                 MotifAbbr, "sample.csv")))
        if args.plot:
            # Chromosome level abundance of the first base.
            BoxDf = DataDf.groupby(["Unit",
                                    "Base",
                                    "WhichSample",
                                    "WhichGroup",
                                    "WhichEnd"],
                                   observed=True)["RelAbSamp"
                                                  ].mean().reset_index()
            Plot = CatPlot(BoxDf,
                           "box",
                           "Unit",
                           "RelAbSamp",
                           "WhichGroup",
                           ChrList,
                           "v",
                           "Base",
                           None,
                           Palette,
                           3,
                           3,
                           True,
                           True,
                           True,
                           True,
                           None,
                           "Fraction of fragments",
                           45,
                           "",
                           True,
                           Groups,
                           ControlGroup)

            SavePlot(Plot.get("CatFig"), "".join((PlotOutPath, "Fig_CM__1.svg")))
            Plot.get("HypoDf").to_csv("".join((DataOutPath, "Dat_C",
                                               MotifAbbr, "hyptest.csv")))

    else:
        n = 1
        if ControlGroup is not None:
            # Save data file.
            DatOutDf = DataDf.groupby(["Unit",
                                       "Base",
                                       "WhichSample",
                                       "WhichGroup"],
                                      observed=True)["RelAbSamp"
                                                     ].mean().reset_index()
            DatOutDf = DatOutDf.pivot(index=["Unit",
                                             "WhichSample",
                                             "WhichGroup"],
                                      columns="Base",
                                      values="RelAbSamp")
            DatOutDf = DatOutDf[Order]
            DatOutDf.to_csv("".join((DataOutPath, "Dat_C",
                                     MotifAbbr, "sample.csv")))

            if args.plot:
                Groups.remove(ControlGroup)
                for affected in Groups:
                    # The log2fold change by chromosome.
                    Log2FCDf = Log2FC(DataDf[(DataDf["WhichGroup"] == affected) |
                                             (DataDf["WhichGroup"] == ControlGroup)
                                             ],
                                      ControlGroup,
                                      affected,
                                      "Unit",
                                      "RelAbGr")
                    Log2FCDf = Log2FCDf.reindex(ChrList, axis=1)
                    Plot = Heatmap(Log2FCDf,
                                   Order,
                                   False,
                                   True,
                                   10,
                                   0.5,
                                   None,
                                   "vlag",
                                   None,
                                   "Log2Ratio "+affected+"/"+ControlGroup)
                    SavePlot(Plot, "".join((PlotOutPath, "Fig_C",
                                            MotifAbbr, str(n), ".svg")))
                    # Set figure nr.
                    n += 1
                Groups.append(ControlGroup)

        # Motif diversity score per length.
        # Iterate through Groups separately as MDSDf
        # calculation would take forever otherwise,
        # because of the .apply() function.
        MDSDf = ThreadMDS(args,
                          DataDf,
                          "Unit",
                          "RelAbSamp")

        # Save data file.
        MDSDat = MDSDf.rename(columns={"Unit": "Chr",
                                       "WhichSample": "sample_name",
                                       "WhichGroup": "group"})
        MDSDat.to_csv("".join((DataOutPath, "Dat_C",
                               MotifAbbr, "MDS.csv")),
                      columns=["Chr", "sample_name", "group", "MDS"],
                      index=False)
        if args.plot:
            Plot = CatPlot(MDSDf,
                           "box",
                           "Unit",
                           "MDS",
                           "WhichGroup",
                           ChrList,
                           "v",
                           None,
                           None,
                           Palette,
                           3,
                           3,
                           True,
                           True,
                           True,
                           True,
                           None,
                           "Motif diversity score",
                           45,
                           "Chromosome level MDS",
                           True,
                           Groups,
                           ControlGroup)
            SavePlot(Plot.get("CatFig"), "".join((PlotOutPath, "Fig_C",
                                                  MotifAbbr, str(n), ".svg")))
            Plot.get("HypoDf").to_csv("".join((DataOutPath, "Dat_C",
                                               MotifAbbr, "MDS_hyptest.csv")))


def SubChrLevelAnalysis(DataDf, sampTDf, Groups, ControlGroup, Palette,
                        Order, prefix, lvl, args, MotifAbbr):
    PlotOutPath = "".join((args.ResultsFolder, "4_FrEIA/4_Compare/Plots/"))
    DataOutPath = "".join((args.ResultsFolder, "4_FrEIA/4_Compare/Data/"))
    n = 1
    if prefix == "P__":
        pass
    elif prefix == "M__":
        # Save data file.
        DatOutDf = DataDf.groupby(["Unit", "Base",
                                   "WhichSample", "WhichGroup"],
                                  observed=True)["RelAbSamp"
                                                 ].mean().reset_index()
        DatOutDf = DatOutDf.pivot(index=["Unit", "WhichSample", "WhichGroup"],
                                  columns="Base",
                                  values="RelAbSamp")
        DatOutDf = DatOutDf[Order]
        DatOutDf.to_csv("".join((DataOutPath, "Dat_S",
                                 MotifAbbr, "sample.csv")))
        if args.plot:
            if ControlGroup is not None:
                Groups.remove(ControlGroup)
                for affected in Groups:
                    n = 1
                    # The log2fold change.
                    Log2FCDf = Log2FC(DataDf[(DataDf["WhichGroup"] == affected) |
                                             (DataDf["WhichGroup"] == ControlGroup)
                                             ],
                                      ControlGroup,
                                      affected,
                                      "Unit",
                                      "RelAbGr").reset_index()
                    Log2FCDf = pd.melt(Log2FCDf,
                                       id_vars="Base",
                                       value_vars=list(Log2FCDf.columns[1:]),
                                       var_name="Unit",
                                       value_name="Log2FC")
                    Log2FCDf[["Chr", "Unit"]] = Log2FCDf["Unit"
                                                         ].str.split("_",
                                                                     expand=True)
                    Log2FCDf["Unit"] = Log2FCDf["Unit"].astype("int32")
                    Plot = RelPlotSubchr(Log2FCDf,
                                         "scatter",
                                         False,
                                         "Unit",
                                         "Log2FC",
                                         "Base",
                                         "Chr",
                                         None,
                                         ChrList,
                                         colBase,
                                         3,
                                         .3,
                                         "auto",
                                         "",
                                         "Log2Ratio",
                                         "Log2Ratio "+affected+"/"+ControlGroup)
                    SavePlot(Plot, "".join((PlotOutPath, "Fig_S",
                                            MotifAbbr, str(n), ".svg")))
                    n += 1
                Groups.append(ControlGroup)

    else:
        MDSDf = ThreadMDS(args, DataDf, "Unit", "RelAbSamp")

        MDSDf[["Chr", "Unit"]] = MDSDf["Unit"].str.split("_", expand=True)
        MDSDf["Unit"] = MDSDf["Unit"].astype("int32")

        # Save data file.
        MDSDat = MDSDf.rename(columns={"WhichSample": "sample_name",
                                       "WhichGroup": "group"})
        MDSDat.to_csv("".join((DataOutPath, "Dat_S",
                               MotifAbbr, "MDS.csv")),
                      columns=["Chr", "Unit", "sample_name", "group", "MDS"],
                      index=False)

        if args.plot:
            Plot = RelPlotSubchr(MDSDf,
                                 "line",
                                 True,
                                 "Unit",
                                 "MDS",
                                 "WhichGroup",
                                 "Chr",
                                 None,
                                 ChrList,
                                 Palette,
                                 3,
                                 .3,
                                 "auto",
                                 "",
                                 "MDS",
                                 "MDS by region")
            SavePlot(Plot, "".join((PlotOutPath, "Fig_S",
                                    MotifAbbr, str(n), ".svg")))


def SummaryStatistics(DataDf, prefix, lvl, args):
    sumStatFracDf = pd.DataFrame()
    sumStatDiversDf = pd.DataFrame()
    DataOutPath = "".join((args.ResultsFolder, "4_FrEIA/4_Compare/Data/"))

    sumStatFracDf["Mean_Fraction"] = DataDf.groupby(["Unit",
                                                     "WhichGroup",
                                                     "Base"],
                                                    observed=True
                                                    ).mean()["RelAbSamp"]
    sumStatFracDf["Median_Fraction"] = DataDf.groupby(["Unit",
                                                       "WhichGroup",
                                                       "Base"],
                                                      observed=True
                                                      ).median()["RelAbSamp"]
    sumStatFracDf["STD_Fraction"] = DataDf.groupby(["Unit",
                                                    "WhichGroup",
                                                    "Base"],
                                                   observed=True
                                                   ).std()["RelAbSamp"]
    sumStatFracDf["Min_Fraction"] = DataDf.groupby(["Unit",
                                                    "WhichGroup",
                                                    "Base"],
                                                   observed=True
                                                   ).min()["RelAbSamp"]
    sumStatFracDf["Max_Fraction"] = DataDf.groupby(["Unit",
                                                    "WhichGroup",
                                                    "Base"],
                                                   observed=True
                                                   ).max()["RelAbSamp"]

    sumStatFracDf.to_csv("".join((DataOutPath, "Summary_fraction_",
                                  prefix, ".csv")))

    if prefix in ["T__", "JM__", "JT__"]:
        DiversDf = ThreadMDS(args, DataDf, "Unit", "RelAbSamp")
        # DiversDf = MotifDiversityScore(DataDf, "Unit", "RelAbSamp")
        sumStatDiversDf["Mean_MDS"] = DiversDf.groupby("WhichGroup",
                                                       observed=True
                                                       ).mean()["MDS"]
        sumStatDiversDf["Median_MDS"] = DiversDf.groupby("WhichGroup",
                                                         observed=True
                                                         ).median()["MDS"]
        sumStatDiversDf["STD_MDS"] = DiversDf.groupby("WhichGroup",
                                                      observed=True
                                                      ).std()["MDS"]
        sumStatDiversDf["Min_MDS"] = DiversDf.groupby("WhichGroup",
                                                      observed=True
                                                      ).min()["MDS"]
        sumStatDiversDf["Max_MDS"] = DiversDf.groupby("WhichGroup",
                                                      observed=True
                                                      ).max()["MDS"]

        sumStatDiversDf["Mean_Gini"] = DiversDf.groupby("WhichGroup",
                                                        observed=True
                                                        ).mean()["Gini"]
        sumStatDiversDf["Median_Gini"] = DiversDf.groupby("WhichGroup",
                                                          observed=True
                                                          ).median()["Gini"]
        sumStatDiversDf["STD_Gini"] = DiversDf.groupby("WhichGroup",
                                                       observed=True
                                                       ).std()["Gini"]
        sumStatDiversDf["Min_Gini"] = DiversDf.groupby("WhichGroup",
                                                       observed=True
                                                       ).min()["Gini"]
        sumStatDiversDf["Max_Gini"] = DiversDf.groupby("WhichGroup",
                                                       observed=True
                                                       ).max()["Gini"]

        sumStatFracDf.to_csv("".join((DataOutPath, "Summary_diversity_",
                                      prefix, ".csv")))


def Main():
    # Parse arguments.
    args = ParsingArguments()

    WhichLvl = ["Genome_Lvl"]  # Options: "Genome_Lvl", "Chromosome_Lvl", "SubChr_Lvl"
    Prefixes = ["M__", "T__"]  # Options: "M__", "T__", "P__", "JM__", "JT__"

    sampTDf = pd.read_csv(args.SampleTable, delim_whitespace=True)
    InputGroups = DetectControl(sampTDf).get("SampleGroups")
    InputControlGroup = DetectControl(sampTDf).get("ControlGroup")

    Order = pd.Series(dtype="category")
    CreateFolders(WhichLvl, args)

    for prefix in Prefixes:
        for lvl in WhichLvl:
            # print(prefix, lvl)
            DataDf = ReadData(args, sampTDf, InputGroups, prefix, lvl)
            DataDf = CalcRelAbGroup(DataDf, prefix)
            if args.Regroup.endswith(".csv"):
                regrTDf = pd.read_csv(args.Regroup, delim_whitespace=True)
                ControlGroup = DetectControl(regrTDf).get("ControlGroup")
                Groups = DetectControl(regrTDf).get("SampleGroups")
            else:
                ControlGroup = InputControlGroup
                Groups = InputGroups
            Palette = GeneratePalette(Groups, ControlGroup)

            # print(DataDf.info(memory_usage="deep"), "\n")
            if lvl == "Genome_Lvl":
                # SummaryStatistics(DataDf, prefix, lvl, args)
                Order = MotifOrder(DataDf,
                                   "RelAbGr",
                                   args)
                GenomeLevelAnalysis(DataDf,
                                    sampTDf,
                                    Groups,
                                    ControlGroup,
                                    Palette,
                                    Order,
                                    prefix,
                                    lvl,
                                    args,
                                    prefix)
            elif lvl == "Chromosome_Lvl":
                # SummaryStatistics(DataDf,
                #                    prefix,
                #                    lvl,
                #                    args)
                Order = MotifOrder(DataDf,
                                   "RelAbGr",
                                   args)
                ChrLevelAnalysis(DataDf,
                                 sampTDf,
                                 Groups,
                                 ControlGroup,
                                 Palette,
                                 Order,
                                 prefix,
                                 lvl,
                                 args,
                                 prefix)
            elif (len(DataDf) != 0) & (lvl == "SubChr_Lvl"):
                Order = MotifOrder(DataDf,
                                   "RelAbGr",
                                   args)
                SubChrLevelAnalysis(DataDf,
                                    sampTDf,
                                    Groups,
                                    ControlGroup,
                                    Palette,
                                    Order,
                                    prefix,
                                    lvl,
                                    args,
                                    prefix)


if __name__ == "__main__":
    Main()
