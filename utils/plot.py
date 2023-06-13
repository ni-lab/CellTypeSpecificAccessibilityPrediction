import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr


#helper functions
def scatter_lims(vals1, vals2=None, buffer=.05):
    if vals2 is not None:
        vals = np.concatenate((vals1, vals2))
    else:
        vals = vals1
    vmin = np.nanmin(vals)
    vmax = np.nanmax(vals)

    buf = .05 * (vmax - vmin)

    if vmin == 0:
        vmin -= buf / 2
    else:
        vmin -= buf
    vmax += buf

    return vmin, vmax

def scatter_plot(targets, preds, name, ax, alpha=0.5, nonzero_pearson=False, xlim=None, ylim=None):
    sns.set(font_scale=1.2, style='ticks')
    gold = sns.color_palette('husl', 8)[1]
    sns.regplot(targets,preds, color='black',
                order=1,
                scatter_kws={'s': 10,
                             'alpha': alpha},
                line_kws={'color': gold},
                ax=ax)
    
    if not xlim:
        xmin, xmax = scatter_lims(targets)
    else:
        xmin, xmax = xlim
    if not ylim:
        ymin, ymax = scatter_lims(preds)
    else:
        ymin, ymax = ylim
        
    ax.set_title(name, fontsize=24)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel('log2 Experiment')
    ax.set_ylabel('log2 Prediction')
    
    corr, csig = pearsonr(targets, preds)
    corr_str = 'PearsonR: %.3f' % corr
    xlim_eps = (xmax - xmin) * .03
    ylim_eps = (ymax - ymin) * .05

    ax.text(xmin + xlim_eps,
                             ymax - ylim_eps,
                             corr_str,
                             horizontalalignment='left',
                             fontsize=20)
    if nonzero_pearson:
        nonzero_inds = np.where(targets > 0)[0]
        nonzero_corr, nonzero_csig = pearsonr(targets[nonzero_inds], preds[nonzero_inds])
        nonzero_corr_str = 'PearsonR (nonzero): %.3f' % nonzero_corr
        ax.text(xmin + xlim_eps,
                             ymax + 3 * ylim_eps,
                             nonzero_corr_str,
                             horizontalalignment='left',
                             fontsize=20)
    return ax

def joint_plot(targets, preds, name, alpha=0.5, nonzero_pearson=False, xlim=None, ylim=None):
    if not xlim:
        xmin, xmax = scatter_lims(targets)
    else:
        xmin, xmax = xlim
    if not ylim:
        ymin, ymax = scatter_lims(preds)
    else:
        ymin, ymax = ylim
        
    sns.set(font_scale=1.2, style='ticks')
    gold = sns.color_palette('husl', 8)[1]
    h = sns.jointplot(targets,preds, color='black',
                order=1,
                scatter_kws={'s': 10,
                             'alpha': alpha},
                line_kws={'color': gold},
                kind='reg',
                xlim=xlim, ylim=ylim)
    
    
        
    h.fig.suptitle(name, fontsize=16)

    h.ax_joint.set_xlabel('log2 Experiment')
    h.ax_joint.set_ylabel('log2 Prediction')
    
    corr, csig = pearsonr(targets, preds)
    corr_str = 'PearsonR: %.3f' % corr
    xlim_eps = (xmax - xmin) * .03
    ylim_eps = (ymax - ymin) * .05

    h.ax_joint.text(xmin + xlim_eps,
                             ymin + ylim_eps,
                             corr_str,
                             horizontalalignment='left',
                             fontsize=20)
    if nonzero_pearson:
        nonzero_inds = np.where(targets > 0)[0]
        nonzero_corr, nonzero_csig = pearsonr(targets[nonzero_inds], preds[nonzero_inds])
        nonzero_corr_str = 'PearsonR (nonzero): %.3f' % nonzero_corr
        h.ax_joint.text(xmin + xlim_eps,
                             ymin + 3 * ylim_eps,
                             nonzero_corr_str,
                             horizontalalignment='left',
                             fontsize=20)
    return ax