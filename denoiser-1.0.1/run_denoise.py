#!/usr/bin/env python
"""
Nuisance regression for fMRI data -- low-memory version
(nilearn 0.9.2 / numpy 1.23 / scipy 1.9 / pandas 1.5 / matplotlib 3.6 / seaborn 0.12).

Saves the cleaned 4D nifti, diagnostic figures, and an html report.

Pipeline: confounds are z-scored and an intercept is added. If --hp_filter and/or
--lp_filter are given, a zero-phase Butterworth bandpass (order --bpf_order, default 2,
as in XCP-D) is applied to the data AND to every confound column, and the regression is
then run on the filtered data (no cosine drift regressors are used). The 'Clean Pre-filter'
carpet comes from a regression without filtering; the T-stat maps and the R-squared map come
from the regression that produced the saved output.

Memory strategy
---------------
* The 4D image is held ONCE, as float32. The cleaned data is written back into that same
  array slab by slab, so there is never a second full-size copy.
* Everything voxelwise (filtering, OLS, R-squared, T-stats, DVARS, carpet sampling) is done
  on chunks of z-slices; the working set is ~6-8x --chunk_mb instead of several full float64
  copies of the whole run.
* OLS is done directly with the pseudo-inverse of the (small) design matrix, so no full-size
  whitened / residual arrays are kept by nilearn's results object.
* DVARS is accumulated from running sums across chunks (no extra nifti files are written to
  disk and re-loaded), and carpet plots use an evenly subsampled set of voxels.
"""
import argparse
import os
import warnings
from os.path import join as pjoin

import matplotlib
matplotlib.use('Agg')  # headless backend, safe on clusters
import matplotlib.pyplot as plt

import jinja2
import nibabel as nb
import numpy as np
import pandas
import scipy.ndimage
import seaborn as sns
from nilearn import plotting
from nilearn.image import load_img
from nilearn.masking import compute_epi_mask
from nilearn.signal import butterworth

NII_EXT = '.nii.gz'
CONSTANT = 'constant'
FD_THR = [.5]                 # framewise displacement threshold(s) for scrubbing preview
SC_RANGE = np.arange(-1, 3)   # volumes flagged around each high-motion volume: -1, 0, +1, +2
TEMPLATE_NAME = 'report_template.html'
CARPET_MAX_VOX = 20000        # max voxels (evenly subsampled) shown in each carpet plot
DVARS_MIN_SD = 1e-3           # voxels with robust SD below this are ignored in DVARS (as nipype)


def strip_nii_ext(name):
    """Remove .nii.gz / .nii from a file name."""
    for ext in ('.nii.gz', '.nii'):
        if name.endswith(ext):
            return name[:-len(ext)]
    return name.split('.')[0]


def expand_scrub_mask(high_motion, buffer=SC_RANGE):
    """Flag each high-motion volume plus its neighbours, staying inside the scan."""
    n = len(high_motion)
    out = high_motion.copy()
    for i in np.where(high_motion)[0]:
        for b in buffer:
            j = i + b
            if 0 <= j < n:
                out[j] = True
    return out


def shade_clusters(ax, mask, bound=.4):
    """Shade contiguous runs of True in `mask` on a time axis."""
    labeled, n_clusters = scipy.ndimage.label(mask)
    for c in range(1, n_clusters + 1):
        idx = np.where(labeled == c)[0]
        ax.axvspan(idx.min() - bound, idx.max() + bound, alpha=.5, color='red')


class DvarsAccumulator:
    """Non-standardized DVARS (std over voxels of the temporal difference), built up chunk by
    chunk from running sums. Matches nipype's compute_dvars(...)[1] with
    remove_zerovariance=True (voxels with robust SD <= 1e-3 are dropped)."""

    def __init__(self, ntrs):
        self.s1 = np.zeros(ntrs - 1)
        self.s2 = np.zeros(ntrs - 1)
        self.n = 0

    def update(self, x):
        """x: voxels x time."""
        if x.shape[0] == 0:
            return
        q1, q3 = np.percentile(x, [25, 75], axis=1, method='lower')
        keep = (q3 - q1) / 1.349 > DVARS_MIN_SD
        if not keep.any():
            return
        d = np.diff(x[keep], axis=1)
        self.s1 += d.sum(axis=0)
        self.s2 += np.square(d).sum(axis=0)
        self.n += d.shape[0]

    def result(self):
        if self.n == 0:
            return np.zeros_like(self.s1)
        mean = self.s1 / self.n
        return np.sqrt(np.maximum(self.s2 / self.n - mean ** 2, 0))


def packed_carpet_img(rows, tr):
    """Pack subsampled voxel rows (n x time) into a tiny 4D image + matching all-ones mask
    so nilearn's plot_carpet can be used without touching the full-size data."""
    arr = np.concatenate(rows, axis=0)
    n, nt = arr.shape
    c_img = nb.Nifti1Image(arr.reshape(n, 1, 1, nt), np.eye(4))
    c_img.header.set_zooms((1, 1, 1, tr))
    c_mask = nb.Nifti1Image(np.ones((n, 1, 1), dtype='uint8'), np.eye(4))
    return c_img, c_mask


def denoise(img_file, tsv_file, out_path, col_names=None, hp_filter=None,
            lp_filter=None, out_figure_path=None, bpf_order=2, chunk_mb=64):
    os.makedirs(out_path, exist_ok=True)

    # ------------------------------------------------------------------ read
    img = load_img(img_file)
    if len(img.shape) != 4:
        raise ValueError('Expected a 4D image, got shape %s' % (img.shape,))
    img_name = os.path.basename(img_file)
    file_base = strip_nii_ext(img_name)
    save_img_file = pjoin(out_path, file_base + '_NR' + NII_EXT)
    Nx, Ny, Nz, Ntrs_img = img.shape

    df_orig = pandas.read_csv(tsv_file, sep='\t', na_values='n/a')
    Ntrs = df_orig.shape[0]
    print('# of TRs: ' + str(Ntrs))
    assert Ntrs == Ntrs_img, \
        'TSV has %d rows but the image has %d volumes' % (Ntrs, Ntrs_img)

    # ----------------------------------------------- select nuisance columns
    if col_names:
        df = df_orig[col_names].copy()
        str_append = '  [SELECTED regressors in CSV]'
    else:
        col_names = df_orig.columns.tolist()
        df = df_orig.copy()
        str_append = '  [ALL regressors in CSV]'

    # fill in missing nuisance values with the mean of that variable
    for col in df.columns:
        n_nan = int(df[col].isnull().sum())
        if n_nan > 0:
            print('Filling in ' + str(n_nan) + ' NaN value for ' + col)
            df[col] = df[col].fillna(df[col].mean())
    print('# of Confound Regressors: ' + str(len(df.columns)) + str_append)

    # z-score confounds (a constant is always in the model, so the fitted/cleaned data
    # are unchanged; this only improves the conditioning of the design matrix)
    df = df.astype(float)
    col_sd = df.std(ddof=0)
    varying_cols = col_sd[col_sd > 0].index
    df[varying_cols] = (df[varying_cols] - df[varying_cols].mean()) / col_sd[varying_cols]
    if len(varying_cols) < df.shape[1]:
        print('Note: %d confound column(s) are constant and were left unscaled: %s'
              % (df.shape[1] - len(varying_cols), [c for c in df.columns if c not in varying_cols]))

    # ------------------------------ design matrix: confounds + intercept
    TR = img.header.get_zooms()[-1]
    Fs = 1. / TR
    do_filter = bool(hp_filter or lp_filter)
    if hp_filter and lp_filter and not (0 < hp_filter < lp_filter):
        raise ValueError('hp_filter must be > 0 and smaller than lp_filter')
    if lp_filter and lp_filter >= Fs / 2:
        raise ValueError('Low pass filter cutoff is too close to the Nyquist frequency (%s)' % (Fs / 2))
    if do_filter:
        print('Butterworth temporal filter (order %d): high-pass %s Hz, low-pass %s Hz'
              % (bpf_order, hp_filter or 'none', lp_filter or 'none'))
    else:
        print('No temporal filtering applied')

    df[CONSTANT] = 1.0

    def report_rank(mat, label):
        rank = np.linalg.matrix_rank(mat)
        print('%s: %d columns, rank %d, condition number %.3g'
              % (label, mat.shape[1], rank, np.linalg.cond(mat)))
        if rank < mat.shape[1]:
            print('WARNING: %d column(s) are linearly dependent on the others '
                  '(e.g. cosine columns or duplicated regressors)' % (mat.shape[1] - rank))
        return rank

    def bandpass(x):
        return butterworth(x, sampling_rate=Fs, low_pass=lp_filter or None,
                           high_pass=hp_filter or None, order=bpf_order, copy=True)

    dm = df.values.astype(float)
    report_rank(dm, 'Design matrix')
    pinv_dm = np.linalg.pinv(dm)

    if do_filter:
        # Apply the SAME Butterworth bandpass to every confound column (Hallquist et al., 2013).
        # The intercept is added back unfiltered (a high-pass would otherwise reduce it to ~0).
        dm_fit = np.hstack((bandpass(dm[:, :-1]), dm[:, -1:]))
        rank_fit = report_rank(dm_fit, 'Filtered design matrix')
        print('Approx. residual degrees of freedom (nominal, before accounting for filtering): %d'
              % (Ntrs - rank_fit))
    else:
        dm_fit = dm
        rank_fit = np.linalg.matrix_rank(dm)
        print('Approx. residual degrees of freedom: %d' % (Ntrs - rank_fit))
    df_resid = Ntrs - rank_fit
    if df_resid <= 0:
        raise ValueError('No residual degrees of freedom left (%d regressors, %d volumes)'
                         % (dm_fit.shape[1], Ntrs))
    pinv_fit = np.linalg.pinv(dm_fit)
    cov_diag = np.einsum('ij,ij->i', pinv_fit, pinv_fit)   # diag((X'X)^-1)

    # ------------------------------------------------ load data ONCE, as float32
    data = img.get_fdata(dtype=np.float32, caching='unchanged')   # (X, Y, Z, T)
    data_mean = data.mean(axis=-1, dtype=np.float64)              # (X, Y, Z)
    varying = np.ptp(data, axis=-1) > 0     # constant voxels make the DVARS stats meaningless
    mean_img = nb.Nifti1Image(data_mean.astype('float32'), img.affine)

    # EPI mask from the mean image (same result as masking the 4D image, without loading it again)
    epi_mask = compute_epi_mask(mean_img).get_fdata() > 0
    mask3d = epi_mask & varying
    print('DVARS/carpet mask: %d voxels (%d constant voxels excluded)'
          % (mask3d.sum(), int(epi_mask.sum() - mask3d.sum())))
    mask_idx = np.flatnonzero(mask3d)
    if mask_idx.size == 0:
        raise ValueError('Empty brain mask')
    step = max(1, int(np.ceil(mask_idx.size / CARPET_MAX_VOX)))
    sel3d = np.zeros(mask3d.size, dtype=bool)
    sel3d[mask_idx[::step]] = True
    sel3d = sel3d.reshape(mask3d.shape)
    del epi_mask, varying, mask_idx

    # which columns get a T-stat map
    plot_cols = [(i, c) for i, c in enumerate(df.columns) if 'drift' not in c and CONSTANT not in c]
    plot_idx = np.array([i for i, _ in plot_cols], dtype=int)
    tmaps = np.zeros((len(plot_cols), Nx, Ny, Nz), dtype='float32')
    rsquare = np.zeros((Nx, Ny, Nz), dtype='float32')

    series = ['in'] + (['pre'] if do_filter else []) + ['out']
    dv = {k: DvarsAccumulator(Ntrs) for k in series}
    rows = {k: [] for k in series}

    # ------------------------------------------------ chunked regression
    dz = int(max(1, min(Nz, (chunk_mb * 2 ** 20) // (Nx * Ny * Ntrs * 8))))
    print('Processing %d slices per chunk (%d chunks)...' % (dz, int(np.ceil(Nz / dz))))
    for z0 in range(0, Nz, dz):
        z1 = min(z0 + dz, Nz)
        dzc = z1 - z0
        y = data[:, :, z0:z1, :].reshape(-1, Ntrs).astype(np.float64)   # voxels x time
        mean_c = data_mean[:, :, z0:z1].reshape(-1)
        m_c = mask3d[:, :, z0:z1].reshape(-1)
        s_c = sel3d[:, :, z0:z1].reshape(-1)

        dv['in'].update(y[m_c])
        rows['in'].append(y[s_c].astype(np.float32))

        if do_filter:
            # regression without temporal filtering (for the 'Clean Pre-filter' carpet / DVARS)
            resid0 = y.T - dm @ (pinv_dm @ y.T)                          # time x voxels
            pre = resid0.T + mean_c[:, None]
            dv['pre'].update(pre[m_c])
            rows['pre'].append(pre[s_c].astype(np.float32))
            del resid0, pre
            y_fit = bandpass(y.T)                                        # time x voxels
        else:
            y_fit = y.T
        del y

        beta = pinv_fit @ y_fit                                          # regressors x voxels
        resid = y_fit - dm_fit @ beta                                    # time x voxels

        # per-voxel sums of squares (on the data that was actually regressed)
        sse = np.sum(np.square(resid), axis=0)
        sst = np.sum(np.square(y_fit - y_fit.mean(axis=0)), axis=0)
        del y_fit
        disp = sse / df_resid

        # T-stat maps for the nuisance regressors
        if len(plot_idx):
            sd = np.sqrt(disp[None, :] * cov_diag[plot_idx][:, None])
            with np.errstate(divide='ignore', invalid='ignore'):
                t = beta[plot_idx] / sd
            t = np.nan_to_num(t, nan=0., posinf=0., neginf=0.)
            tmaps[:, :, :, z0:z1] = t.reshape(-1, Nx, Ny, dzc)
            del sd, t
        del beta

        # R-squared map
        zero_idx = np.logical_and(sst == 0, sse == 0)
        sse[zero_idx] = 1
        sst[zero_idx] = 1          # would be NaN -> becomes rsquared = 0
        with np.errstate(divide='ignore', invalid='ignore'):
            rsq = 1 - np.true_divide(sse, sst)
        rsq[~np.isfinite(rsq)] = 0
        rsquare[:, :, z0:z1] = rsq.reshape(Nx, Ny, dzc)

        # cleaned data (mean added back), written in place over the input slab
        clean = resid.T + mean_c[:, None]
        del resid
        dv['out'].update(clean[m_c])
        rows['out'].append(clean[s_c].astype(np.float32))
        data[:, :, z0:z1, :] = clean.reshape(Nx, Ny, dzc, Ntrs)
        del clean

    if do_filter:
        print('Bandpass applied to data and regressors')

    # ------------------------------------------------------------------ save
    print('Saving output file...')
    hdr = img.header.copy()
    hdr.set_data_dtype(np.float32)      # input may be int16 etc.; never write floats into that
    new_img = nb.Nifti1Image(data, img.affine, header=hdr)
    new_img.to_filename(save_img_file)

    # ------------------------------------------------------------- figures
    fontsize = 12
    fontsize_title = 14
    def_img_size = 8

    if not out_figure_path:
        out_figure_path = pjoin(out_path, file_base + '_NR_figures')
    os.makedirs(out_figure_path, exist_ok=True)
    png_append = '_' + file_base + '.png'
    print('Output directory: ' + out_figure_path)

    # --- correlation matrix of regressors
    cm = df[df.columns[0:-1]].corr()      # last column is the constant
    curr_sz = def_img_size
    if cm.shape[0] > def_img_size:
        curr_sz = curr_sz + ((cm.shape[0] - curr_sz) * .3)
    mtx_scale = curr_sz * 100

    mask = np.zeros_like(cm, dtype=bool)
    mask[np.triu_indices_from(mask)] = True

    fig, ax = plt.subplots(figsize=(curr_sz, curr_sz))
    cmap = sns.diverging_palette(220, 10, as_cmap=True)
    sns.heatmap(cm, mask=mask, cmap=cmap, center=0,
                vmax=cm[cm < 1].max().max(), vmin=cm[cm < 1].min().min(),
                square=True, linewidths=.5, cbar_kws={"shrink": .6}, ax=ax)
    plt.setp(ax.get_xticklabels(), rotation=60, ha='right', fontsize=fontsize)
    plt.setp(ax.get_yticklabels(), rotation=-30, va='bottom', fontsize=fontsize)
    ax.set_title('Nuisance Corr. Matrix', fontsize=fontsize_title)
    fig.tight_layout()
    file_corr_matrix = 'Corr_matrix_regressors' + png_append
    fig.savefig(pjoin(out_figure_path, file_corr_matrix))
    plt.close(fig)

    # --- design matrix
    fig, ax = plt.subplots(figsize=(curr_sz - 4.1, def_img_size))
    x_scale_html = ((curr_sz - 4.1) / def_img_size) * 890
    plotting.plot_design_matrix(df, ax=ax)
    ax.set_title('Nuisance Design Matrix', fontsize=fontsize_title)
    plt.setp(ax.get_xticklabels(), rotation=60, ha='right', fontsize=fontsize)
    plt.setp(ax.get_yticklabels(), fontsize=fontsize)
    ax.set_ylabel('TR (Volume #)', fontsize=fontsize)
    fig.tight_layout()
    file_design_matrix = 'Design_matrix' + png_append
    fig.savefig(pjoin(out_figure_path, file_design_matrix))
    plt.close(fig)

    # --- framewise displacement timeseries
    FD = 'FD'
    bound = .4
    fd_name = next((n for n in ('framewise_displacement', 'FD', 'fd') if n in df_orig.columns), None)
    has_fd = fd_name is not None
    file_fd_plot = None
    y = None
    sc_idx = []
    Nremove = []

    if has_fd:
        fd_series = df_orig[fd_name]
        if fd_series.isnull().any():
            fd_series = fd_series.fillna(fd_series.mean())
        y = fd_series.values
        for thr in FD_THR:
            sc = expand_scrub_mask(y >= thr)
            sc_idx.append(sc)
            Nremove.append(int(sc.sum()))

        Nplots = len(FD_THR)
        sns.set_theme(style='ticks', font_scale=1.5)
        fig, axes = plt.subplots(Nplots, 1, figsize=(def_img_size * 1.5, def_img_size / 2),
                                 squeeze=False)
        sns.despine()
        fd_mean = np.mean(y)
        for curr in range(Nplots):
            a = axes[curr, 0]
            a.plot(y)
            a.plot((-bound, Ntrs + bound), (FD_THR[curr], FD_THR[curr]), '--', color='black')
            a.scatter(np.arange(0, Ntrs), y, s=20)
            if Nremove[curr] > 0:
                shade_clusters(a, sc_idx[curr], bound)
            a.set_ylabel('Framewise Disp. (' + FD + ')')
            a.set_title('%s: %.1f%% of scan (%d volumes) would be scrubbed (FD thr.= %s)'
                        % (FD, 100 * Nremove[curr] / Ntrs, Nremove[curr], FD_THR[curr]))
            a.text(Ntrs + 1, FD_THR[curr] - .01, FD + ' = ' + str(FD_THR[curr]), fontsize=fontsize)
            a.text(Ntrs, fd_mean - .01, 'avg = %.3f' % fd_mean, fontsize=fontsize)
            a.set_xlim((-bound, Ntrs + 8))
        fig.tight_layout()
        axes[-1, 0].set_xlabel('TR (Volume #)')
        file_fd_plot = FD + '_timeseries' + png_append
        fig.savefig(pjoin(out_figure_path, file_fd_plot))
        plt.close(fig)
        print(FD + ' timeseries plot saved')
    else:
        print(FD + ' not found: ' + FD + ' timeseries not plotted')

    # --- carpet and DVARS plots (before & after nuisance regression)
    Ncarpet = 2
    total_sz = 16
    carpet_scale = 840
    y_labels = ['Input (voxels)', "Output 'cleaned'"]
    color = ['red', 'salmon']
    labels = ['input', 'cleaned']
    if do_filter:
        Ncarpet = 3
        total_sz = 20
        carpet_scale = carpet_scale * (9 / 8)
        y_labels = ['Input', 'Clean Pre-filter', 'Clean bandpass']
        color.insert(1, 'firebrick')
        labels.insert(1, 'clean pre-filter')
        labels[-1] = 'clean bandpass'

    dvars = []
    for k in series:
        d = dv[k].result()
        dvars.append(np.hstack((d.mean(), d)))
    del dv

    carpet_imgs = []
    carpet_mask = None
    for k in series:
        c_img, carpet_mask = packed_carpet_img(rows[k], TR)
        carpet_imgs.append(c_img)
    del rows
    print('Carpet plots use %d of %d mask voxels (every %d-th)'
          % (carpet_mask.shape[0], int(mask3d.sum()), step))

    small_sz = 2
    fig = plt.figure(figsize=(def_img_size * 1.5, def_img_size + ((Ncarpet - 2) * 1)))
    row_used = 0
    if has_fd:
        ax0 = plt.subplot2grid((total_sz, 1), (0, 0), rowspan=small_sz)
        ax0.plot(y)
        ax0.scatter(np.arange(0, Ntrs), y, s=10)
        if Nremove[0] > 0:
            shade_clusters(ax0, sc_idx[0], bound)
        ax0.set_ylabel(FD)
        for side in ("top", "right", "bottom"):
            ax0.spines[side].set_visible(False)
        ax0.set_xticks([])
        ax0.set_xlim((-.5, Ntrs - .5))
        ax0.spines["left"].set_position(('outward', 10))
        row_used += small_sz

    ax_d = plt.subplot2grid((total_sz, 1), (row_used, 0), rowspan=small_sz)
    for dvs, c, lab in zip(dvars, color, labels):
        ax_d.plot(dvs, color=c, label=lab)
    ax_d.set_ylabel('DVARS')
    for side in ("top", "right", "bottom"):
        ax_d.spines[side].set_visible(False)
    ax_d.set_xticks([])
    ax_d.set_xlim((-.5, Ntrs - .5))
    ax_d.spines["left"].set_position(('outward', 10))
    ax_d.legend(fontsize=fontsize - 2)
    row_used += small_sz

    st = 0
    carpet_each = int((total_sz - row_used) / Ncarpet)
    for idx, img_curr in enumerate(carpet_imgs):
        ax_curr = plt.subplot2grid((total_sz, 1), (row_used + st, 0), rowspan=carpet_each)
        plotting.plot_carpet(img_curr, mask_img=carpet_mask, figure=fig, axes=ax_curr)
        ax_curr.set_ylabel(y_labels[idx])
        for side in ("bottom", "left"):
            ax_curr.spines[side].set_position(('outward', 10))
        if idx < len(carpet_imgs) - 1:
            ax_curr.spines["bottom"].set_visible(False)
            ax_curr.set_xticklabels('')
            ax_curr.set_xlabel('')
            st += carpet_each

    file_carpet_plot = 'Carpet_plots' + png_append
    fig.savefig(pjoin(out_figure_path, file_carpet_plot))
    plt.close(fig)
    del dvars, carpet_imgs, carpet_mask
    print('Carpet/DVARS plots saved')

    # --- T-stat maps for the nuisance regressors
    mx = None
    if len(plot_cols):
        mx = .8 * max(abs(float(tmaps.min())), abs(float(tmaps.max())))

    t_png = 'Tstat_'
    file_tstat = []
    for k, (idx, col) in enumerate(plot_cols):
        m_img = nb.Nifti1Image(tmaps[k], img.affine)
        title_str = col + ' Tstat'
        display = plotting.plot_stat_map(m_img, mean_img, threshold=3, colorbar=True,
                                         display_mode='z', vmax=mx, title=title_str,
                                         cut_coords=7)
        file_temp = t_png + col + png_append
        display.savefig(pjoin(out_figure_path, file_temp))
        display.close()
        file_tstat.append({'name': col, 'file': file_temp})
        print(title_str + ' map saved')
    del tmaps

    # --- R-squared map
    m_img = nb.Nifti1Image(rsquare, img.affine)
    title_str = 'Nuisance Rsq'
    display = plotting.plot_stat_map(m_img, mean_img, threshold=.2, colorbar=True,
                                     display_mode='z', vmax=.95 * rsquare.max(),
                                     title=title_str, cut_coords=7)
    file_rsq_map = 'Rsquared' + png_append
    display.savefig(pjoin(out_figure_path, file_rsq_map))
    display.close()
    print(title_str + ' map saved')

    # ----------------------------------------------------------- html report
    template_dir = os.path.dirname(os.path.abspath(__file__))
    if not os.path.isfile(pjoin(template_dir, TEMPLATE_NAME)):
        template_dir = os.getcwd()
    if os.path.isfile(pjoin(template_dir, TEMPLATE_NAME)):
        env = jinja2.Environment(loader=jinja2.FileSystemLoader(searchpath=template_dir))
        template = env.get_template(TEMPLATE_NAME)
        output_text = template.render({
            "img_file": img_file,
            "save_img_file": save_img_file,
            "Ntrs": Ntrs,
            "tsv_file": tsv_file,
            "col_names": col_names,
            "hp_filter": hp_filter if hp_filter else False,
            "lp_filter": lp_filter if lp_filter else False,
            "bpf_order": bpf_order,
            "file_design_matrix": file_design_matrix,
            "file_corr_matrix": file_corr_matrix,
            "file_fd_plot": file_fd_plot,
            "file_rsq_map": file_rsq_map,
            "file_tstat": file_tstat,
            "x_scale": x_scale_html,
            "mtx_scale": mtx_scale,
            "file_carpet_plot": file_carpet_plot,
            "carpet_scale": carpet_scale,
        })
        html_file = pjoin(out_figure_path, file_base + '.html')
        with open(html_file, "w") as f:
            f.write(output_text)
        print('')
        print('HTML report: ' + html_file)
    else:
        warnings.warn('%s not found next to the script or in the working directory; '
                      'HTML report skipped (all other outputs were saved).' % TEMPLATE_NAME)

    return new_img


def main():
    parser = argparse.ArgumentParser(
        description='Function for performing nuisance regression. Saves resulting output nifti file, '
                    'information about nuisance regressors and motion (html report), and outputs '
                    'nibabel object containing clean data')
    parser.add_argument('img_file', help='4d nifti img file path')
    parser.add_argument('tsv_file', help='tsv file containing nuisance regressors to be removed')
    parser.add_argument('out_path', help='output directory for saving new data file')
    parser.add_argument('--col_names', nargs="+",
                        help='which columns of TSV file to include as nuisance regressors. '
                             'defaults to ALL columns.')
    parser.add_argument('--hp_filter', type=float,
                        help='frequency cut-off for high pass filter (removing low frequencies). '
                             'Butterworth high-pass edge. Recommend .009 Hz')
    parser.add_argument('--lp_filter', type=float,
                        help='frequency cut-off for low pass filter (removing high frequencies). '
                             'Butterworth low-pass edge. Recommend .08-.1 Hz for non-task data')
    parser.add_argument('--bpf_order', type=int, default=2,
                        help='order of the Butterworth filter (default 2, the XCP-D default). '
                             'Applied forwards and backwards (zero phase).')
    parser.add_argument('--chunk_mb', type=float, default=64,
                        help='size (MB, float64) of the voxel chunk processed at a time; peak extra '
                             'memory is roughly 6-8x this. Lower it to save RAM (default 64).')
    parser.add_argument('--out_figure_path',
                        help='output directory for saving figures. Defaults to <out_path>/<name>_NR_figures')
    args = parser.parse_args()

    denoise(args.img_file, args.tsv_file, args.out_path, args.col_names,
            args.hp_filter, args.lp_filter, args.out_figure_path, args.bpf_order,
            args.chunk_mb)


if __name__ == '__main__':
    main()