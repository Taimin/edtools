import os
import subprocess
import io
import numpy as np
from numpy import fft
import dask
import dask as da
from dask.distributed import Client, LocalCluster, TimeoutError
from scipy.optimize import least_squares, minimize
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from glob import glob
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

from diffractem import version, proc2d, pre_proc_opts, io, tools, proc_peaks
from diffractem.dataset import Dataset
from diffractem.stream_parser import StreamParser, augment_stream, make_substream

def main():
    import argparse

    description = "Program to convert results from DIALS to files required by Crystfel."
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("args",
                        type=str, nargs="*", metavar="FILE",
                        help="List of dials_process.bat files or list of directories. If a list of directories is given "
                        "the program will find all dials_process.bat files in the subdirectories. If no arguments are given "
                        "the current directory is used as a starting point.")

    parser.add_argument("-r", "--read",
                        action="store_true", dest="read",
                        help="Read hdf5 files in h5 directory.")

    parser.add_argument("-comp", "--compute",
                        action="store_true", dest="compute",
                        help="Compute and fit the center beam using Lorentz function.")

    parser.add_argument("-ref_g", "--refine_geometry",
                        action="store_true", dest="refine_geometry",
                        help="Refine experimental data collection geometry.")

    parser.add_argument("-ref_c", "--refine_cell",
                        action="store_true", dest="refine_cell",
                        help="Refine unit cell using powder pattern.")

    parser.add_argument("-plt", "--plot",
                        action="store_true", dest="plot",
                        help="Plot images.")

    parser.add_argument("-uc", "--unit_cell",
                        action="store", type=str, nargs=8, dest="unit_cell",
                        help="Input unit cell, format: lattice cell_type a b c alpha beta gamma.")

    parser.add_argument("-m", "--merge_hkl",
                        action="store_true", dest="merge_hkl",
                        help="Merge hkls from each frame to one file (crystfel format).")

    parser.add_argument("-hit", "--hit_rule",
                        action="store", type=float, nargs=3, dest="hit_rule",
                        help="Rules to judge whether a crystal is hitted by the beam or not, format: #peaks #peaks_lorentz resolution.")

    parser.add_argument("-ind", "--index",
                        action="store_true", dest="index",
                        help="Generate indexing command.")

    parser.add_argument("-int", "--integrate",
                        action="store_true", dest="integrate",
                        help="Rules to judge whether a crystal is hitted by the beam or not, format: #peaks #peaks_lorentz resolution.")

    parser.set_defaults(read=False, compute=False, hit_rule=[15, 10, 0.8], refine_geometry=False, refine_cell=False, 
                        plot=False, merge_hkl=False, index=False, integrate=False)

    options = parser.parse_args()
    fns = options.args
    read = options.read
    compute = options.compute
    refine_geometry = options.refine_geometry
    refine_cell = options.refine_cell
    plot = options.plot
    unit_cell = options.unit_cell
    merge_hkl = options.merge_hkl
    hit_rule = options.hit_rule
    index = options.index
    integrate = options.integrate

    opts = pre_proc_opts.PreProcOpts('preproc.yaml')
    opts.load()
    cluster = LocalCluster(host=None, n_workers=4)
    client = Client(address=cluster)
    #client = Client()

    if read:
        raw_files = io.expand_files('./h5/*data.h5', validate=False)
        print(f'Found {len(raw_files)} raw files. Have fun pre-processing!')
        ds = Dataset.from_files(raw_files, chunking=50, )

    if compute:
        #ds.shots.loc[:, 'selected'] = False
        #ds.shots.loc[1200:, 'selected'] = True
        #ds = ds.get_selection()
        ds.compute_pattern_info(opts='preproc.yaml', client=client, output_file='image_info.h5', calc_center=False)
        selection = f'num_peaks > {hit_rule[0]}'
        # this block is optional... a bit of manual mangling required
        lores_limit = hit_rule[2] # in Angstroms
        pk_D = opts.wavelength/(opts.pixel_size/opts.cam_length)/\
        (((ds.peakXPosRaw - ds.shots.center_x.values.reshape(-1,1))**2
        +(ds.peakYPosRaw - ds.shots.center_y.values.reshape(-1,1))**2)**.5).compute()
        ds.shots['num_lores_peaks'] = (pk_D > lores_limit).sum(axis=1)
        selection += f' and num_lores_peaks > {hit_rule[1]}'
        ds_hit = ds.get_selection(selection, file_suffix='_hit.h5')
        ds_compute = ds_hit # extra step just in case you want to work on a different set
        img_final = proc2d.correct_image(ds_compute.raw_counts, opts,
                                        ds_compute.shots.lor_x.values,
                                        ds_compute.shots.lor_y.values,
                                        ds_compute.peak_data) # keep in mind, that this a lazy computation, so nothing is actually done yet
        ds_compute.add_stack('corrected', img_final, overwrite=True, set_diff_stack=True)

        ds_compute.close_files() # just for safety
        ds_compute.compute_and_save(diff_stack_label='corrected', list_file='hits_agg.lst', exclude_stacks='raw_counts',
                                    client=client, overwrite=True)
        ds_compute.close_files()

    if refine_geometry:
        ds = Dataset.from_files('image_info.h5', chunking=-1)
        pkd = da.compute({k: v for k, v in ds.stacks.items() if 'peak' in k.lower()})[0] # for convenience
        # ellipticity checker
        # radial range to show (in pixels)
        rad_range = (20, 200)
        peakdata = proc_peaks.get_pk_data(pkd['nPeaks'], pkd['peakXPosRaw'], pkd['peakYPosRaw'], 
                                          ds.shots.center_x.values, 
                                          ds.shots.center_y.values, 
                                          opts=opts,
                                          return_vec=False)

        az = np.arctan2((opts.y_scale*peakdata['peakYPosCor']), peakdata['peakXPosCor']).ravel()
        tt = (((opts.y_scale*peakdata['peakYPosCor'])**2 + peakdata['peakXPosCor']**2)**.5).ravel()
        # powder pattern in polar coordinates
        powder_polar = np.histogram2d(tt, az*180/np.pi, 
                                      bins=[np.linspace(*rad_range, 200), np.linspace(-180, 180, 20)])
        corr = powder_polar[0] * np.roll(powder_polar[0], powder_polar[0].shape[1]//4, axis=1)

        if plot:
            fh, ax = plt.subplots(1, 2, figsize=(8,5), sharey=True)
            ax[0].pcolormesh(powder_polar[2][:-1], powder_polar[1][:-1], powder_polar[0])
            ax[0].set_xlabel('Azimuth (deg)')
            ax[0].set_ylabel('Corrected radius')
            ax[1].plot(np.nanmean(powder_polar[0]**2, axis=1), powder_polar[1][:-1], label='Mean squared pattern')
            ax[1].plot(np.nanmean(corr, axis=1), powder_polar[1][:-1], label='Quadrant correlation')
            ax[1].legend()
            ax[1].set_xlabel('Mean squared counts')
            plt.title(f'Median ellipticity variance: {np.nanmedian((np.nanvar(powder_polar[0], axis=1)/np.nanmean(powder_polar[0],axis=1))):.2f} \n'
                     f'Rel. quadrant correlation: {np.mean(corr)/np.mean(powder_polar[0]**2):.3g}')

        angles = np.arange(80, 90, 1)
        ratio = np.arange(1.02, 1.03, 0.0005)
        def cost(p):
            peakdata = proc_peaks.get_pk_data(pkd['nPeaks'], pkd['peakXPosRaw'], pkd['peakYPosRaw'], 
                                          ds.shots.center_x.values, 
                                          ds.shots.center_y.values, 
                                          opts=opts, 
                                          return_vec=False, 
                                              el_rat=p[0], el_ang=p[1])
            az = np.arctan2((opts.y_scale*peakdata['peakYPosCor']), peakdata['peakXPosCor']).ravel()
            tt = (((opts.y_scale*peakdata['peakYPosCor'])**2 + peakdata['peakXPosCor']**2)**.5).ravel()
            powder_polar = np.histogram2d(tt, az*180/np.pi, 
                                          bins=[np.linspace(*rad_range, 200), np.linspace(-180, 180, 20)])
            return np.mean(powder_polar[0]**2)/np.mean(powder_polar[0] * np.roll(powder_polar[0], powder_polar[0].shape[1]//4, axis=1)) - 1

        X, Y = np.meshgrid(ratio, angles)
        with ThreadPoolExecutor() as exc:
            foms = exc.map(cost, zip(X.ravel(), Y.ravel()))
        foms = np.array(list(foms)).reshape(X.shape)
        if plot:
            # show result
            plt.figure()
            plt.pcolormesh(ratio, angles, foms)
            plt.colorbar()
            plt.ylabel('Elliptical axis')
            plt.xlabel('Ellipticity')
            plt.title('Ellipticity correction')
        tools.make_geometry(opts, 'refined.geom')
        ds.update_det_shift('preproc.yaml') # we have image_info.h5 open already... so not use tools
        ds.store_tables(shots=True)
        tools.update_det_shift('hits_agg.lst', 'preproc.yaml')

        if plot:
            #virtual powder
            plt.figure()
            peakdata = proc_peaks.get_pk_data(pkd['nPeaks'], pkd['peakXPosRaw'], pkd['peakYPosRaw'], 
                                              ds.shots.center_x.values, 
                                              ds.shots.center_y.values, pk_I=pkd['peakTotalIntensity'], 
                                              opts=opts, return_vec=True)
            powder, svec = np.histogram(10/peakdata['peakD'].ravel(), bins=np.linspace(0.3,3,1000))
            plt.plot(svec[:-1], powder, 'k');
            plt.xlabel('Scattering vector (1/nm)')
            plt.ylabel('Frequency')
            ax2 = plt.twinx()
            powder, svec = np.histogram(10/peakdata['peakD'].ravel(), bins=np.linspace(0.3,3,1000), 
                                        weights=peakdata['peakTotalIntensity'].ravel(), density=True)
            ax2.plot(svec[:-1], powder, 'r')
            ax2.set_ylabel('Rel. Intensity')
            ax2.yaxis.label.set_color('r')
            plt.title('Scattering vector distribution')

            # Generate peak pair distribution
            d_min = 4
            oversample = 4
            shot_selection = peakdata['nPeaks'] > 20
            out_rad = int((opts.wavelength / d_min) / (opts.pixel_size / opts.cam_length)) + 1
            acfs, r_acfs = proc_peaks.get_acf(peakdata['nPeaks'][shot_selection], 
                                              peakdata['peakXPosCor'][shot_selection,:], 
                                              peakdata['peakYPosCor'][shot_selection,:],
                                              I = None,
                                              output_radius=out_rad, roi_length=512, 
                                              oversample=oversample)

            # scattering vector axis for distances
            s_d = np.arange(r_acfs.shape[1]) * opts.pixel_size/opts.cam_length/(.1*opts.wavelength)/oversample
            fh, axh = plt.subplots(2,1, figsize=(8,5), sharex=True)
            powder, svec = np.histogram(10/peakdata['peakD'].ravel(), bins=np.linspace(0.3,s_d.max(),len(s_d)))
            axh[0].plot(svec[1:]/2 + svec[:-1]/2, powder, color='b')
            axh[0].legend(['Peak resolution distribution'])
            axh[0].grid(True)
            axh[1].plot(s_d,r_acfs.mean(axis=0), color='r')
            axh[1].set_xlabel('Scattering vector (1/nm)')
            axh[1].legend(['Peak pair distance distribution'])
            axh[1].grid(True)
            axh[0].set_title('Virtual powder')

    if refine_cell:
        # refine peaks using cross-correlation, derivative, or distance metric
        # set initial cell
        if unit_cell[0] == 'triclinic':
            C0 = proc_peaks.Cell.triclinic(unit_cell[2], unit_cell[3], unit_cell[4], unit_cell[5], unit_cell[6], unit_cell[7], unit_cell[1])
        if unit_cell[0] == 'monoclinic':
            C0 = proc_peaks.Cell.monoclinic(unit_cell[2], unit_cell[3], unit_cell[4], unit_cell[5], unit_cell[1])
        if unit_cell[0] == 'orthorhombic':
            C0 = proc_peaks.Cell.orthorhombic(unit_cell[2], unit_cell[3], unit_cell[4], unit_cell[1])
        if unit_cell[0] == 'tetragonal':
            C0 = proc_peaks.Cell.tetragonal(unit_cell[2], unit_cell[4], unit_cell[1])
        if unit_cell[0] == 'hexagonal':
            C0 = proc_peaks.Cell.hexagonal(unit_cell[2], unit_cell[4], unit_cell[1])
        if unit_cell[0] == 'rhombohedral':
            C0 = proc_peaks.Cell.rhombohedral(unit_cell[2], unit_cell[5], unit_cell[1])
        if unit_cell[0] == 'cubic':
            C0 = proc_peaks.Cell.cubic(unit_cell[2], unit_cell[1])
        # define minimum d-spacings for the refinement
        dmin_d = 7 # for distances
        dmin_p = 7 # for powder
        # get powder pattern. Tweak histogram bins such that it looks smooth but detailed.
        powder, s_p = np.histogram(10/peakdata['peakD'].ravel(), bins=np.linspace(0.3,10/dmin_p,300), density=True)
        s_p = s_p[1:]/2 + s_p[:-1]/2
        # get average distance distribution
        distance = r_acfs.mean(axis=0)
        C0.init_hkl(dmin_d)
        C_d, par_d = C0.refine_powder(s_d, distance, method='distance', min_prom=0.2, length_bound=2)
        print('Initial cell: ', C0)
        print(f'Distance refinement: {par_d["lsq_result"].x.round(2)}: {par_d["initial_cost"]:.2g} -> {par_d["lsq_result"].cost:.2g}')
        C_d.init_hkl(dmin_p)
        C_p, par_p = C_d.refine_powder(s_p, powder, method='distance', min_prom=0.4, length_bound=2)
        print(f'Powder refinement: {par_p["lsq_result"].x.round(2)}: {par_p["initial_cost"]:.2g} -> {par_p["lsq_result"].cost:.2g}')
        # make a plot for powder-refinement results
        if plot:
            plt.close('all')
            fh, axh = plt.subplots(1,1, figsize=(8,5), sharex=True)
            plt.vlines(10/C0.d(unique=True), 0, powder.max(), color='g', alpha=0.2)
            plt.vlines(10/C_d.d(unique=True), 0, powder.max(), color='r', alpha=0.2)
            plt.vlines(10/C_p.d(unique=True), 0, powder.max(), color='b', alpha=0.2)
            plt.plot(s_p, powder, color='b')
            plt.plot(par_p['peak_position'], par_p['peak_height'], 'bx')
            ax2 = axh.twinx()
            ax2.plot(s_d,distance, color='r')
            ax2.plot(par_d['peak_position'], par_d['peak_height'], 'rx')
            plt.xlim(0, max(10/dmin_d, 10/dmin_p))
            plt.xlabel('Scattering vector (1/nm)')
        C_p.export('refined.cell')

    if index:
        stream_fields = [ 'center_x', 'center_y'] 
        stream_fields = [f'/%/shots/{f}' for f in  stream_fields]
        # generate geometry file for virtual geometry from yaml file parameters.
        tools.make_geometry(opts, 'refined.geom', image_name='corrected', write_mask=False)
        cfcall = tools.call_indexamajig('hits_agg.lst', 'refined.geom', script='im_run_index.sh', 
                           output='master.stream',  cell='refined.cell', im_params=opts.indexing_params, 
                           copy_fields=stream_fields, procs=8)
        print('----------Indexing----------')
        print(cfcall)
        print('----------Indexing----------')

    if integrate:
        # generate sol file
        dsname = 'hits_agg'
        ds = Dataset.from_files(dsname + '.lst', open_stacks=False)
        # Integration
        copy_fields = ['adf1', 'adf2', 'lor_hwhm', 'center_x', 'center_y']
        copy_fields = [f'/%/shots/{cf}' for cf in copy_fields]
        cfcall = tools.call_indexamajig(f'hits_agg.lst', 'refined.geom', 
                                        output=f'hits_agg.stream', 
                                        cell='refined.cell', script='im_run_integrate.sh', 
                                        im_params=opts.integration_params, 
                                        procs=8, exc='indexamajig',
                                        fromfile_input_file = f'indexed.sol',
                                        copy_fields=copy_fields)
        print('----------Integrate----------')
        with open('im_run_integrate.sh', 'r') as f:
            print(f.read())
        print('----------Integrate----------')

    if merge_hkl:
        # get a list of stream files
        stream_list = ['/cygdrive/i/Box Sync/Serial3DED/Only_stage_tilt/large_cover/ZSM-5-exp3-20230917/hits_agg.stream']
        # reject stream files with "all" or "cum" in them, which contain many shots
        # per pattern.
        stream_list = [st for st in stream_list if not (('all' in st) or ('cum' in st))]
        popts = {'no-polarisation': True, 'no-Bscale': False, 'no-scale': False, 
                'force-bandwidth': 2e-5,  'force-radius': False, 'force-lambda': 0.01968,
                    'push-res': 0.8,  'min-measurements': [3, ], 'model': ['unity', 'xsphere'],
                    'symmetry': 'mmm', 'stop-after': list(range(200, 1147, 200)) + [1147],
                    'no-logs': False, 'iterations': 3, 'j': 10}
        cfcall = tools.call_partialator(stream_list, popts, par_runs=4, 
                               split=False, out_dir='merged',
                               slurm=False, cache_streams=False,) 
        print('----------Merging----------')
        print(cfcall)
        print('----------Merging----------')


    cluster.close()