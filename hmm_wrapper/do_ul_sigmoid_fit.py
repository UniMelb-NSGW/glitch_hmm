#!/usr/bin/env python

import pulsar_hmm.HMM as HMM
import libstempo
import libstempo.toasim as toasim
import numpy as np
import configparser
import glob
import sys
import argparse
import tempfile
import subprocess
import matplotlib.pyplot as plt
import run_hmm
import scipy.optimize as opt
from lmfit import Model
from copy import deepcopy

def refit_ignore_tn(par, tim, refit_f2=False):
    with tempfile.NamedTemporaryFile(mode='w') as notn_par:
        tn_details = []
        with open(par) as orig_par:
            for line in orig_par.readlines():
                if "TNRed" in line or "MODEL" in line:
                    tn_details.append(line.strip())
                print(line.replace("TNRed", "C TNRed").strip(), file=notn_par)
                print(line.replace("TNRed", "C TNRed").strip())

        notn_par.flush()
        print(notn_par.name)
        lt_psr = libstempo.tempopulsar(notn_par.name, tim)
            
        for p in lt_psr.pars():
            lt_psr[p].fit = False

        lt_psr['F0'].fit = True
        lt_psr['F1'].fit = True
        lt_psr['F2'].fit = refit_f2
        lt_psr.fit()
        lt_psr.savepar(par)

        with open(par, 'a') as new_par:
            for line in tn_details:
                print(line, file=new_par)


def gen_realisation(psr, config, i, glitch_size, tn_params, refit=True, refit_f2=False):
    new_psr = libstempo.tempopulsar(parfile = psr.parfile, timfile = psr.timfile)
    new_psr.name = f'{psr.name}_SYNTH_{i}'
    new_psr['GLF0_1'].val = 10**glitch_size
    new_psr['GLF0_1'].fit = False
    psr_toas = sorted(psr.toas())
    toas = psr_toas

    if 'toas' in config:
        if 'min_toa_gap' in config['toas']:
            toas = [psr_toas[0]]
            for i in range(1, len(psr.toaerrs)):
                if psr_toas[i] - toas[-1] >= float(config['toas']['min_toa_gap']):
                    toas = np.append(toas, psr_toas[i])
        
        if 'mjd_min' in config['toas'] and 'mjd_max' in config['toas']:
            toas = [toa for toa in toas if toa > float(config['toas']['mjd_min']) and toa < float(config['toas']['mjd_max'])]

    #TODO reinstate
    new_psr['GLEP_1'].val = np.random.uniform(toas[2], toas[-2])
    #new_psr['GLEP_1'].val = toas[2] + (toas[-2] - toas[2])/2

    num_preglitch_toas = np.sum([toa < new_psr['GLEP_1'].val for toa in toas])
    toasim.make_ideal(new_psr)
    toasim.add_efac(new_psr)
    print(tn_params)
    if tn_params['red_amp'] is not None:
        toasim.add_rednoise(new_psr, tn_params['red_amp'], tn_params['red_idx'], components=tn_params['red_comp'])

    fake_par = f"{config['ul']['working_dir']}/{new_psr.name}.par"
    fake_tim = f"{config['ul']['working_dir']}/{new_psr.name}.tim"
    glep = new_psr['GLEP_1'].val
    new_psr['GLEP_1'].val = 80000
    new_psr.fit()
    new_psr.fit()
    new_psr.fit()
    new_psr['GLEP_1'].val = glep
    new_psr.savetim(fake_tim)
    new_psr.savepar(fake_par)

    if refit:
        refit_ignore_tn(fake_par, fake_tim, refit_f2)

    return fake_par, fake_tim, num_preglitch_toas

def compute_det_rate(psr, tn_params, config, glitch_size, refit=True, refit_f2=False):
    iters = 0
    working_prefix = config['matlab']['working_prefix']
    num_realisations = int(config['ul']['num_realisations'])
    config['matlab']['working_prefix'] = working_prefix
    glitches = []
    f_fids = []
    fd_fids = []
    fdd_fids = []
    for i in range(num_realisations):
        num_detected = 0
        config['matlab']['working_prefix'] = working_prefix + f"{i}_"
        par, tim, num_preglitch_toas = gen_realisation(psr, config, i, glitch_size, tn_params, refit, refit_f2)
        glitches.append(num_preglitch_toas)
        hmm, _, sigma, _ = run_hmm.setup_hmm(par, tim, config)
        f_fids.append(hmm.f_fiducial)
        fd_fids.append(-hmm.fd_fiducial)
        fdd_fids.append(hmm.fdd_fiducial)
        run_hmm.save_hmm_files(hmm, config)

    extra_matlab_cmd = f"num_realisations={config['ul']['num_realisations']}; glitches = {glitches};\
            f_fids = {f_fids}; fd_fids = {fd_fids}; fdd_fids = {fdd_fids};"
    config['matlab']['working_prefix'] = working_prefix
    run_hmm.do_psr(hmm, sigma, config, extra_matlab_cmd=extra_matlab_cmd, ul=True)
    with open(f"{working_prefix}res.dat", 'r') as f:
        num_detected = int(f.readlines()[0].strip())

    return num_detected/num_realisations

def fit_sigmoid(psr, config, detection_rates):

    glitch_min = float(config['ul']['glitch_min'])
    glitch_max = float(config['ul']['glitch_max'])

    print(detection_rates)

    sizes = [x[0] for x in detection_rates]
    rates = [x[1] for x in detection_rates]
    
    def sigmoid(size, centre, shape):
            return 1/(1 + np.exp(-shape*((size-centre))))

    sigmoid_model = Model(sigmoid)

    result = sigmoid_model.fit(rates, size=sizes, centre=glitch_min + (glitch_max - glitch_min)/2, shape=5)
    result.plot()
    #plt.savefig(f"{config['out']['out_prefix']sigmoid_fit.pdf")
    
    best_fit_fun = lambda size: sigmoid(size, result.params['centre'].value, result.params['shape'].value)

    plt.clf()
    xs = np.linspace(min(sizes), max(sizes), 100)
    sigmoid_curve = [best_fit_fun(x) for x in xs]
    plt.plot(xs, sigmoid_curve)
    plt.plot(sizes, rates, '.', markersize=20)
    plt.xlabel(r"$\log_{10}(\Delta f/\mathrm{Hz})$")
    plt.ylabel(r"Detection rate")
    plt.title(f"{psr.name}")
    plt.tight_layout()
    #plt.savefig(f"sigmoid_fit_custom.png")
    plt.savefig(f"{config['out']['out_prefix']}sigmoid_fit.pdf")

    upper_limit = opt.root_scalar(lambda size: best_fit_fun(size)-0.9, bracket=[glitch_min, glitch_max]).root
    return upper_limit

def compute_ul(psr, tn_params, config):
    iters = 0
    glitch_min = float(config['ul']['glitch_min'])
    glitch_max = float(config['ul']['glitch_max'])
    glitch_size = glitch_min
    working_prefix = config['matlab']['working_prefix']
    num_realisations = int(config['ul']['num_realisations'])
    num_points = int(config['ul']['sigmoid_points'])
    detection_rates = []
    while iters < num_points:
        detection_rates.append((glitch_size, compute_det_rate(psr, tn_params, config, glitch_size)))
        with open(f"{config['ul']['results']}", 'a') as f:
            print(f"{glitch_size} {detection_rates[-1][1]}", file=f)

        if detection_rates[-1][1] > 0.95:
            break
        iters += 1
        glitch_size += (glitch_max - glitch_min)/num_points

    ul = fit_sigmoid(psr, config, detection_rates)
    with open(f"{config['ul']['results']}", 'a') as f:
        print(f"{ul}", file=f)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--par', '-P', type=str, help='The .par file to be used', required=True)
    parser.add_argument('--tim', '-T', type=str, help='The .tim file to be used', required=True)
    parser.add_argument('--ini', '-I', type=str, help='The .ini file containing HMM parameters to be used', required=True)
    parser.add_argument('--split-uls', type=str, help='Directory containing files with individual detection rates to be collated and fit', required=False)

    args = parser.parse_args()
    config = configparser.ConfigParser()
    config.read(args.ini)
    if not args.split_uls:
        hmm, psr, sigma, tn_params = run_hmm.setup_hmm(args.par, args.tim, config)
        compute_ul(psr, tn_params, config)
    else:
        detection_rates = []
        det_rate_files = sorted(glob.glob(f"{args.split_uls}/*dat"))
        for fname in det_rate_files:
            with open(fname) as f:
                line = f.readlines()[0].split()
                detection_rates.append((float(line[0]), float(line[1])))

        with open(config['ul']['results'], 'w') as f:
            for size, rate in detection_rates:
                print(f"{size} {rate}", file=f)

            ul = fit_sigmoid(libstempo.tempopulsar(args.par, args.tim), config, detection_rates)
            print(f"{ul}", file=f)
