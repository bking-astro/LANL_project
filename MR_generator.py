import numpy as np
import pandas as pd
import os
import time

import TOVsolver
import EOS_extension

def make_folders(datapath, ext_type, nsamp_EOS):
    EOSdir_name = ext_type + str(nsamp_EOS) + 'EOS'
    MRLdir_name = ext_type + str(nsamp_EOS) + 'MRL'

    if not(EOSdir_name in os.listdir(datapath)) and not(MRLdir_name in os.listdir(datapath)):
        os.makedirs(datapath+EOSdir_name)
        os.makedirs(datapath+MRLdir_name)

    return EOSdir_name, MRLdir_name


def generate_tables(datapath, number_tables, nsamp_EOS, ext_type, MRL_size=100, maxm_thresh=1.8, meanm_thresh=0.8):
    # make directories for the EOS and MRL if they don't exist already
    EOSdir_name, MRLdir_name, filenumstart = make_folders(datapath, ext_type, nsamp_EOS)

    EOS_start = pd.read_table(datapath + '/EOSCEFTVE1.dat', header=None).to_numpy()

    start_time = time.time()

    i = 0
    done_25 = False
    done_50 = False
    done_75 = False
    num_failed = 0
    while i < number_tables:

        if ext_type == 'poly':
            EOS, ns, gammas, Ks = EOS_extension.extend(EOS_start, nsamp=nsamp_EOS,
                                                       ext_type=ext_type, max_gamma=9)
            param_string = "ns =" + str(ns) + ' gammas =' + str(gammas) + ' Ks =' + str(Ks)
        else:
            EOS, ns, cs = EOS_extension.extend(EOS_start, nsamp=nsamp_EOS, ext_type=ext_type)
            param_string = "ns =" + str(ns) + ' cs =' + str(cs)

        # extend EOS and store parameters

        MRL_table = TOVsolver.solve(EOS, MRL_size)  # solve tov
        raw_mass = MRL_table[:, 0]
        raw_radius = MRL_table[:, 1]
        raw_Lambda = MRL_table[:, 2]

        # create boolean arrays to test if points are good
        m2big = raw_mass < 4
        r2small = raw_radius > 7
        # the bool array of points we will keep
        keep = np.logical_and(m2big, r2small)
        # define new arrays we will keep
        radius = raw_radius[keep]
        mass = raw_mass[keep]
        Lambda = raw_Lambda[keep]

        # check if maximum mass is realistic
        maxm = np.max(mass)
        meanm = np.mean(mass)

        if maxm < maxm_thresh or meanm < meanm_thresh:
            #             print("maximum mass is too low")
            num_failed += 1
            i -= 1
        else:
            leng = len(radius)  # get number of physical points
            MRL = np.zeros((leng, 3))  # initialize MRL table
            MRL[:, 0], MRL[:, 1], MRL[:, 2] = mass, radius, Lambda  # put into table

            EOSname = datapath + EOSdir_name + '/' + str(filenumstart + i) + '.dat'  # make names for file
            MRLname = datapath + MRLdir_name + '/' + str(filenumstart + i) + '.dat'

            #             print("saving file...")

            np.savetxt(EOSname, EOS, header=param_string)  # save files
            np.savetxt(MRLname, MRL)

        if i * 100 / number_tables > 25 and done_25 == False:
            print("The run is 25% complete")
            done_25 = True
        elif i * 100 / number_tables > 50 and done_50 == False:
            print("The run is 50% complete")
            done_50 = True
        elif i * 100 / number_tables > 75 and done_75 == False:
            print("The run is 75% complete")
            done_75 = True

        i += 1
    #     print(str(num_failed)+' iterations failed for '+str(number_tables)+' successful tables.')
    print("Generating "+str(number_tables)+" tables took "+str(time.time()-start_time)+" seconds to complete.")


def generate(datapath, nsamp_EOS, ext_type, MRL_size=100, maxm_thresh=1.8, meanm_thresh=0.8, newsavename=None):
    # make directories for the EOS and MRL if they don't exist already
    EOSdir_name, MRLdir_name, filenumstart = make_folders(datapath, ext_type, nsamp_EOS)

    EOS_start = pd.read_table(datapath + '/EOSCEFTVE1.dat', header=None).to_numpy()

    # when creation started
    start_time = time.time()

    # initialize mean and max M
    meanm = 0
    maxm = 0

    while meanm < meanm_thresh or maxm < maxm_thresh:

        # make correct
        if ext_type == 'poly':
            EOS, ns, gammas, Ks = EOS_extension.extend(EOS_start, nsamp=nsamp_EOS,
                                                       ext_type=ext_type, max_gamma=9)
            param_string = "ns =" + str(ns) + ' gammas =' + str(gammas) + ' Ks =' + str(Ks)
        else:
            EOS, ns, cs = EOS_extension.extend(EOS_start, nsamp=nsamp_EOS, ext_type=ext_type)
            param_string = "ns =" + str(ns) + ' cs =' + str(cs)

        # extend EOS and store parameters

        MRL_table = TOVsolver.solve(EOS, MRL_size)  # solve tov
        raw_mass = MRL_table[:, 0]
        raw_radius = MRL_table[:, 1]
        raw_Lambda = MRL_table[:, 2]

        # create boolean arrays to test if points are good
        m2big = raw_mass < 4
        r2small = raw_radius > 7
        # the bool array of points we will keep
        keep = np.logical_and(m2big, r2small)
        # define new arrays we will keep
        radius = raw_radius[keep]
        mass = raw_mass[keep]
        Lambda = raw_Lambda[keep]

        # check if maximum mass and mean mass is realistic
        maxm = np.max(mass)
        meanm = np.mean(mass)

    leng = len(radius)  # get number of physical points
    MRL = np.zeros((leng, 3))  # initialize MRL table
    MRL[:, 0], MRL[:, 1], MRL[:, 2] = mass, radius, Lambda  # put into table

    if newsavename != None:
        EOSdir_name, MRLdir_name, filenumstart = make_folders(datapath, ext_type, nsamp_EOS, newsavename)

    EOSname = datapath + EOSdir_name + '/' + str(filenumstart) + '.dat'  # make names for file
    MRLname = datapath + MRLdir_name + '/' + str(filenumstart) + '.dat'

    np.savetxt(EOSname, EOS, header=param_string)  # save files
    np.savetxt(MRLname, MRL)

    print("Generation took " + str(time.time() - start_time) + " seconds to complete.")