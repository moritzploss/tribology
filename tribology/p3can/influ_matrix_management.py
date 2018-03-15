import glob
import os
import pickle
import shutil

from Constants import SubDir, NpDBs, PrintOpts
from system_functions import make_directory, print_it


def make_infl_mat_db(results_master_folder):
    """Checks if influenc matrix database exists (pickle file). If no,
    creates pickle file with empty dictionary"""
    dir_path = make_directory(results_master_folder,
                              SubDir.infl_mat_db_folder.value)
    db_fname = os.sep.join([dir_path, NpDBs.infl_mat_db.value])
    if not os.path.exists(db_fname):
        with open(db_fname, 'wb') as handle:
            pickle.dump({}, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return db_fname


def cache_influ_mat(ui, influ_mats, res_dir):
    """Checks if inluence matrix is in cache (and database) already. If no,
    influence matrix is added"""
    infl_mat_folder = os.sep.join(
        [os.path.dirname(res_dir), SubDir.infl_mat_db_folder.value])
    infl_data_dict_path = os.sep.join(
        [infl_mat_folder, NpDBs.infl_mat_db.value])

    with open(infl_data_dict_path, 'rb') as handle:
        infl_data_dict = pickle.load(handle)
    for idx, mat in enumerate(influ_mats):
        value = '{}{}'.format(ui['parameter_id'], idx * '-')
        if value in infl_data_dict.values():
            print_it('influence matrix already in cache', PrintOpts.lvl1.value)
        else:
            print_it("caching influence matrix '{}'".format(mat),
                     PrintOpts.lvl1.value)
            infl_data_dict.update({mat: value})
            shutil.copy(os.sep.join([res_dir, SubDir.np_db.value, mat]),
                        os.sep.join([infl_mat_folder, mat]))

    with open(infl_data_dict_path, 'wb') as handle:
        pickle.dump(infl_data_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)


def load_influ_mat(ui, res_dir, number_matrices):
    """Load cached influence matrix and copy it to the current results folder"""
    infl_mat_folder = os.sep.join(
        [os.path.dirname(res_dir), SubDir.infl_mat_db_folder.value])
    infl_data_dict_path = os.sep.join(
        [infl_mat_folder, NpDBs.infl_mat_db.value])

    with open(infl_data_dict_path, 'rb') as handle:
        infl_data_dict = pickle.load(handle)
    key_list = number_matrices * ['empty']
    for i in range(number_matrices):
        for key, value in infl_data_dict.items():
            if '{}{}'.format(ui['parameter_id'], i * '-') == value:
                print_it("loading cached influence matrix '{}'".format(key),
                         PrintOpts.lvl1.value)
                make_directory(res_dir, SubDir.np_db.value)
                shutil.copy(os.sep.join([infl_mat_folder, key]),
                            os.sep.join([res_dir, SubDir.np_db.value, key]))
                key_list[i] = key
    return key_list


def manage_influ_mat_cache(res_dir):
    """Check if the influence matrix cache size exceeds the size limit. If yes,
    remove oldest influence matrix files from cache until size limit is no
    longer exceeded"""
    print_it('clearing influence matrix cache')
    infl_mat_folder = os.sep.join(
        [os.path.dirname(res_dir), SubDir.infl_mat_db_folder.value])
    infl_data_dict_path = os.sep.join(
        [infl_mat_folder, NpDBs.infl_mat_db.value])
    size_cache = 0
    deletable_files = []

    # get size of cache
    for f in sorted(glob.iglob('{}{}*'.format(infl_mat_folder, os.sep)),
                    key=os.path.getctime):
        size_cache += os.path.getsize(f)
        if f != os.sep.join([infl_mat_folder, NpDBs.infl_mat_db.value]):
            deletable_files.extend(f)
    counter = 0

    # remove file from cache if cache size limit is exceeded
    if size_cache / 1024 ** 2 < NpDBs.max_cache_size.value:
        print_it('cache size smaller than cache size limit ({} MB)'.format(
            NpDBs.max_cache_size.value),
                 PrintOpts.lvl1.value)
        print_it('nothing to clear'.format(NpDBs.max_cache_size.value),
                 PrintOpts.lvl1.value)
    else:
        print_it('cache size larger than cache size limit ({} MB)'.format(
            NpDBs.max_cache_size.value),
                 PrintOpts.lvl1.value)
    while size_cache / 1024 ** 2 > NpDBs.max_cache_size.value:
        print_it(
            'deleting influence matrix {}'.format(deletable_files[counter]),
            PrintOpts.lvl1.value)
        os.remove(str(deletable_files[counter]))
        with open(infl_data_dict_path, 'rb') as handle:
            infl_data_dict = pickle.load(handle)
        del infl_data_dict[os.path.basename(str(deletable_files[counter]))]
        with open(infl_data_dict_path, 'wb') as handle:
            pickle.dump(infl_data_dict, handle,
                        protocol=pickle.HIGHEST_PROTOCOL)

        size_cache = 0
        for f in sorted(glob.iglob('{}{}*'.format(infl_mat_folder, os.sep)),
                        key=os.path.getctime):
            size_cache += os.path.getsize(f)
            if f != NpDBs.infl_mat_db.value:
                deletable_files.append(f)
        counter += 1
