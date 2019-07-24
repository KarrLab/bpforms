""" Parse PDB database and analyze non-canonical amino acid

:Author: Mike Zheng <xzheng20@colby.edu>
:Date: 2019-07-24
:Copyright: 2019, Karr Lab
:License: MIT
"""

from ete3 import NCBITaxa
from ftplib import FTP
import gzip
from io import BytesIO
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import random
from ruamel import yaml

PDB_FTP_HOST = 'ftp.wwpdb.org'
PDB_ROOT = '/pub/pdb/data/structures/divided/pdb/'
PROTEIN_YML_PATH = '../bpforms/alphabet/protein.yml'

# get org_ids for ecoli and descendants
ncbi = NCBITaxa()
ecoli_ids = ncbi.get_descendant_taxa(562, intermediate_nodes=True)
ecoli_ids.append(562)

rank_distance = {'species':0,
                 'genus':1,
                 'family':2,
                 'order':3,
                 'class':4,
                 'phylum':5,
                 'kingdom':6,
                 'superkingdom':7,
                 'no rank':8}

class Entry(object):
    """ Simple class to hold information about a parsed PDB entry

    """

    def __init__(self, id=None, org_taxid=None, exp_org_taxid=None, het=None):
        self.id = id
        if org_taxid is None:
            self.org_taxid = set()
        if exp_org_taxid is None:
            self.exp_org_taxid = set()
        if het is None:
            self.het = set()

    def to_dict(self):
        d = {}
        d['id'] = self.id
        d['org_taxid'] = list(self.org_taxid)
        d['exp_org_taxid'] = list(self.exp_org_taxid)
        d['het'] = list(self.het)
        return d

def to_yaml(data, path):
    """ Write a data dictionary to yaml file

    """
    yaml_writer = yaml.YAML()
    yaml_writer.default_flow_style = False
    with open(path, 'wb') as file:
        yaml_writer.dump(data, file)

def read_pdb_yaml():
    """ Read amino acid pdb-ccd from bpforms to a set

    """
    print('reading amino acid set from bpforms')
    yaml_reader = yaml.YAML()
    with open(PROTEIN_YML_PATH, 'rb') as file:
        data = yaml_reader.load(file)

    pdb_monomers = set()
    for monomer in data['monomers'].values():
        for identifier in monomer['identifiers']:
            if identifier['ns'] == 'pdb-ccd':
                pdb_monomers.add(identifier['id'])

    return pdb_monomers

def load_pdb(max_entries):
    """ read the PDB database and parse the entries

    """
    # get amino acid set (canonical and non-canonical) from bpforms
    pdb_monomers = read_pdb_yaml()

    print('loading PDB entries from FTP server')
    # connect to PDB FTP
    ftp = FTP(PDB_FTP_HOST)
    ftp.login()

    ftp.cwd(PDB_ROOT)

    entries_native = []
    entries_engineered_ecoli = []

    # get folders
    total_i = 0
    non_engineered_i = 0
    engineered_ecoli_i = 0
    for dir in ftp.nlst():

        ftp.cwd(os.path.join(PDB_ROOT, dir))
        # get files
        for file_pdb in ftp.nlst():

            # print(file_pdb)
            f = BytesIO()
            ftp.retrbinary('RETR '+file_pdb, f.write)
            f.seek(0)

            # read file
            entry = Entry()
            engineered = False

            # parse file
            with gzip.GzipFile(fileobj=f, mode='rb') as gz:

                line = gz.readline().decode('utf-8').strip()

                while line != '':
                    record_type = line[:6]

                    # get id
                    if record_type == 'HEADER':
                        entry.id = line[62:66]

                    # get only non-engineered protein
                    elif record_type == 'COMPND' and line[11:26] == 'ENGINEERED: YES':
                        engineered = True

                    # get organism
                    elif record_type == 'SOURCE':
                        if line[11:25] == 'ORGANISM_TAXID':
                            if line[-1] == ';':
                                taxid = line[27:-1]
                            else:
                                taxid = line[27:]

                            if taxid.isdigit():
                                entry.org_taxid.add(int(taxid))

                        elif line[11:34] == 'EXPRESSION_SYSTEM_TAXID':
                            # print(file_pdb, line)
                            if line[-1] == ';':
                                exp_taxid = line[36:-1]
                            else:
                                exp_taxid = line[36:]
                            if exp_taxid.isdigit():
                                entry.exp_org_taxid.add(int(exp_taxid))

                    # get heterogen
                    elif record_type == 'HET   ':
                        het = line[7:10].strip()
                        if het in pdb_monomers:
                            entry.het.add(het)

                    line = gz.readline().decode('utf-8').strip()

                gz.close()

            if not engineered:
                # for now, only get entries that are exclusively from one organism
                if len(entry.org_taxid) == 1:
                    # does not include those whose taxid is 'unidentified' or 'synthetic'
                    taxid = list(entry.org_taxid)[0]
                    if taxid != 32630 and taxid != 32644:

                        entries_native.append(entry)
                        # print(entry.id, entry.org_taxid, entry.het)
                        non_engineered_i += 1
            else:
                if len(entry.exp_org_taxid) == 1:
                    taxid = list(entry.exp_org_taxid)[0]
                    if taxid in ecoli_ids:
                        entries_engineered_ecoli.append(entry)
                        engineered_ecoli_i += 1


            total_i += 1
            # print(total_i)

            if total_i % 100 == 0:
                print('progress:', total_i)

        # this block causes early termination of the reading (for testing)

            if isinstance(max_entries, int) and total_i == max_entries:
                break
        else:
            continue
        break


    print('native', non_engineered_i)
    print('engineered in e coli', engineered_ecoli_i)
    print('total examined', total_i)

    return entries_native, entries_engineered_ecoli


def save_entries(entries_native, entries_engineered_ecoli):
    """ save PDB entry parsing results to yaml files

    """

    # save results to a yaml file
    native_yaml_data = []
    for entry in entries_native:
        native_yaml_data.append(entry.to_dict())
        # print(entry.id, entry.org_taxid, entry.het)

    # print(native_yaml_data)
    to_yaml(native_yaml_data, 'pdb_het_native.yml')

    exp_ecoli_yaml_data = []
    for entry in entries_engineered_ecoli:
        exp_ecoli_yaml_data.append(entry.to_dict())
        # print(entry.id, entry.org_taxid, entry.het)

    # print(native_yaml_data)
    to_yaml(exp_ecoli_yaml_data, 'pdb_het_engineered_ecoli.yml')

def reorganize_entries(entries_native):
    """ Reorganize entries by organism

    """
    entries_org = {}
    for entry in entries_native:
        taxid = str(list(entry.org_taxid)[0])
        if taxid not in entries_org:
            entries_org[taxid] = [entry]
        else:
            entries_org[taxid].append(entry)

    return entries_org

def get_ecoli_heterogen(entries_org, entries_engineered_ecoli):
    """ Get the native and full heterogen set of ecoli and descendants

    """

    ecoli_native_het_set = set()
    for ecoli_id in ecoli_ids:
        if str(ecoli_id) in entries_org:
            for ecoli_entry in entries_org[str(ecoli_id)]:
                ecoli_native_het_set.update(ecoli_entry.het)

    ecoli_full_het_set = set()
    for entry in entries_engineered_ecoli:
        ecoli_full_het_set.update(entry.het)
    ecoli_full_het_set = ecoli_full_het_set | ecoli_native_het_set

    return (ecoli_native_het_set, ecoli_full_het_set)

def calc_perc_transformable(entries_org, ecoli_het_set, filename):
    """ Calculate and write percent transformable

    """
    fp = open(filename,'w')
    fp.write('org_id,possible_entries,total_entries,percent_transformable\n')

    df = pd.DataFrame(columns=['org_id','possible_entries','total_entries','percent_transformable'])

    for org_id, org_entries_native in entries_org.items():
        possible_entries = 0
        total_entries = 0
        for entry in org_entries_native:
            total_entries += 1
            if entry.het.issubset(ecoli_het_set):
                possible_entries += 1
        fp.write('{},{},{},{:.4f}\n'.format(org_id, possible_entries, total_entries, possible_entries/total_entries))
        df = df.append(pd.Series([org_id, possible_entries, total_entries, possible_entries/total_entries], index=df.columns), ignore_index=True)

    fp.close()
    df = df.astype({'org_id':'int64', 'possible_entries':'int64', 'total_entries':'int64', 'percent_transformable':'float64' })
    return df

def get_hets_by_entry(data_native, data_engineered_ecoli):

    # get ecoli hets data
    ecoli_hets_by_entry = []
    for entry in data_native:
        if int(entry['org_taxid'][0]) in ecoli_ids:
            # print(int(entry['org_taxid'][0]))
            ecoli_hets_by_entry.append(set(entry['het']))

    print('ecoli native entries', len(ecoli_hets_by_entry))

    # get ecoli hets data
    ecoli_full_hets_by_entry = []
    for entry in data_engineered_ecoli:
        ecoli_full_hets_by_entry.append(set(entry['het']))
    ecoli_full_hets_by_entry.extend(ecoli_hets_by_entry)

    print('ecoli total entries', len(ecoli_full_hets_by_entry))

    return (ecoli_hets_by_entry, ecoli_full_hets_by_entry)


def analyze_rarefaction(ecoli_hets_by_entry, num_iter, filename):
    """ Rarefaction analysis + write figure

    """

    # get all permutations of the list and calculate cumulative discovery
    rarefaction_mat = np.linspace(1, len(ecoli_hets_by_entry), len(ecoli_hets_by_entry)).reshape(-1, 1)

    permutation_list = []

    for i in range(num_iter):
        # print(i)

        random.shuffle(ecoli_hets_by_entry)

        # get rarefaction curve data
        hets_set = set()
        rarefaction_list = []
        for i in range(len(ecoli_hets_by_entry)):
            hets_set = hets_set | ecoli_hets_by_entry[i]
            # print(hets_set)
            rarefaction_list.append(len(hets_set))

        permutation_list.append(rarefaction_list)

    permutation_means = np.mean(np.array(permutation_list), axis=0).reshape(-1,1)

    rarefaction_mat = np.hstack((rarefaction_mat, permutation_means))

    # print(rarefaction_mat)


    plt.figure()
    plt.plot(rarefaction_mat[:,0], rarefaction_mat[:,1])
    plt.savefig(filename, dpi=300)

def analyze_taxonomy(df_perc_transformable, type_suffix):
    print('Analyzing taxonomy of the',type_suffix,'dataset')

    org_ids = df_perc_transformable['org_id'].tolist()

    # get organism scientific names
    org_names = ncbi.get_taxid_translator(org_ids)
    df_perc_transformable['org_name'] = df_perc_transformable['org_id'].map(org_names)


    df_valid = df_perc_transformable[df_perc_transformable['org_name'].notnull()].copy()

    # print(df_valid)
    org_ids_valid = df_valid['org_id'].tolist()

    # get tree
    tree = ncbi.get_topology(org_ids_valid)

    node_ecoli = tree&'562'

    # for each node, find common ancestor between it and the ecoli node.
    # Iterate back from that node back to the root, and get the first node
    # that is not 'no rank'. That rank approximates the distance.

    # The exception is when the two organisms are so far apart that their common
    # ancestor is almost the root of the tree. In that case, the distance rank
    # is 'no rank'

    # cannot handle ids that are merged in ncbi taxonomy database
    errors = []

    common_ancestor_ranks = {}
    for org_id in org_ids_valid:
        try:
            common = node_ecoli.get_common_ancestor(tree&str(org_id))

            if common.rank != 'no rank':
                common_ancestor_ranks[org_id] = common.rank
                continue
            else:
                for ancestor_node in common.iter_ancestors():
                    if ancestor_node.rank != 'no rank':
                        common_ancestor_ranks[org_id] = ancestor_node.rank
                        continue
                common_ancestor_ranks[org_id] = 'no rank'

        except Exception as error:
            errors.append((org_id, error))

    # print(errors)

    # print(common_ancestor_ranks)
    df_valid['ancestor_rank'] = df_valid['org_id'].map(common_ancestor_ranks)
    df_ranksuccess = df_valid[df_valid['ancestor_rank'].notnull()].copy()
    # print(df_ranksuccess)

    # map rank distances
    df_ranksuccess['ancestor_rank_val'] = df_ranksuccess['ancestor_rank'].map(rank_distance)

    # print(df_ranksuccess)

    df_filtered = df_ranksuccess[df_ranksuccess['total_entries']>5]

    # print(df_filtered)

    # plot
    plt.figure()
    plt.scatter(df_filtered['ancestor_rank_val'], df_filtered['percent_transformable'], 10)
    plt.ylim(0.75,1.05)
    # plt.show()
    plt.savefig("perc_vs_taxdist_by_species_"+type_suffix+".png", dpi=300)

    # calculate percentage transformable over each distance

    percent_by_dist_list = []

    for dist in rank_distance.values():
        subset = df_filtered[df_filtered['ancestor_rank_val']==dist]
        subset_sum = subset.sum(axis=0)
        print(dist, subset_sum['total_entries'])
        perc = subset_sum['possible_entries'] / subset_sum['total_entries']
        percent_by_dist_list.append((dist, perc))

    percent_by_dist_mat = np.array(percent_by_dist_list)

    plt.figure()
    plt.scatter(percent_by_dist_mat[:,0], percent_by_dist_mat[:,1], 10)
    plt.ylim(0.75,1.05)
    # plt.show()
    plt.savefig("perc_vs_taxdist_by_dist_"+type_suffix+".png", dpi=300)




def run(max_entries=None):

    (entries_native, entries_engineered_ecoli) = load_pdb(max_entries)

    save_entries(entries_native, entries_engineered_ecoli)

    entries_org = reorganize_entries(entries_native)

    (ecoli_native_het_set, ecoli_full_het_set) = get_ecoli_heterogen(entries_org, entries_engineered_ecoli)
    print('e coli native heterogen set', ecoli_native_het_set)
    print('e coli full heterogen set', ecoli_full_het_set)

    df_perc_transformable_native = calc_perc_transformable(entries_org, ecoli_native_het_set, 'ecoli_transformable_based_on_native.csv')
    df_perc_transformable_full = calc_perc_transformable(entries_org, ecoli_full_het_set, 'ecoli_transformable_based_on_full.csv')

    data_native = []
    for entry in entries_native:
        data_native.append(entry.to_dict())
    data_engineered_ecoli = []
    for entry in entries_engineered_ecoli:
        data_engineered_ecoli.append(entry.to_dict())

    # # if from yaml file, instead of the code above, use the block below:
    # yaml_reader = yaml.YAML()
    # with open(NATIVE_YML_FILE_PATH, 'rb') as file:
    #     data_native = yaml_reader.load(file)
    # with open(ENGINEERED_ECOLI_YML_FILE_PATH, 'rb') as file:
    #     data_engineered_ecoli = yaml_reader.load(file)

    (ecoli_native_hets_by_entry, ecoli_full_hets_by_entry) = get_hets_by_entry(data_native, data_engineered_ecoli)

    analyze_rarefaction(ecoli_native_hets_by_entry, 1000, 'rarefaction_native.png')
    analyze_rarefaction(ecoli_full_hets_by_entry, 10, 'rarefaction_full.png')

    # # if from csv file, use the two lines below
    # df_perc_transformable_native = pd.read_csv(CSV_FILE_PATH_NATIVE)
    # df_perc_transformable_full = pd.read_csv(CSV_FILE_PATH_FULL)

    analyze_taxonomy(df_perc_transformable_native, 'native')
    analyze_taxonomy(df_perc_transformable_native, 'full')

if __name__ == '__main__':
    # run(max_entries=1000)
    run()
