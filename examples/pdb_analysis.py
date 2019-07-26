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
LOCAL_PDB_PATH = '/root/karrlab/pdb/'

# Escherichia coli
query_org_id = 562

# Gammaproteobacteria
# query_org_id = 1236

# Saccharomyces cerevisiae
# query_org_id = 4932

# Saccharomycetes
# query_org_id = 4891

# get org_ids for query organism and descendants
ncbi = NCBITaxa()
query_ids = ncbi.get_descendant_taxa(query_org_id, intermediate_nodes=True)
query_ids.append(query_org_id)

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

def load_full_pdb_hets_by_entry_from_ftp(max_entries):
    """ read the PDB database and get heterogen sets for all entries

    """
    # get amino acid set (canonical and non-canonical) from bpforms
    pdb_monomers = read_pdb_yaml()

    print('loading PDB entries from FTP server')
    # connect to PDB FTP
    ftp = FTP(PDB_FTP_HOST)
    ftp.login()

    ftp.cwd(PDB_ROOT)

    hets_by_entry = []
    total_i = 0


    for dir in ftp.nlst():

        ftp.cwd(os.path.join(PDB_ROOT, dir))
        # get files
        for file_pdb in ftp.nlst():

            # print(file_pdb)
            f = BytesIO()
            ftp.retrbinary('RETR '+file_pdb, f.write)
            f.seek(0)

            with gzip.GzipFile(fileobj=f, mode='rb') as gz:

                het_set = set()

                line = gz.readline().decode('utf-8').strip()

                while line != '':
                    record_type = line[:6]

                    if record_type == 'HET   ':
                        het = line[7:10].strip()
                        if het in pdb_monomers:
                            het_set.add(het)

                    line = gz.readline().decode('utf-8').strip()

                hets_by_entry.append(het_set)
                gz.close()

            total_i += 1
            # print(total_i)

            if total_i % 100 == 0:
                print('progress:', total_i)

            if isinstance(max_entries, int) and total_i == max_entries:
                break
        else:
            continue
        break

    ftp.close()
    print('total examined', total_i)

    return hets_by_entry


def load_full_pdb_hets_by_entry_from_local(max_entries):
    """ read the PDB database and get heterogen sets for all entries from local directory

    """
    # get amino acid set (canonical and non-canonical) from bpforms
    pdb_monomers = read_pdb_yaml()

    print('loading PDB entries from local directory', LOCAL_PDB_PATH)
    # connect to PDB FTP

    hets_by_entry = []
    total_i = 0

    for root, directories, filenames in os.walk(LOCAL_PDB_PATH):

        for filename in filenames:
            if filename[0] != '.':
                f = os.path.join(root,filename)

                with gzip.open(f, 'rb') as gz:

                    het_set = set()

                    line = gz.readline().decode('utf-8').strip()

                    while line != '':
                        record_type = line[:6]

                        if record_type == 'HET   ':
                            het = line[7:10].strip()
                            if het in pdb_monomers:
                                het_set.add(het)

                        line = gz.readline().decode('utf-8').strip()

                    hets_by_entry.append(het_set)
                    gz.close()

                total_i += 1
                # print(total_i)

                if total_i % 100 == 0:
                    print('progress:', total_i)

                if isinstance(max_entries, int) and total_i == max_entries:
                    break
        else:
            continue
        break

    print('total examined', total_i)

    return hets_by_entry


def load_pdb_from_ftp(max_entries):
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
    entries_engineered_in_query = []

    # get folders
    total_i = 0
    non_engineered_i = 0
    engineered_in_query_i = 0
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
                    if taxid in query_ids:
                        entries_engineered_in_query.append(entry)
                        engineered_in_query_i += 1


            total_i += 1
            # print(total_i)

            if total_i % 100 == 0:
                print('progress:', total_i)

            if isinstance(max_entries, int) and total_i == max_entries:
                break
        else:
            continue
        break

    ftp.close()

    print('native', non_engineered_i)
    print('engineered in query', engineered_in_query_i)
    print('total examined', total_i)

    return entries_native, entries_engineered_in_query

def load_pdb_from_local(max_entries):
    """ read the PDB database and parse the entries from local directory

    """
    # get amino acid set (canonical and non-canonical) from bpforms
    pdb_monomers = read_pdb_yaml()

    print('loading PDB entries from local directory', LOCAL_PDB_PATH)
    # connect to PDB FTP

    entries_native = []
    entries_engineered_in_query = []

    # get folders
    total_i = 0
    non_engineered_i = 0
    engineered_in_query_i = 0

    for root, directories, filenames in os.walk(LOCAL_PDB_PATH):

        for filename in filenames:
            if filename[0] != '.':
                f = os.path.join(root,filename)

                # read file
                entry = Entry()
                engineered = False

                # parse file
                with gzip.open(f, 'rb') as gz:

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
                        if taxid in query_ids:
                            entries_engineered_in_query.append(entry)
                            engineered_in_query_i += 1


                total_i += 1
                # print(total_i)

                if total_i % 100 == 0:
                    print('progress:', total_i)

                if isinstance(max_entries, int) and total_i == max_entries:
                    break
        else:
            continue
        break


    print('native', non_engineered_i)
    print('engineered in query', engineered_in_query_i)
    print('total examined', total_i)

    return entries_native, entries_engineered_in_query




def save_entries(entries_native, entries_engineered_in_query):
    """ save PDB entry parsing results to yaml files

    """

    # save results to a yaml file
    native_yaml_data = []
    for entry in entries_native:
        native_yaml_data.append(entry.to_dict())
        # print(entry.id, entry.org_taxid, entry.het)

    # print(native_yaml_data)
    to_yaml(native_yaml_data, 'pdb_het_native_of_'+str(query_org_id)+'.yml')

    exp_query_yaml_data = []
    for entry in entries_engineered_in_query:
        exp_query_yaml_data.append(entry.to_dict())
        # print(entry.id, entry.org_taxid, entry.het)

    # print(native_yaml_data)
    to_yaml(exp_query_yaml_data, 'pdb_het_engineered_in_'+str(query_org_id)+'.yml')

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

def get_query_heterogen(entries_org, entries_engineered_in_query):
    """ Get the native and full heterogen set of the query organism and descendants

    """

    query_native_het_set = set()
    for query_id in query_ids:
        if str(query_id) in entries_org:
            for query_entry in entries_org[str(query_id)]:
                query_native_het_set.update(query_entry.het)

    query_full_het_set = set()
    for entry in entries_engineered_in_query:
        query_full_het_set.update(entry.het)
    query_full_het_set = query_full_het_set | query_native_het_set

    return (query_native_het_set, query_full_het_set)

def calc_perc_transformable(entries_org, query_het_set, out_filename):
    """ Calculate and write percent transformable

    """
    fp = open(out_filename,'w')
    fp.write('org_id,possible_entries,total_entries,percent_transformable\n')

    df = pd.DataFrame(columns=['org_id','possible_entries','total_entries','percent_transformable'])

    for org_id, org_entries_native in entries_org.items():
        possible_entries = 0
        total_entries = 0
        for entry in org_entries_native:
            total_entries += 1
            if entry.het.issubset(query_het_set):
                possible_entries += 1
        fp.write('{},{},{},{:.4f}\n'.format(org_id, possible_entries, total_entries, possible_entries/total_entries))
        df = df.append(pd.Series([org_id, possible_entries, total_entries, possible_entries/total_entries], index=df.columns), ignore_index=True)

    fp.close()
    df = df.astype({'org_id':'int64', 'possible_entries':'int64', 'total_entries':'int64', 'percent_transformable':'float64' })
    return df

def get_hets_by_entry(data_native, data_engineered_in_query):

    # get query hets data
    query_hets_by_entry = []
    for entry in data_native:
        if int(entry['org_taxid'][0]) in query_ids:
            # print(int(entry['org_taxid'][0]))
            query_hets_by_entry.append(set(entry['het']))

    print('query native entries', len(query_hets_by_entry))

    # get query hets data
    query_full_hets_by_entry = []
    for entry in data_engineered_in_query:
        query_full_hets_by_entry.append(set(entry['het']))
    query_full_hets_by_entry.extend(query_hets_by_entry)

    print('query total entries', len(query_full_hets_by_entry))

    return (query_hets_by_entry, query_full_hets_by_entry)


def analyze_rarefaction(query_hets_by_entry, num_iter, out_filename):
    """ Rarefaction analysis + write figure

    """

    # get all permutations of the list and calculate cumulative discovery
    rarefaction_mat = np.linspace(1, len(query_hets_by_entry), len(query_hets_by_entry)).reshape(-1, 1)

    permutation_list = []

    for i in range(num_iter):
        # print(i)

        random.shuffle(query_hets_by_entry)

        # get rarefaction curve data
        hets_set = set()
        rarefaction_list = []
        for i in range(len(query_hets_by_entry)):
            hets_set = hets_set | query_hets_by_entry[i]
            # print(hets_set)
            rarefaction_list.append(len(hets_set))

        permutation_list.append(rarefaction_list)

    permutation_means = np.mean(np.array(permutation_list), axis=0).reshape(-1,1)

    rarefaction_mat = np.hstack((rarefaction_mat, permutation_means))

    # print(rarefaction_mat)


    plt.figure()
    plt.plot(rarefaction_mat[:,0], rarefaction_mat[:,1])
    plt.xlabel('# of PDB entries')
    plt.ylabel('# of n.c. aa')
    plt.savefig(out_filename, dpi=300)

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
    tree = ncbi.get_topology(org_ids_valid, intermediate_nodes=True)

    node_query = tree&(str(query_org_id))

    # for each node, find common ancestor between it and the query organism node.
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
            common = node_query.get_common_ancestor(tree&str(org_id))

            if common.rank != 'no rank':
                common_ancestor_ranks[org_id] = common.rank
            else:
                for ancestor_node in common.iter_ancestors():
                    if ancestor_node.rank != 'no rank':
                        common_ancestor_ranks[org_id] = ancestor_node.rank
                        break
                else:
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

    # df_filtered = df_ranksuccess[df_ranksuccess['total_entries']>5]
    df_filtered = df_ranksuccess.copy()

    # print(df_filtered)

    # # plot
    # plt.figure()
    # plt.scatter(df_filtered['ancestor_rank_val'], df_filtered['percent_transformable'], 10)
    # plt.ylim(0.75,1.05)
    # # plt.show()
    # plt.savefig('perc_vs_taxdist_by_species_'+type_suffix+'_'+str(query_org_id)+'.png', dpi=300)

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
    plt.xlabel('Taxonomic Distance')
    plt.ylabel('Transformable')
    plt.ylim(0.7,1.05)
    # plt.show()
    plt.savefig('perc_vs_taxdist_by_dist_'+type_suffix+'_'+str(query_org_id)+'.png', dpi=300)


def run_analyze_full_pdb_rarefaction(max_entries=None, use_local=False):

    print('rarefaction analysis of full PDB database')

    if use_local:
        hets_by_entry = load_full_pdb_hets_by_entry_from_local(max_entries)
    else:
        hets_by_entry = load_full_pdb_hets_by_entry_from_ftp(max_entries)

    analyze_rarefaction(hets_by_entry, 10, 'rarefaction_pdb.png')



def run_analyze_org(max_entries=None, use_local=False):

    print('analyze organism', query_org_id)

    if use_local:
        (entries_native, entries_engineered_in_query) = load_pdb_from_local(max_entries)
    else:
        (entries_native, entries_engineered_in_query) = load_pdb_from_ftp(max_entries)

    save_entries(entries_native, entries_engineered_in_query)

    entries_org = reorganize_entries(entries_native)

    (query_native_het_set, query_full_het_set) = get_query_heterogen(entries_org, entries_engineered_in_query)
    print('query organism native heterogen set', query_native_het_set)
    print('query organism full heterogen set', query_full_het_set)

    df_perc_transformable_native = calc_perc_transformable(entries_org, query_native_het_set, 'query_transformable_based_on_native_'+str(query_org_id)+'.csv')
    df_perc_transformable_full = calc_perc_transformable(entries_org, query_full_het_set, 'query_transformable_based_on_full_'+str(query_org_id)+'.csv')

    data_native = []
    for entry in entries_native:
        data_native.append(entry.to_dict())
    data_engineered_in_query = []
    for entry in entries_engineered_in_query:
        data_engineered_in_query.append(entry.to_dict())

    # # if from yaml file, instead of the code above, use the block below:
    # yaml_reader = yaml.YAML()
    # with open(NATIVE_YML_FILE_PATH, 'rb') as file:
    #     data_native = yaml_reader.load(file)
    # with open(ENGINEERED_IN_QUERY_YML_FILE_PATH, 'rb') as file:
    #     data_engineered_in_query = yaml_reader.load(file)

    (query_native_hets_by_entry, query_full_hets_by_entry) = get_hets_by_entry(data_native, data_engineered_in_query)

    analyze_rarefaction(query_native_hets_by_entry, 1000, 'rarefaction_native_'+str(query_org_id)+'.png')
    analyze_rarefaction(query_full_hets_by_entry, 10, 'rarefaction_full_'+str(query_org_id)+'.png')

    # # if from csv file, use the two lines below
    # df_perc_transformable_native = pd.read_csv(CSV_FILE_PATH_NATIVE)
    # df_perc_transformable_full = pd.read_csv(CSV_FILE_PATH_FULL)

    analyze_taxonomy(df_perc_transformable_native, 'native')
    analyze_taxonomy(df_perc_transformable_full, 'full')

if __name__ == '__main__':
    # run_analyze_org(max_entries=1000)
    run_analyze_full_pdb_rarefaction(max_entries=1000)

    # run_analyze_org(use_local=True)
    # run_analyze_full_pdb_rarefaction(use_local=True)
