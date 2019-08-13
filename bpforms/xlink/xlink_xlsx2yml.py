""" Read the manually-curated xlink xlsx file and export it to an xlink ontology

:Author: Mike Zheng <xzheng20@colby.edu>
:Date: 2019-08-05
:Copyright: 2019, Karr Lab
:License: MIT
"""

import pandas as pd
from ruamel import yaml

XLSX_PATH = './xlink.xlsx'
YML_PATH = './xlink.yml'

def main():

    df = pd.read_excel(XLSX_PATH)

    # print(df)

    yml_data = {}
    for index, row in df.iterrows():
        row_dict = {}

        row_dict['synonyms'] = []
        row_dict['synonyms'].append(row['resid_name'])

        if not pd.isnull(row['bond_common_name']):
            row_dict['common_name'] = row['bond_common_name']

        row_dict['l_monomer_alphabet'] = row['left_monomer_alphabet']
        row_dict['l_monomer'] = row['left_monomer']
        row_dict['l_bond_atoms'] = [atom.strip() for atom in row['left_bond_atoms'].split(',')]
        row_dict['l_displaced_atoms'] = [atom.strip() for atom in row['left_displace_atoms'].split(',')]
        row_dict['r_monomer_alphabet'] = row['right_monomer_alphabet']
        row_dict['r_monomer'] = row['right_monomer']
        row_dict['r_bond_atoms'] = [atom.strip() for atom in row['right_bond_atoms'].split(',')]
        row_dict['r_displaced_atoms'] = [atom.strip() for atom in row['right_displace_atoms'].split(',')]

        yml_data[row['xlink_id']] = row_dict

    yaml_writer = yaml.YAML()
    yaml_writer.default_flow_style = False
    with open(YML_PATH, 'wb') as file:
        yaml_writer.dump(yml_data, file)

if __name__ == '__main__':
    main()
