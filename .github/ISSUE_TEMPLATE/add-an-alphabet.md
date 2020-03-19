---
name: Add an alphabet
about: Add an additional alphabet of residues
title: Add alphabet
labels: alphabet
assignees: ''

---

### Please suggest an id and name for the alphabet

### Please describe the residues in the alphabet
Please provide as much information as possible about each residue
- Id
- Name
- List of synonyms
- List of identifiers (namespace and id)
- Comments
- List of the ides of the originating residues (e.g., the origin of m2A is A)
- Structure in SMILES format
- Bonding sites

Please use the following YAML template
```
id: rna
name: RNA nucleotide monophosphates
description: The canonical RNA nucleotide monophosphates, plus non-canonical RNA nucleotide
  monophosphates based on <a href="http://modomics.genesilico.pl/modifications">MODOMICS</a>
  and the <a href="https://mods.rna.albany.edu/mods/">RNA Modification Database</a>
monomers:
  01A:
    id: 01A
    name: 1,2′-O-dimethyladenosine
    synonyms:
    - 1,2'-O-dimethyladenosine
    - m1Am
    - œ
    identifiers:
    - ns: cas
      id: 91101-00-7
    - ns: rnamods
      id: '97'
    - ns: modomics.short_name
      id: m1Am
    - ns: modomics.new_nomenclature
      id: 01A
    comments: <p>This compound is also listed under CA registry number 59867-24-2
      as the monohydriodide salt, reflecting the form in which the first reported
      synthesis was carried out.</p> <p>Phylogenetic distribution<ul><li>eukarya tRNA</li></ul></p>  
    base_monomers:
    - A
    structure: COC1C(O)C(OC1n1cnc2c1ncn(c2=N)C)COP(=O)([O-])[O-]
    r_bond_atoms:
    - molecule: Monomer
      element: O
      position: 5
    l_bond_atoms:
    - molecule: Monomer
      element: P
      position: 22
    r_displaced_atoms:
    - molecule: Monomer
      element: H
      position: 5
    l_displaced_atoms:
    - molecule: Monomer
      element: O
      position: 25
      charge: -1
```
