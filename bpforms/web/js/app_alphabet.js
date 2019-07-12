$(document).foundation()

set_alphabet = function(data, status, jqXHR) {
    i_hash = document.URL.indexOf('#')
    alphabet_id = document.URL.substr(i_hash + 1)

    $('#alphabet_heading').html(data['name'] + ' alphabet')
    $('#alphabet_description').html(data['description'])

    var html = ''

    codes = Object.keys(data['monomers'])
    codes.sort()
    for (i_monomer in codes) {
        code = codes[i_monomer]
        monomer = data['monomers'][code]
        html += '<div class="cell monomer">'

        // code
        html += '<div class="code">' + code + '</div>'

        // structure
        if (monomer['structure'] == null) {
            img = ''
        } else {
            img = '<img class="lazy" data-src="img/alphabet/' + alphabet_id + '/' + code + '.png">'
        }
        html += '<div class="structure">' + img + '</div>'

        // show/hide details
        html += '<div class="toggle" id="monomer_' + i_monomer.toString() + '_toggle">'
        html += '<a onclick="javascript: toggle_monomer_details(\'' + i_monomer.toString() + '\')">Show details</a>'
        html += '</div>'

        /////////////////////
        // details
        /////////////////////
        html += '<div id="monomer_' + i_monomer.toString() + '_details" class="details hide">'
        html += '<table width="100%"><tbody>'

        // id
        if (monomer['id'] != null)
            html += '<tr><th>Id</th><td>' + monomer['id'] + '</td></tr>'

        // name
        if (monomer['name'] != null)
            html += '<tr><th>Name</th><td>' + monomer['name'] + '</td></tr>'

        // synonyms
        if (monomer['synonyms'] != null)
            html += '<tr><th>Synonyms</th><td><p>' + monomer['synonyms'].join('</p><p>') + '</p></td></tr>'

        // identifiers
        identifiers = monomer['identifiers']
        if (monomer['identifiers'] != null) {
            identifiers_str = ''
            for (i in identifiers) {
                ns_long = identifiers[i]['ns']
                switch (identifiers[i]['ns']){
                    case 'cas':
                        ns_short = 'CAS'
                        url = 'https://www.chemicalbook.com/Search_EN.aspx?keyword=' + identifiers[i]['id']
                        break;
                    case 'chebi':
                        ns_short = 'ChEBI'
                        url = 'https://www.ebi.ac.uk/chebi/searchId.do?chebiId=' + identifiers[i]['id']
                        break;
                    case 'dnamod':
                        ns_short = 'DNAmod'
                        url = 'https://dnamod.hoffmanlab.org/' + identifiers[i]['id'] + '.html'
                        break;
                    case 'go':
                        ns_short = 'GO'
                        url = 'http://amigo.geneontology.org/amigo/term/' + identifiers[i]['id']
                        break;
                    case 'metacyc.compound':
                        ns_short = 'MetaCyc'
                        url = 'https://biocyc.org/compound?orgid=META&id=' + identifiers[i]['id']
                        break;
                    case 'mod':
                        ns_short = 'MOD'
                        url = 'https://www.ebi.ac.uk/ols/ontologies/mod/terms?obo_id=' + identifiers[i]['id']
                        break;
                    case 'modomics.short_name':
                        ns_short = 'MODOMICS'
                        url = 'http://modomics.genesilico.pl/modifications/' + identifiers[i]['id']
                        break;
                    case 'modomics.new_nomenclature':
                        ns_short = 'MODOMICS'
                        url = ''
                        break;
                    case 'pdb.ligand':
                        ns_short = 'PDB'
                        url = 'http://www.rcsb.org/ligand/' + identifiers[i]['id']
                        break;
                    case 'pdb-ccd':
                        ns_short = 'PDB'
                        url = 'http://www.rcsb.org/ligand/' + identifiers[i]['id']
                        break;
                    case 'pubchem.compound':
                        ns_short = 'PubChem'
                        url = 'https://pubchem.ncbi.nlm.nih.gov/compound/' + identifiers[i]['id']
                        break;
                    case 'resid':
                        ns_short = 'RESID'
                        url = 'https://annotation.dbi.udel.edu/cgi-bin/resid?id=' + identifiers[i]['id']
                        break;

                }

                if (url)
                    identifiers_str += '<p>' + ns_short + ': <a href="' + url + '">' + identifiers[i]['id'] + '</a></p>'
                else
                    identifiers_str += '<p>' + ns_short + ': ' + identifiers[i]['id'] + '</p>'
            }
            html += '<tr><th>Links</th><td>' + identifiers_str + '</td></tr>'
        }

        // base
        if (monomer['base_monomers'] != null)
            html += '<tr><th>Base<br/>monomer(s)</th><td>' + monomer['base_monomers'].join(', ') + '</td></tr>'

        // binds backbone?
        html += '<tr><th>Binds<br/>backbone?</th><td>' + (monomer['binds_backbone'] ? 'Yes' : 'No') + '</td></tr>'

        // binds left?
        html += '<tr><th>Binds<br/>left?</th><td>' + (monomer['binds_left'] ? 'Yes' : 'No') + '</td></tr>'

        // binds right?
        html += '<tr><th>Binds<br/>right?</th><td>' + (monomer['binds_right'] ? 'Yes' : 'No') + '</td></tr>'

        // structure
        html += '<tr><th>Structure<br/>(SMILES)</th><td>' + monomer['structure'] + '</td></tr>'

        // delta mass
        if (monomer['delta_mass'])
            html += '<tr><th>Mass<br/>uncertainty</th><td>' + monomer['delta_mass'] + '</td></tr>'

        // delta charge
        if (monomer['delta_charge'])
            html += '<tr><th>Charge<br/>uncertainity</th><td>' + monomer['delta_charge'] + '</td></tr>'

        // start/end position
        if (monomer['start_position'] || monomer['end_position'])
            html += '<tr><th>Position<br/>uncertainty</th><td>' + monomer['start_position'] + '-' + monomer['end_position'] + '</td></tr>'

        // formula
        formula = ''
        for (el in monomer['formula'])
            formula += el + '<sub>' + monomer['formula'][el] + '</sub>'
        html += '<tr><th>Formula</th><td>' + formula + '</td></tr>'

        // molecular weight
        html += '<tr><th>Mass</th><td>' + (Math.round(monomer['mol_wt'] * 1000) / 1000) + '</td></tr>'

        // charge
        html += '<tr><th>Charge</th><td>' + monomer['charge'] + '</td></tr>'

        // comments
        if (monomer['comments'])
            html += '<tr><th>Comments</th><td>' + monomer['comments'] + '</td></tr>'

        html += '</tbody></table>'
        html += '</div>'
        html += '</div>'
    }

    $("#alphabet_container").html(html)
    $(function() {
        $('.lazy').Lazy();
    });
}

get_alphabet_error = function( jqXHR, textStatus, errorThrown ) {
    $("#alphabet_container").html('<span style="color: red">Error loading alphabet</span>')
}

get_alphabet = function() {
    i_hash = document.URL.indexOf('#')
    if (i_hash == -1) {
    }else{
        hash = document.URL.substr(i_hash + 1)
        $.ajax({
          url: 'data/alphabet/' + hash + '.json',
          dataType: 'json',
          success: set_alphabet
        })
        .fail(get_alphabet_error);
    }
}

get_alphabet()
window.onhashchange = get_alphabet

toggle_monomer_details = function(i_monomer) {
    details = $('#monomer_' + i_monomer + '_details')
    details.toggleClass('hide');

    if (details.hasClass('hide'))
        $('#monomer_' + i_monomer + '_toggle>a').html('Show details');
    else
        $('#monomer_' + i_monomer + '_toggle>a').html('Hide details');
}

// https://swisnl.github.io/jQuery-contextMenu
$(function(){
    $.contextMenu({
        selector: '.context-menu-one',
        items: {
            key: {
                name: "Copy as SMILES",
                callback: function(data, status, jqXHR) {
                    structure = decodeURI($(status.$trigger[0]).attr('structure'))
                    copyToClipboard(structure)
                }
            }
        },
    });
});

// https://hackernoon.com/copying-text-to-clipboard-with-javascript-df4d4988697f
const copyToClipboard = str => {
  const el = document.createElement('textarea');  // Create a <textarea> element
  el.value = str;                                 // Set its value to the string that you want copied
  el.setAttribute('readonly', '');                // Make it readonly to be tamper-proof
  el.style.position = 'absolute';
  el.style.left = '-9999px';                      // Move outside the screen to make it invisible
  document.body.appendChild(el);                  // Append the <textarea> element to the HTML document
  const selected =
    document.getSelection().rangeCount > 0        // Check if there is any content selected previously
      ? document.getSelection().getRangeAt(0)     // Store selection if found
      : false;                                    // Mark as false to know no selection existed before
  el.select();                                    // Select the <textarea> content
  document.execCommand('copy');                   // Copy - only works as a result of a user action (e.g. click events)
  document.body.removeChild(el);                  // Remove the <textarea> element
  if (selected) {                                 // If a selection existed before copying
    document.getSelection().removeAllRanges();    // Unselect everything on the HTML document
    document.getSelection().addRange(selected);   // Restore the original selection
  }
};