$(document).foundation()

set_alphabet = function(data, status, jqXHR) {
    var html = '<h2>' + data['name'] + ' alphabet</h2>'
    html += '<p>' + data['description'] + '</p>'

    html += '<table class="alphabet">'
    html += '  <thead>'
    html += '    <tr>'
    html += '      <th>Code</th>'
    html += '      <th>Id</th>'
    html += '      <th>Name</th>'
    html += '      <th>Synonyms</th>'
    html += '      <th>Identifiers</th>'
    html += '      <th>Structure</th>'
    html += '    </tr>'
    html += '  </thead>'
    html += '  <tbody>'

    for (chars in data['monomers']) {
        monomer = data['monomers'][chars]
        html += '<tr>'
        html += '<td class="chars"><div>' + chars + '</div></td>'

        id = monomer['id']
        if (id == null)
            id = ''
        html += '<td class="id"><div>' + id + '</div></td>'

        name = monomer['name']
        if (name == null)
            name = ''
        html += '<td class="name"><div>' + name + '</div></td>'

        synonyms = monomer['synonyms']
        if (synonyms == null)
            synonyms = ''
        else
            synonyms = '<ul><li>' + synonyms.join('</li><li>') + '</li></ul>'
        html += '<td class="synonyms"><div>' + synonyms + '</div></td>'


        identifiers = monomer['identifiers']
        if (identifiers == null)
            identifiers = ''
        else{
            identifiers_strs = '<ul>'
            for (i in identifiers) {
                identifiers_strs += '<li>' + identifiers[i]['ns'] + '/' + identifiers[i]['id'] + '</li>'
            }
            identifiers_strs += '</ul>'
            identifiers = identifiers_strs
        }
        html += '<td class="identifiers"><div>' + identifiers + '</div></td>'

        structure = monomer['structure']
        if (structure == null)
            structure = ''
        html += '<td class="structure"><div>' + structure + '</div></td>'

        html += '</tr>'
    }
    html += '  </tbody>'
    html += '</table>'

    $("#alphabet_container").html(html)
}

get_alphabet_error = function( jqXHR, textStatus, errorThrown ) {
    $("#alphabet_container").html('<span style="color: red">Error loading alphabet</span>')
}

get_alphabet = function() {
    i_hash = document.URL.indexOf('#')
    if (i_hash == -1) {
        console.log('Hash must be the id of an alphabet')
    }else{
        hash = document.URL.substr(i_hash + 1)
        $.ajax({
          url: '/api/alphabet/' + hash,
          success: set_alphabet
        })
        .fail(get_alphabet_error);
    }
}

get_alphabet()
window.onhashchange = get_alphabet