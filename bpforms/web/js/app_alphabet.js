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
        if (structure == null) {
            img = ''
        } else {
            // https://cactus.nci.nih.gov/blog/?tag=png
            img = '<img src="https://cactus.nci.nih.gov/chemical/structure/' 
                  + encodeURI(structure)
                  + '/image?format=gif&bgcolor=transparent&antialiasing=0" class="context-menu-one" structure="' 
                  + encodeURI(structure) 
                  + '"/>'
        }
        html += '<td class="structure"><div>' + img + '</div></td>'

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