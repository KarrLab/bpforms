$(document).foundation()

jQuery.fn.extend({
 scrollToMe: function () {
   var x = jQuery(this).offset().top - 100;
   jQuery('html,body').animate({scrollTop: x}, 400);
}});

set_alphabets = function(data, status, jqXHR) {
    var html = ''
    for (alphabet_id in data) {
        alphabet = data[alphabet_id]
        html += '<name="' + alphabet['id'] + '"></a>'
        html += '<div class="grid-x grid-padding-x blue" id="' + alphabet['id'] + '">'
        html += '  <div class="large-12 cell">'
        html += '    <h2>' + alphabet['id'] + ': ' + alphabet['name'] + '</h2>'
        html += '    <p>' + alphabet['description'] + '</p>'
        html += '    <table class="alphabet">'
        html += '      <thead>'
        html += '        <tr>'
        html += '          <th>Code</th>'
        html += '          <th>Id</th>'
        html += '          <th>Name</th>'
        html += '          <th>Synonyms</th>'
        html += '          <th>Identifiers</th>'
        html += '          <th>Structure</th>'
        html += '        </tr>'
        html += '      </thead>'
        html += '      <tbody id="alphabet_' + alphabet['id'] + '">'
        html += '      </tbody>'
        html += '    </table>'
        html += '  </div>'
        html += '</div>'
        
        $.ajax({
          url: '/api/alphabet/' + alphabet['id'] + '/',
          success: set_alphabet
        })
        .fail(get_alphabet_error);
    }
    
    $("#alphabets_container").html(html)
    
    scroll_to_hash()
}

get_alphabets_error = function( jqXHR, textStatus, errorThrown ) {
    $("#alphabets_container").html('<span style="color: red">Error loading alphabets</span>')
}

set_alphabet = function(data, status, jqXHR) {
    var html = ''
    
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
    
    $("#alphabet_" + data['id']).html(html)
    
    scroll_to_hash()
}

get_alphabet_error = function( jqXHR, textStatus, errorThrown ) {
    $("#alphabets_container").html('<span style="color: red">Error loading alphabets</span>')
}

 $.ajax({
  url: '/api/alphabet/',
  success: set_alphabets
})
.fail(get_alphabets_error);

scroll_to_hash = function() {
    i_hash = document.URL.indexOf('#')
    if (i_hash >= 0) {
        hash = document.URL.substr(i_hash + 1)
        $("#" + hash).scrollToMe();
    }
}