$(document).foundation()

set_alphabets = function(data, status, jqXHR) {
    for (alphabet_id in data) {
        alphabet = data[alphabet_id]
        $("#alphabet").append('<option value="' + alphabet['id'] + '">' + alphabet['name'] + '</option>')
        $("#alphabet_descriptions").append('<li><a href="alphabet.html#' + alphabet['id'] + '"><i>' + alphabet['name'] + '</i></a>: ' + alphabet['description'] + '</li>')
    }
}

get_alphabets_error = function( jqXHR, textStatus, errorThrown ) {
    $("#alphabet_descriptions").html('<span style="color: red">Error loading alphabets</span>')
}

 $.ajax({
  url: '/api/alphabet/',
  success: set_alphabets
})
.fail(get_alphabets_error);

set_properties = function(data, status, jqXHR){
    formula = ''
    for (el in data['formula']) {
        formula += el + data['formula'][el].toString()
    }

    $("#errors").html('')
    $("#monomer_seq_out").val(data['monomer_seq'])
    $("#length").val(data['length'])

    console.log(data['structure'])
    if (data['structure'] != null && data['structure'] != '') {
        structure = data['structure']
        img = '<img src="https://cactus.nci.nih.gov/chemical/structure/'
              + encodeURI(data['structure'])
              + '/image?format=gif&bgcolor=transparent&antialiasing=0" class="context-menu-one"'
              + '/>'

    } else {
        structure = ''
        img = ''
    }
    $("#structure").val(structure)
    $("#structure-img").html(img)
    $("#formula").val(formula)
    $("#mol_wt").val(data['mol_wt'])
    $("#charge").val(data['charge'])
}

display_error = function( jqXHR, textStatus, errorThrown ) {
    error = '<b>' + jqXHR['responseJSON']['message'] + '</b>'
    if ('errors' in jqXHR['responseJSON']) {
        error += '<ul>'
        for (field in jqXHR['responseJSON']['errors'])
            error += '<li>'
            if (field != '')
                error += '<span style="text-decoration: underline;">' + field + '</span>: '
            error += jqXHR['responseJSON']['errors'][field]
            error += '</li>'
        error += '</ul>'
    }
    console.log(jqXHR['responseJSON'])
    $("#errors").html(error)
}

$('#submit').click(function (evt) {
    alphabet = $("#alphabet").val()
    monomer_seq = $("#monomer_seq_in").val().trim()
    circular = $("#circular").val()
    ph = $("#ph").val()
    major_tautomer = $("#major_tautomer").val()

    if (monomer_seq == null || monomer_seq == '') {
        return;
    }

    data = {
        'alphabet': alphabet,
        'monomer_seq': monomer_seq
    }
    if (circular != null && circular != '') {
        data['circular'] = circular == '1'
    }
    if (ph != null && ph != '') {
        data['ph'] = parseFloat(ph)
    }
    if (major_tautomer != null && major_tautomer != '') {
        data['major_tautomer'] = major_tautomer == '1'
    }

    $.ajax({
      type: 'post',
      url: '/api/bpform/',
      data: JSON.stringify(data),
      contentType : 'application/json',
      dataType: 'json',
      success: set_properties
    })
    .fail(display_error);
});
