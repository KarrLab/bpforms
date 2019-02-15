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

    $("#monomer_seq_out").val(data['monomer_seq'])
    $("#length").val(data['length'])
    $("#formula").val(formula)
    $("#mol_wt").val(data['mol_wt'])
    $("#charge").val(data['charge'])
    $("#monomer_seq_out").css('border-color', '#cacaca');
    $("#monomer_seq_out").css('border-width', '1px');
    $("#monomer_seq_out").css('background-color', '#e6e6e6');
}

display_error = function( jqXHR, textStatus, errorThrown ) {
    error = jqXHR['responseJSON']['message']
    if ('details' in jqXHR['responseJSON']) {
        error += '\n\n' + jqXHR['responseJSON']['details']
    }
    $("#monomer_seq_out").val(error)
    $("#monomer_seq_out").css('border-color', '#ff0000');
    $("#monomer_seq_out").css('border-width', '2px');
    $("#monomer_seq_out").css('background-color', '#ffebeb');
}

$('#submit').click(function (evt) {
    alphabet = $("#alphabet").val()
    monomer_seq = $("#monomer_seq_in").val().trim()
    ph = $("#ph").val()
    major_tautomer = $("#major_tautomer").val()

    if (monomer_seq == null || monomer_seq == '') {
        return;
    }

    url = '/api/bpform/' + alphabet + '/' + monomer_seq
    if (ph != null && ph != "") {
         url += '/' + ph
         url += '/' + major_tautomer
    }

    $.ajax({
      url: url,
      success: set_properties
    })
    .fail(display_error);
});