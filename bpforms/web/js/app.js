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

    $("#base_seq_out").val(data['base_seq'])
    $("#length").val(data['length'])
    $("#formula").val(formula)
    $("#mol_wt").val(data['mol_wt'])
    $("#charge").val(data['charge'])
    $("#base_seq_out").css('border-color', '#cacaca');
    $("#base_seq_out").css('border-width', '1px');
    $("#base_seq_out").css('background-color', '#e6e6e6');
}

display_error = function( jqXHR, textStatus, errorThrown ) {
    error = jqXHR['responseJSON']['message']
    if ('details' in jqXHR['responseJSON']) {
        error += '\n\n' + jqXHR['responseJSON']['details']
    }
    $("#base_seq_out").val(error)
    $("#base_seq_out").css('border-color', '#ff0000');
    $("#base_seq_out").css('border-width', '2px');
    $("#base_seq_out").css('background-color', '#ffebeb');
}

$('#submit').click(function (evt) {
    alphabet = $("#alphabet").val()
    base_seq = $("#base_seq_in").val().trim()
    ph = $("#ph").val()

    if (base_seq == null || base_seq == '') {
        return;
    }

    url = '/api/bpform/' + alphabet + '/' + base_seq
    if (ph != null && ph != "") {
         url += '/' + ph
    }

    $.ajax({
      url: url,
      success: set_properties
    })
    .fail(display_error);
});