$(document).foundation()

set_crosslink = function(data, status, jqXHR) {

    var html = ''

    codes = Object.keys(data)
    codes.sort()

    for (i in codes) {
        code = codes[i]
        xlink = data[code]
        html += '<div class="cell monomer">'

        // code
        html += '<div class="code">' + code + '</div>'

        // structure
        img = '<img class="lazy" data-src="img/crosslink/' + code + '.png">'
        html += '<div class="structure">' + img + '</div>'

        // show/hide details
        html += '<div class="toggle" id="monomer_' + i.toString() + '_toggle">'
        html += '<a onclick="javascript: toggle_monomer_details(\'' + i.toString() + '\')">Show details</a>'
        html += '</div>'

        /////////////////////
        // details
        /////////////////////
        html += '<div id="monomer_' + i.toString() + '_details" class="details hide">'
        html += '<table width="100%"><tbody>'

        // id
        html += '<tr><th>Id</th><td>' + code + '</td></tr>'

        // name
        if (xlink['common_name'] != null)
            html += '<tr><th>Common Name</th><td>' + xlink['common_name'].replace(/_/g, ' ') + '</td></tr>'

        // synonyms
        if (xlink['synonyms'] != null)
            html += '<tr><th>Synonyms</th><td><p>' + xlink['synonyms'].join('</p><p>') + '</p></td></tr>'

        // l_monomer
        if (xlink['l_monomer'] != null)
            html += '<tr><th>Left residue</th><td><p>' + xlink['l_monomer'] + '</p></td></tr>'

        // r_monomer
        if (xlink['r_monomer'] != null)
            html += '<tr><th>Right residue</th><td><p>' + xlink['r_monomer'] + '</p></td></tr>'

        html += '</tbody></table>'
        html += '</div>'
        html += '</div>'
    }

    $("#crosslink_container").html(html)
    $(function() {
        $('.lazy').Lazy();
    });

}

get_crosslink_error = function( jqXHR, textStatus, errorThrown ) {
    $("#crosslink_container").html('<span style="color: red">Error loading crosslinks</span>')
}

get_crosslink = function() {
    $.ajax({
      url: 'data/xlink.json',
      dataType: 'json',
      success: set_crosslink
    })
    .fail(get_crosslink_error);
}

get_crosslink()

toggle_monomer_details = function(i_monomer) {
    details = $('#monomer_' + i_monomer + '_details')
    details.toggleClass('hide');

    if (details.hasClass('hide'))
        $('#monomer_' + i_monomer + '_toggle>a').html('Show details');
    else
        $('#monomer_' + i_monomer + '_toggle>a').html('Hide details');
}
