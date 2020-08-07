rfmake = function(threshold){
  // hide alignment rows below threshold
  var show_up_to = null;
  $("tr").each(function(row_index) {
    $this = $(this);
    var td = $this.find(".bit_score");
    if (td.length > 0) {
      var bit_score = parseFloat($(td).html());
      if (bit_score && bit_score < threshold) {
        $this.hide();
      } else {
        $this.show();
        show_up_to = row_index;
      }
    }
  });

  // hide empty alignment columns
  var alignment_length = $($("tr>td.text-monospace")[1]).find('span').length;
  for (var i = 1; i <= alignment_length; i++) {
    var column = $('tr>td.text-monospace>span:nth-child(' + i + ')');
    var all_gaps = 0;
    for (var j = 2; j < show_up_to; j++) {
      if (column[j].innerHTML.match(/\w/)) {
        all_gaps = 0;
        break;
      }
      all_gaps = 1;
    }
    if (all_gaps) {
      column.hide();
    } else {
      column.show();
    }
  }
};

ralee_conservation = function(threshold){
  reset_ralee_conservation();
  // find rows above threshold
  var show_up_to = null;
  $("tr").each(function(row_index) {
    $this = $(this);
    var td = $this.find(".bit_score");
    if (td.length > 0) {
      var bit_score = parseFloat($(td).html());
      if (bit_score && bit_score >= threshold) {
        show_up_to = row_index;
      } else {
        return false;
      }
    }
  });

  // analyse conservation
  var alignment_length = $($("tr>td.text-monospace")[1]).find('span').length;
  for (var i = 1; i <= alignment_length; i++) {
    var column = $('tr>td.text-monospace>span:nth-child(' + i + ')');
    var guanine = 1, adenine = 1, cytosine = 1, uracil = 1, other = 1;
    for (var j = 2; j <= show_up_to; j++) {
      var nt = column[j].innerHTML;

      if (nt === 'G' || nt === 'g') {
        guanine += 1;
      } else if (nt === 'A' || nt === 'a') {
        adenine += 1;
      } else if (nt === 'C' || nt === 'c') {
        cytosine += 1;
      } else if (nt === 'U' || nt === 'u') {
        uracil += 1;
      } else {
        other += 1;
      }
    }

    adenine = adenine / (show_up_to - other);
    cytosine = cytosine / (show_up_to - other);
    guanine = guanine / (show_up_to - other);
    uracil = uracil / (show_up_to - other);

    var colours = {
      adenine: (adenine >= 0.8) ? 'ralee_80' : (adenine >= 0.6) ? 'ralee_60' : (adenine >= 0.4) ? 'ralee_40' : null,
      cytosine: (cytosine >= 0.8) ? 'ralee_80' : (cytosine >= 0.6) ? 'ralee_60' : (cytosine >= 0.4) ? 'ralee_40' : null,
      guanine: (guanine >= 0.8) ? 'ralee_80' : (guanine >= 0.6) ? 'ralee_60' : (guanine >= 0.4) ? 'ralee_40' : null,
      uracil: (uracil >= 0.8) ? 'ralee_80' : (uracil >= 0.6) ? 'ralee_60' : (uracil >= 0.4) ? 'ralee_40' : null,
    }

    var to_paint = '';
    if (colours.adenine) {
      to_paint += 'Aa';
    }
    if (colours.guanine) {
      to_paint += 'Gg';
    }
    if (colours.cytosine) {
      to_paint += 'Cc';
    }
    if (colours.uracil) {
      to_paint = 'Uu';
    }

    for (var j = 2; j <= show_up_to; j++) {
      var nt = column[j].innerHTML;
      if (to_paint.indexOf(nt) === -1) {
        continue;
      }
      if (nt === 'G' || nt === 'g') {
        column[j].classList.add(colours.guanine);
      } else if (nt === 'A' || nt === 'a') {
        column[j].classList.add(colours.adenine);
      } else if (nt === 'C' || nt === 'c') {
        column[j].classList.add(colours.cytosine);
      } else if (nt === 'U' || nt === 'u') {
        column[j].classList.add(colours.uracil);
      }
    }

  }

};

reset_ralee_conservation = function(){
  $('span.ralee_80, span.ralee_60, span.ralee_40').removeClass('ralee_80 ralee_60 ralee_40');
}

$(function () {
  $('[data-toggle="tooltip"]').tooltip()

  $('#ga_threshold').on('change input', function() {
      $('#ga-label').html($('#ga_threshold').val());
  });

  $('#set-threshold').on('click', function() {
      var old_title = document.title;
      document.title = 'Setting threshold...';
      var threshold = parseFloat($('#ga_threshold').val());
      reset_ralee_conservation();
      rfmake(threshold);
      document.title = old_title;
  });

  $('#reset-threshold').on('click', function() {
      var threshold = 25;
      reset_ralee_conservation();
      rfmake(threshold);
      $('#ga-label').html(threshold);
      $('#ga_threshold').val(threshold);
  });

  $('#ralee-conservation').on('click', function() {
    var threshold = parseFloat($('#ga_threshold').val());
    ralee_conservation(threshold);
  });

  $('#reset-ralee-conservation').on('click', function() {
    reset_ralee_conservation();
  });

  $('.bit_score').on('click', function(){
      var threshold = $(this).html();
      reset_ralee_conservation();
      rfmake(parseFloat(threshold));
      $('#ga-label').html(threshold);
      $('#ga_threshold').val(threshold);
  });

  $('#toggle-columns').on('click', function(e) {
    if ($('.uracil').length) {
      $('.adenine').toggleClass("adenine adenine-inactive");
      $('.guanine').toggleClass("guanine guanine-inactive");
      $('.cytosine').toggleClass("cytosine cytosine-inactive");
      $('.uracil').toggleClass("uracil uracil-inactive");
    } else {
      $('.adenine-inactive').toggleClass("adenine-inactive adenine");
      $('.guanine-inactive').toggleClass("guanine-inactive guanine");
      $('.cytosine-inactive').toggleClass("cytosine-inactive cytosine");
      $('.uracil-inactive').toggleClass("uracil-inactive uracil");
    }
    e.preventDefault();
  });

});
