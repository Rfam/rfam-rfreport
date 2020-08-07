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
  for (i = 1; i <= alignment_length; i++) {
    var column = $('tr>td.text-monospace>span:nth-child(' + i + ')');
    var all_gaps = 0;
    for (j = 0; j <= show_up_to; j++) {
      if (column[j].innerHTML.match(/\w/)) {
        all_gaps = 0;
        break;
      }
      all_gaps = 1;
    }
    if (all_gaps) {
      var a = 1;
      // column.hide();
    } else {
      var a = 1;
      // column.show();
    }
  }
};

$(function () {
  $('[data-toggle="tooltip"]').tooltip()

  $('#ga_threshold').on('change input', function() {
      $('#ga-label').html($('#ga_threshold').val());
  });

  $('#set-threshold').on('click', function() {
      rfmake(parseFloat($('#ga_threshold').val()));
  });

  $('#reset-threshold').on('click', function() {
      var threshold = 25;
      rfmake(threshold);
      $('#ga-label').html(threshold);
      $('#ga_threshold').val(threshold);
  });

  $('.bit_score').on('click', function(){
      var threshold = $(this).html();
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
