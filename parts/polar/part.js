// Ring of test masses beside the line segment that stands for the same wave.
// Stage 0: the ring alone. Stage 1: the segment appears, on the same clock.
// The segment is the "stretch whisker" of the workshop visualizer: a headless bar centred
// on the point, along psi = 0.5 atan2(h_x, h_+), with length proportional to hypot(h_+, h_x).
(function () {
  "use strict";

  // Constants of the gw-ring part, which this part continues.
  var R0 = 300;             // ring radius
  var N = 24;               // number of test masses
  var DOT = 11;             // test-mass radius, px
  var EPS = 0.26;           // strain amplitude h
  var PERIOD = 4000;        // orbital period, ms

  var CY = 590;             // height of both centres, clear of the title
  var CXR = 960 - 450;      // ring centre
  var CXS = 960 + 450;      // segment centre
  var SEG = R0;             // segment half-length at |h| = EPS
  var WIDTH = 18;           // segment width, px

  // World coordinates have Y up, as in gw-ring; screen Y is down.
  function draw(ctx, phi, circular, showSegment) {
    var hp = EPS * Math.cos(2 * phi), hx = circular ? EPS * Math.sin(2 * phi) : 0;
    ctx.setTransform(2, 0, 0, 2, 0, 0);
    ctx.clearRect(0, 0, 1920, 1080);

    // Each mass moves by dx_i = (1/2) h_ij x_j.
    ctx.fillStyle = "#000";
    for (var k = 0; k < N; k++) {
      var th = 2 * Math.PI * k / N;
      var x = R0 * Math.cos(th), y = R0 * Math.sin(th);
      ctx.beginPath();
      ctx.arc(CXR + x + 0.5 * (hp * x + hx * y), CY - (y + 0.5 * (hx * x - hp * y)),
              DOT, 0, 2 * Math.PI);
      ctx.fill();
    }

    if (!showSegment) return;
    var psi = 0.5 * Math.atan2(hx, hp), half = SEG * Math.hypot(hp, hx) / EPS;
    var dx = half * Math.cos(psi), dy = half * Math.sin(psi);
    ctx.strokeStyle = "#000";
    ctx.lineWidth = WIDTH;
    ctx.lineCap = "butt";   // no caps, so that the bar has exactly this length and vanishes at h = 0
    ctx.beginPath();
    ctx.moveTo(CXS - dx, CY + dy);
    ctx.lineTo(CXS + dx, CY - dy);
    ctx.stroke();
  }

  function widget(id, circular) {
    var raf = null, ctx = null, start = 0, stage = 0;
    function loop(now) {
      draw(ctx, 2 * Math.PI * (now - start) / PERIOD, circular, stage >= 1);
      raf = requestAnimationFrame(loop);
    }
    function stop() { if (raf !== null) cancelAnimationFrame(raf); raf = null; }
    Deck.widget(id, {
      steps: 1,
      enter: function (slide) {
        ctx = slide.querySelector("canvas").getContext("2d");
        start = performance.now();
        stage = 0;
        stop();
        raf = requestAnimationFrame(loop);
      },
      leave: stop,
      // The clock runs on through the stages: only `stage` changes here.
      step: function (slide, k) { stage = k; }
    });
  }

  widget("polar-circular", true);
  widget("polar-plus", false);
})();
