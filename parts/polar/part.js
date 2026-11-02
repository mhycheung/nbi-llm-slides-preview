// Ring of test masses beside the line segment that stands for the same wave.
// Stage 0: the ring alone. Stage 1: the segment appears, on the same clock.
// A slider below sets the strain amplitude of both; at its left end the ring is a circle and the
// segment has zero length.
// The segment is the "stretch whisker" of the workshop visualizer: a headless bar centred
// on the point, along psi = 0.5 atan2(h_x, h_+), with length proportional to hypot(h_+, h_x).
(function () {
  "use strict";

  // Constants of the gw-ring part, which this part continues.
  var R0 = 300;             // ring radius
  var N = 24;               // number of test masses
  var DOT = 11;             // test-mass radius, px
  var EPS = 0.26;           // strain amplitude h on arrival
  var EPS_MAX = 0.36;       // right end of the slider; the ring and the segment stay clear of the slider and title
  var PERIOD = 4000;        // orbital period, ms

  var CY = 590;             // height of both centres, clear of the title
  var CXR = 960 - 450;      // ring centre
  var CXS = 960 + 450;      // segment centre
  var SEG = R0;             // segment half-length at |h| = EPS; the length is proportional to |h|
  var WIDTH = 18;           // segment width, px

  // World coordinates have Y up, as in gw-ring; screen Y is down.
  function draw(ctx, phi, circular, showSegment, amp) {
    var hp = amp * Math.cos(2 * phi), hx = circular ? amp * Math.sin(2 * phi) : 0;
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

    if (!showSegment || amp === 0) return;
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
    var raf = null, ctx = null, start = 0, stage = 0, amp = EPS, slider = null, knob = null;
    function loop(now) {
      draw(ctx, 2 * Math.PI * (now - start) / PERIOD, circular, stage >= 1, amp);
      raf = requestAnimationFrame(loop);
    }
    function stop() { if (raf !== null) cancelAnimationFrame(raf); raf = null; }
    function setAmp(a) {
      amp = a;
      knob.style.left = 100 * amp / EPS_MAX + "%";
    }
    // Slider: left end is exactly zero amplitude. Pointer position is read against the track's
    // rectangle on screen, so it holds at any window size.
    function wire() {
      function seek(e) {
        var r = slider.getBoundingClientRect();
        setAmp(EPS_MAX * Math.max(0, Math.min(1, (e.clientX - r.left) / r.width)));
      }
      slider.addEventListener("pointerdown", function (e) {
        e.preventDefault();
        slider.setPointerCapture(e.pointerId);
        seek(e);
      });
      slider.addEventListener("pointermove", function (e) {
        if (slider.hasPointerCapture(e.pointerId)) seek(e);
      });
    }
    Deck.widget(id, {
      steps: 1,
      enter: function (slide) {
        ctx = slide.querySelector("canvas").getContext("2d");
        if (!slider) {
          slider = slide.querySelector(".polar-slider");
          knob = slide.querySelector(".polar-knob");
          wire();
        }
        setAmp(EPS);
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
