// Example part: shows a stepped animation and a mouse-driven widget.
(function () {
  "use strict";

  // Damped sinusoid on the canvas, up to fraction `upto` of its length.
  function draw(ctx, tau, freq, upto, color) {
    var w = ctx.canvas.width, h = ctx.canvas.height;
    ctx.strokeStyle = color;
    ctx.lineWidth = 4;
    ctx.beginPath();
    var n = Math.floor(w * upto);
    for (var x = 0; x <= n; x++) {
      var t = x / w;
      var y = h / 2 - 0.45 * h * Math.exp(-t / tau) * Math.cos(2 * Math.PI * freq * t);
      if (x === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }
    ctx.stroke();
  }

  function axis(ctx) {
    var w = ctx.canvas.width, h = ctx.canvas.height;
    ctx.clearRect(0, 0, w, h);
    ctx.strokeStyle = "#000";
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.moveTo(0, h / 2);
    ctx.lineTo(w, h / 2);
    ctx.stroke();
  }

  // Stage 0: axis. Stage 1: first mode drawn over 1 s. Stage 2: second mode added.
  var raf = null;
  function stop() { if (raf) cancelAnimationFrame(raf); raf = null; }
  Deck.widget("example-anim", {
    steps: 2,
    leave: stop,
    step: function (slide, step, dir) {
      stop();
      var ctx = slide.querySelector("canvas").getContext("2d");
      axis(ctx);
      if (step === 0) return;
      if (step === 2 || dir < 0) draw(ctx, 0.3, 8, 1, "#000");
      if (step === 1 && dir > 0) {
        var t0 = performance.now();
        (function frame(now) {
          var f = Math.min(1, (now - t0) / 1000);
          axis(ctx);
          draw(ctx, 0.3, 8, f, "#000");
          raf = f < 1 ? requestAnimationFrame(frame) : null;
        })(t0);
      }
      if (step === 2) draw(ctx, 0.15, 13, 1, "#888");
    }
  });

  // Slider sets the damping time; arrow keys on the focused slider still turn the page.
  Deck.widget("example-widget", {
    enter: function (slide) {
      var ctx = slide.querySelector("canvas").getContext("2d");
      var input = slide.querySelector("input");
      function redraw() { axis(ctx); draw(ctx, Number(input.value), 8, 1, "#000"); }
      input.oninput = redraw;
      redraw();
    }
  });
})();
