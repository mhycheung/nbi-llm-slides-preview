// Shared deck engine: scaling, keyboard navigation, steps, widget hooks, math.
// Parts use only Deck.widget(); see slides/FORMAT.md.
(function () {
  "use strict";

  var widgets = {};   // slide id -> {steps, enter, leave, step}
  var slides = [];
  var totals = [];    // number of steps of each slide
  var cur = -1;       // current slide index
  var s = 0;          // current step on it
  var deckEl = null;

  var NEXT = { ArrowRight: 1, " ": 1, Spacebar: 1, PageDown: 1 };
  var PREV = { ArrowLeft: 1, PageUp: 1 };

  function fit() {
    var k = Math.min(window.innerWidth / 1920, window.innerHeight / 1080);
    var x = (window.innerWidth - 1920 * k) / 2;
    var y = (window.innerHeight - 1080 * k) / 2;
    deckEl.style.transform = "translate(" + x + "px," + y + "px) scale(" + k + ")";
  }

  function call(i, name, a, b) {
    var w = widgets[slides[i].id];
    if (w && typeof w[name] === "function") w[name](slides[i], a, b);
  }

  function showStep(i, step, dir) {
    var els = slides[i].querySelectorAll("[data-step]");
    for (var j = 0; j < els.length; j++) {
      els[j].classList.toggle("shown", Number(els[j].getAttribute("data-step")) <= step);
    }
    call(i, "step", step, dir);
  }

  function go(i, step, dir) {
    if (i < 0 || i >= slides.length) return;
    if (i !== cur) {
      if (cur >= 0) {
        call(cur, "leave");
        slides[cur].classList.remove("active");
      }
      cur = i;
      slides[cur].classList.add("active");
      call(cur, "enter");
      history.replaceState(null, "", "#" + (cur + 1));
    }
    s = step;
    showStep(cur, s, dir);
    var bar = deckEl.querySelector(".deck-progress-fill");
    if (bar) {
      var f = slides.length > 1 ? cur / (slides.length - 1) : 1;
      bar.style.width = (100 * f) + "%";
      deckEl.querySelector(".deck-progress-knob").style.left = (100 * f) + "%";
    }
  }

  // Progress bar (presentation only): pressing or dragging on it jumps to the slide under
  // the pointer, at its arrival stage. This is the one mouse control that changes slides.
  function seekable() {
    var el = deckEl.querySelector(".deck-progress");
    if (!el) return;
    var track = el.querySelector(".deck-progress-track");
    function seek(e) {
      var r = track.getBoundingClientRect();
      var f = Math.min(1, Math.max(0, (e.clientX - r.left) / r.width));
      var i = Math.round(f * (slides.length - 1));
      if (i !== cur) go(i, 0, 1);
    }
    el.addEventListener("pointerdown", function (e) {
      e.preventDefault();
      el.setPointerCapture(e.pointerId);
      el.classList.add("dragging");
      seek(e);
    });
    el.addEventListener("pointermove", function (e) {
      if (el.hasPointerCapture(e.pointerId)) seek(e);
    });
    function end(e) {
      el.classList.remove("dragging");
      if (el.hasPointerCapture(e.pointerId)) el.releasePointerCapture(e.pointerId);
    }
    el.addEventListener("pointerup", end);
    el.addEventListener("pointercancel", end);
  }

  // Time of day, bottom left; present only in the stitched deck.
  function tick() {
    var el = deckEl.querySelector(".deck-clock");
    if (!el) return;
    var d = new Date();
    el.textContent = d.getHours() + ":" + ("0" + d.getMinutes()).slice(-2);
  }

  function next() {
    if (s < totals[cur]) go(cur, s + 1, 1);
    else if (cur + 1 < slides.length) go(cur + 1, 0, 1);
  }

  function prev() {
    if (s > 0) go(cur, s - 1, -1);
    else if (cur > 0) go(cur - 1, totals[cur - 1], -1);
  }

  // Capture phase: runs before any widget sees the key, and cancels the key's default
  // action, so a focused slider or button never takes a navigation key.
  function onKey(e) {
    if (e.ctrlKey || e.altKey || e.metaKey) return;
    var n = NEXT[e.key], p = PREV[e.key];
    if (!n && !p) return;
    e.preventDefault();
    e.stopPropagation();
    if (e.type !== "keydown") return;
    if (n) next(); else prev();
  }

  // Every credit with one arXiv number, "Smith+ (2401.12345)", becomes a link to the
  // paper's arXiv abstract page. The credit's look does not change.
  function linkCites() {
    deckEl.querySelectorAll(".cite").forEach(function (c) {
      if (c.querySelector("a")) return;
      var ids = c.textContent.match(/\b\d{4}\.\d{4,5}\b/g);
      if (!ids || ids.length !== 1) return;
      var a = document.createElement("a");
      a.className = "cite-link";
      a.href = "https://arxiv.org/abs/" + ids[0];
      a.target = "_blank";
      while (c.firstChild) a.appendChild(c.firstChild);
      c.appendChild(a);
    });
  }

  function start() {
    deckEl = document.querySelector(".deck");
    slides = Array.prototype.slice.call(deckEl.querySelectorAll(".slide"));
    linkCites();
    if (window.renderMathInElement) {
      renderMathInElement(deckEl, {
        delimiters: [
          { left: "$$", right: "$$", display: true },
          { left: "\\[", right: "\\]", display: true },
          { left: "$", right: "$", display: false },
          { left: "\\(", right: "\\)", display: false }
        ],
        throwOnError: false
      });
    }
    totals = slides.map(function (sl) {
      var m = 0;
      sl.querySelectorAll("[data-step]").forEach(function (el) {
        m = Math.max(m, Number(el.getAttribute("data-step")));
      });
      var w = widgets[sl.id];
      return Math.max(m, (w && w.steps) || 0);
    });
    fit();
    tick();
    if (deckEl.querySelector(".deck-clock")) setInterval(tick, 5000);
    seekable();
    window.addEventListener("resize", fit);
    window.addEventListener("keydown", onKey, true);
    window.addEventListener("keyup", onKey, true);
    window.addEventListener("keypress", onKey, true);
    var h = parseInt(location.hash.slice(1), 10);
    go(h >= 1 && h <= slides.length ? h - 1 : 0, 0, 1);
  }

  window.Deck = {
    // Register the widget of slide `id`. hooks: {steps, enter(slide), leave(slide),
    // step(slide, step, dir)}; all optional. step 0 is the state on arrival from the
    // previous slide; dir is +1 going forward and -1 going back.
    widget: function (id, hooks) { widgets[id] = hooks; },
    start: start,
    go: function (n) { go(n - 1, 0, 1); }
  };
})();
