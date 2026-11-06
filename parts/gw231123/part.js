// GW231123, slide gw231123-sphere: the NRSur7dq4 maximum-likelihood waveform on a sphere of viewing
// directions, zoomed in and seen straight on along the actual viewing direction, beside the whitened
// data and waveform of both detectors. Data: window.gw231123_data (assets/data.js).
// Each tile has h_+, h_x, in the source's L frame at f_ref = 10 Hz. A bar is along the stretch axis
// cos(psi) e_theta + sin(psi) e_phi, psi = (1/2) atan2(h_x, h_+), length proportional to |h|.
// The surface is one gray level per tile, from log |h|, filled over the tile's Voronoi cell. A click puts a red circle around a direction's bar, in the tangent plane;
// the same direction clears it, another moves it. The actual direction starts highlighted.
(function () {
  "use strict";

  var CX = 640, CY = 540, RAD = 620;   // sphere centre and radius, slide px: zoomed, so it passes the slide's edges
  var BAR = 90;                  // bar length at |h| = 1, slide px
  var FPS = 24, HOLD = 24;       // data frames per second, about 1/41 of real time; frames held at the end
  var H_MIN = 0.03, SH_MAX = 0.36;   // gray level 1 - SH_MAX * s, s = 1 + log|h| / log(1 / H_MIN), clipped to 0..1
  var GROW = 1.04;               // cells drawn this much larger about their tile, to cover the seams between them
  var BASE = "#ebebeb";          // disc under the cells: what shows in the slivers at the sphere's edge
  var COL = "#d00000";           // highlight and time cursor
  var HL_R = 48, HL_W = 4.5;     // highlight circle radius and width, slide px: around the longest bar
  var PICK_MOVE = 6, PICK_R = 40;   // a click moves less than this; reach of a click, slide px
  var PX0 = 1300, PX1 = 1880;    // whitened traces: horizontal extent, and each row's top and bottom
  var ROWS = [{ ifo: "H1", top: 170, bot: 500 }, { ifo: "L1", top: 610, bot: 940 }];
  var BAND = "rgba(100,197,235,0.5)", DATA = "#8a8a8a";
  var TRACE_X = 2 * (PX0 - 20);  // left edge of the part of the traces image that is copied each frame, canvas px

  var D = null, n = 0, frames = null;
  var Tx, Ty, Tz, Px, Py, Nx, Ny, Nz, X, Y, Z, TX, TY, PX, PY, H, order;
  var VN = 0, VS = null, VZ = null, table = null, grays = null, trace = null;
  var yaw0 = 0, pitch0 = 0, yaw = 0, pitch = 0, cyw = 1, syw = 0, cpi = 1, spi = 0, dirty = true;
  var hl = {};

  function load() {
    D = window.gw231123_data;
    n = D.nTiles;
    var raw = atob(D.frames.replace(/-/g, '+').replace(/_/g, '/')), bytes = new Uint8Array(raw.length);
    for (var b = 0; b < raw.length; b++) bytes[b] = raw.charCodeAt(b);
    frames = new Int16Array(bytes.buffer);
    function F() { return new Float32Array(n); }
    Tx = F(); Ty = F(); Tz = F(); Px = F(); Py = F(); Nx = F(); Ny = F(); Nz = F();
    X = F(); Y = F(); Z = F(); TX = F(); TY = F(); PX = F(); PY = F(); H = F();
    order = new Uint16Array(n);
    for (var i = 0; i < n; i++) {
      var ct = Math.cos(D.theta[i]), st = Math.sin(D.theta[i]), cp = Math.cos(D.phi[i]), sp = Math.sin(D.phi[i]);
      Nx[i] = st * cp; Ny[i] = st * sp; Nz[i] = ct;
      Tx[i] = ct * cp; Ty[i] = ct * sp; Tz[i] = -st;
      Px[i] = -sp; Py[i] = cp;
      order[i] = i;
    }
    // Camera straight on along the actual viewing direction: towards-camera vector = the tile's normal.
    var m = D.maxL;
    pitch0 = -D.theta[m];
    yaw0 = Math.atan2(-Math.cos(D.phi[m]), -Math.sin(D.phi[m]));
    buildTable();
    buildTrace();
  }

  // ---- gray level of |h| ----
  function buildTable() {
    table = new Uint8ClampedArray(1025);
    for (var q = 0; q <= 1024; q++) {
      var s = q === 0 ? 0 : 1 + Math.log(q / 1024) / Math.log(1 / H_MIN);
      table[q] = 255 * (1 - SH_MAX * Math.max(0, Math.min(1, s)));
    }
    grays = [];
    for (var g = 0; g < 256; g++) grays.push("rgb(" + g + "," + g + "," + g + ")");
    VN = D.verts.length / 3;
    VS = new Float32Array(2 * VN); VZ = new Float32Array(VN);
  }

  // ---- the whitened traces, drawn once ----
  function buildTrace() {
    trace = document.createElement("canvas");
    trace.width = 3840; trace.height = 2160;
    var c = trace.getContext("2d");
    c.setTransform(2, 0, 0, 2, 0, 0);
    var t = D.tTrace, t0 = t[0], dt = t[t.length - 1] - t0;
    function xOf(v) { return PX0 + (PX1 - PX0) * (v - t0) / dt; }
    ROWS.forEach(function (r) {
      var T = D.traces[r.ifo], m = 0;
      T.data.forEach(function (v) { m = Math.max(m, Math.abs(v)); });
      function yOf(v) { return 0.5 * (r.top + r.bot) - 0.5 * (r.bot - r.top) * v / m; }
      r.yOf = yOf;
      function line(a, col, w) {
        c.strokeStyle = col; c.lineWidth = w;
        c.beginPath();
        a.forEach(function (v, i) { if (i === 0) c.moveTo(xOf(t[i]), yOf(v)); else c.lineTo(xOf(t[i]), yOf(v)); });
        c.stroke();
      }
      line(T.data, DATA, 1.6);
      c.beginPath();
      T.upper.forEach(function (v, i) { if (i === 0) c.moveTo(xOf(t[i]), yOf(v)); else c.lineTo(xOf(t[i]), yOf(v)); });
      for (var i = t.length - 1; i >= 0; i--) c.lineTo(xOf(t[i]), yOf(T.lower[i]));
      c.closePath();
      c.fillStyle = BAND;
      c.fill();
      line(T.ml, "#000", 2.5);
    });
    trace.xOf = xOf;
  }

  // ---- camera: v = Rx(pitch) Rz(yaw) p; depth = v_z ----
  function byDepth(a, b) { return Z[a] - Z[b]; }
  function reproject() {
    cyw = Math.cos(yaw); syw = Math.sin(yaw); cpi = Math.cos(pitch); spi = Math.sin(pitch);
    for (var j = 0; j < n; j++) {
      var y1 = Nx[j] * syw + Ny[j] * cyw;
      X[j] = CX + RAD * (Nx[j] * cyw - Ny[j] * syw);
      Y[j] = CY - RAD * (y1 * cpi - Nz[j] * spi);
      Z[j] = y1 * spi + Nz[j] * cpi;
      y1 = Tx[j] * syw + Ty[j] * cyw;
      TX[j] = Tx[j] * cyw - Ty[j] * syw;
      TY[j] = -(y1 * cpi - Tz[j] * spi);
      y1 = Px[j] * syw + Py[j] * cyw;
      PX[j] = Px[j] * cyw - Py[j] * syw;
      PY[j] = -y1 * cpi;
    }
    order.sort(byDepth);
    // Cell vertices on the screen.
    for (var v = 0; v < VN; v++) {
      var vx = D.verts[3 * v], vy = D.verts[3 * v + 1], vz = D.verts[3 * v + 2];
      var y2 = vx * syw + vy * cyw, sx = vx * cyw - vy * syw, sy = -(y2 * cpi - vz * spi), z2 = y2 * spi + vz * cpi;
      if (z2 < 0) { var f = 1 / Math.hypot(sx, sy); sx *= f; sy *= f; }   // behind the sphere: on its edge
      VS[2 * v] = CX + RAD * sx; VS[2 * v + 1] = CY + RAD * sy; VZ[v] = z2;
    }
    dirty = false;
  }

  // ---- frame ----
  function loadFrame(f) {
    var off = f * n * 2;
    for (var i = 0; i < n; i++) {
      var a = frames[off + 2 * i], b = frames[off + 2 * i + 1];
      H[i] = Math.sqrt(a * a + b * b) / 32767;
    }
  }
  function drawShade(ctx) {
    // No clip, and no stroke round each cell: either one costs far more than the fills.
    ctx.fillStyle = BASE;
    ctx.beginPath();
    ctx.arc(CX, CY, RAD, 0, 2 * Math.PI);
    ctx.fill();
    for (var i = 0; i < n; i++) {
      var j = order[i], cell = D.cells[j], m = cell.length, any = false;
      for (var k = 0; k < m; k++) if (VZ[cell[k]] >= 0) { any = true; break; }
      if (!any) continue;
      ctx.fillStyle = grays[table[Math.min(1024, Math.round(H[j] * 1024))]];
      ctx.beginPath();
      for (k = 0; k < m; k++) {
        var x = X[j] + GROW * (VS[2 * cell[k]] - X[j]), y = Y[j] + GROW * (VS[2 * cell[k] + 1] - Y[j]);
        if (k === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
      }
      ctx.closePath();
      ctx.fill();
    }
    ctx.globalAlpha = 0.35;
    ctx.strokeStyle = "#000";
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.arc(CX, CY, RAD, 0, 2 * Math.PI);
    ctx.stroke();
    ctx.globalAlpha = 1;
  }
  function width(u, z) { return 2 * (1.2 + 2.3 * u) * (0.7 + 0.3 * z); }
  function drawBars(ctx, f) {
    var off = f * n * 2, half = 0.5 * BAR;
    ctx.strokeStyle = "#000";
    ctx.lineCap = "round";
    function bar(j) {
      var psi = 0.5 * Math.atan2(frames[off + 2 * j + 1], frames[off + 2 * j]), c = Math.cos(psi), s = Math.sin(psi);
      var u = H[j], dx = half * u * (c * TX[j] + s * PX[j]), dy = half * u * (c * TY[j] + s * PY[j]);
      ctx.lineWidth = width(u, Z[j]);
      ctx.beginPath();
      ctx.moveTo(X[j] - dx, Y[j] - dy);
      ctx.lineTo(X[j] + dx, Y[j] + dy);
      ctx.stroke();
    }
    for (var i = 0; i < n; i++) {
      var j = order[i];
      if (Z[j] < 0) continue;
      ctx.globalAlpha = 0.55 + 0.45 * Z[j];
      bar(j);
    }
    ctx.globalAlpha = 1;
    ctx.strokeStyle = COL;
    ctx.lineWidth = HL_W;
    Object.keys(hl).forEach(function (key) {
      var j = +key;
      if (Z[j] < 0) return;
      ctx.beginPath();
      for (var i = 0; i < 48; i++) {
        var a = 2 * Math.PI * i / 48, c = Math.cos(a), d = Math.sin(a);
        var x = X[j] + HL_R * (c * TX[j] + d * PX[j]), y = Y[j] + HL_R * (c * TY[j] + d * PY[j]);
        if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
      }
      ctx.closePath();
      ctx.stroke();
    });
  }

  var ctx = null, canvas = null, raf = null, t0 = 0, drag = null, frameNow = 0;
  function draw(now) {
    var f = Math.floor((now - t0) / 1000 * FPS) % (D.nFrames + HOLD);
    f = Math.min(f, D.nFrames - 1);
    frameNow = f;
    if (dirty) reproject();
    ctx.setTransform(2, 0, 0, 2, 0, 0);
    ctx.clearRect(0, 0, 1920, 1080);
    loadFrame(f);
    drawShade(ctx);
    drawBars(ctx, f);
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.drawImage(trace, TRACE_X, 0, 3840 - TRACE_X, 2160, TRACE_X, 0, 3840 - TRACE_X, 2160);
    ctx.setTransform(2, 0, 0, 2, 0, 0);
    var x = trace.xOf(D.tFrames[f]);
    ctx.strokeStyle = COL;
    ctx.lineWidth = 2.5;
    ctx.lineCap = "butt";
    ROWS.forEach(function (r) {
      ctx.beginPath(); ctx.moveTo(x, r.top); ctx.lineTo(x, r.bot); ctx.stroke();
    });
  }
  function loop(now) { draw(now); raf = requestAnimationFrame(loop); }
  function stop() { if (raf !== null) cancelAnimationFrame(raf); raf = null; }
  function endDrag() { drag = null; if (canvas) canvas.classList.remove("gw231123-grabbing"); }

  function pick(px, py) {
    var best = -1, bestD = PICK_R;
    for (var j = 0; j < n; j++) {
      if (Z[j] < 0) continue;
      var d = Math.hypot(X[j] - px, Y[j] - py);
      if (d < bestD) { bestD = d; best = j; }
    }
    if (best < 0) return;
    var had = hl[best];
    hl = {};
    if (!had) hl[best] = true;
  }

  function wire() {
    canvas.addEventListener("pointerdown", function (e) {
      canvas.setPointerCapture(e.pointerId);
      drag = { x: e.clientX, y: e.clientY, k: canvas.getBoundingClientRect().width / 1920, moved: 0 };
      canvas.classList.add("gw231123-grabbing");
    });
    canvas.addEventListener("pointermove", function (e) {
      if (!drag) return;
      var dx = (e.clientX - drag.x) / drag.k, dy = (e.clientY - drag.y) / drag.k;
      drag.x = e.clientX; drag.y = e.clientY;
      drag.moved += Math.hypot(dx, dy);
      if (drag.moved < PICK_MOVE) return;
      // Slower than on the full sphere: the sphere is larger on the slide.
      yaw += dx * 0.004;
      pitch = Math.max(-Math.PI, Math.min(0, pitch + dy * 0.004));
      dirty = true;
    });
    canvas.addEventListener("pointerup", function (e) {
      if (drag && drag.moved < PICK_MOVE) {
        var r = canvas.getBoundingClientRect();
        if (dirty) reproject();
        pick((e.clientX - r.left) / drag.k, (e.clientY - r.top) / drag.k);
      }
      endDrag();
    });
    canvas.addEventListener("pointercancel", endDrag);
  }

  Deck.widget("gw231123-sphere", {
    steps: 0,
    enter: function (slide) {
      if (!D) {
        load();
        canvas = slide.querySelector("canvas");
        ctx = canvas.getContext("2d");
        wire();
      }
      yaw = yaw0; pitch = pitch0; dirty = true;
      hl = {}; hl[D.maxL] = true;
      t0 = performance.now();
      stop();
      raf = requestAnimationFrame(loop);
    },
    leave: function () { stop(); endDrag(); },
    // Not used by the engine: state for the checks.
    _check: {
      frame: function () { return frameNow; },
      highlighted: function () { return Object.keys(hl).map(Number); },
      screenOf: function (j) { if (dirty) reproject(); return [X[j], Y[j], Z[j]]; },
      seek: function (now, f) { t0 = now - (f + 0.5) / FPS * 1000; }
    }
  });
})();
