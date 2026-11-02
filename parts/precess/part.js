// Equal-mass binary on a circular Newtonian orbit of fixed radius, with spin precession.
// Slide precess-compare: stage 0 title only; stage 1 non-precessing (left) and precessing
// (right) side by side, same camera and same orbital rate.
// Slide precess-detectors: the precessing binary alone, smaller, with x, y, z axes, one Earth on
// the z axis (along J, face on) and one on the x axis (perpendicular to J, edge on). Dragging
// turns the camera about the binary, to any direction. A dotted circle marks the holes' current orbit.
(function () {
  "use strict";

  // ---- PHYSICS BEGIN (tests extract this block; no DOM below this line until PHYSICS END) ----
  // Orbit-averaged equal-mass precession, coupling exaggerated, no spin-spin term:
  //   dL/dt = k S x L,  dS_i/dt = k L x S_i,  S = S1 + S2,  J = L + S conserved.
  // L then precesses about J at rate k|J|. World axes: z along J.
  // Units of L, S: |L| = 1. Time in seconds.
  var PERIOD = 4;           // orbital period, s (as in gw-ring)
  var N_PREC = 5;           // orbital periods per precession period
  var THETA = 30 * Math.PI / 180;   // opening angle between L and J
  var SMAG = 1.5;           // |S|
  var SPLIT = 22 * Math.PI / 180;   // angle of each spin from S
  var YAW = 25 * Math.PI / 180;     // square edges against the camera's right axis
  // S on the far side of J from L, so that J = L + S lies at THETA from L: |S| sin(THETA_S) = sin(THETA).
  // SMAG > 1 and SPLIT keep every spin at most THETA_S + SPLIT = 41.5 deg from J, away from the
  // camera direction (52 deg from J), so that no arrow shrinks to a dot on screen.
  var THETA_S = Math.asin(Math.sin(THETA) / SMAG);
  var JMAG = Math.cos(THETA) + SMAG * Math.cos(THETA_S);
  var K = 2 * Math.PI / (N_PREC * PERIOD * JMAG);
  var DT = 1 / 120;         // fixed RK4 step, s
  var OMEGA = 2 * Math.PI / PERIOD;
  var SPIN = SMAG / (2 * Math.cos(SPLIT));   // |S1| = |S2|

  function cross(a, b) {
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  }
  function dot(a, b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
  function combo(a, f, b) { return [a[0] + f * b[0], a[1] + f * b[1], a[2] + f * b[2]]; }

  // State: 15 numbers, L, S1, S2, e1, e2 (three each).
  function vec(y, i) { return [y[i], y[i + 1], y[i + 2]]; }
  function deriv(y) {
    var L = vec(y, 0), S1 = vec(y, 3), S2 = vec(y, 6), e1 = vec(y, 9), e2 = vec(y, 12);
    var S = combo(S1, 1, S2);
    var dL = cross(S, L), n2 = dot(L, L);
    dL = [K * dL[0], K * dL[1], K * dL[2]];
    var dS1 = cross(L, S1), dS2 = cross(L, S2);
    // Transport without twist: de/dt = -Lhat (dLhat/dt . e), and |L| is constant here.
    var a1 = -dot(dL, e1) / n2, a2 = -dot(dL, e2) / n2;
    return [dL[0], dL[1], dL[2],
            K * dS1[0], K * dS1[1], K * dS1[2],
            K * dS2[0], K * dS2[1], K * dS2[2],
            a1 * L[0], a1 * L[1], a1 * L[2],
            a2 * L[0], a2 * L[1], a2 * L[2]];
  }
  function rk4(y, h) {
    var n = y.length, i, t = new Array(n);
    var k1 = deriv(y);
    for (i = 0; i < n; i++) t[i] = y[i] + 0.5 * h * k1[i];
    var k2 = deriv(t);
    for (i = 0; i < n; i++) t[i] = y[i] + 0.5 * h * k2[i];
    var k3 = deriv(t);
    for (i = 0; i < n; i++) t[i] = y[i] + h * k3[i];
    var k4 = deriv(t);
    var out = new Array(n);
    for (i = 0; i < n; i++) out[i] = y[i] + h * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
    return out;
  }
  // L tilted by THETA towards +x, S towards -x; the spins split along +-y.
  function initialState() {
    var L = [Math.sin(THETA), 0, Math.cos(THETA)], Sh = [-Math.sin(THETA_S), 0, Math.cos(THETA_S)];
    var d = SPIN * Math.sin(SPLIT), m = SPIN * Math.cos(SPLIT);
    var S1 = [m * Sh[0], d, m * Sh[2]], S2 = [m * Sh[0], -d, m * Sh[2]];
    // e1: the fixed square's edge direction, made perpendicular to L; e1 x e2 = L.
    var f = [Math.cos(YAW), Math.sin(YAW), 0];
    var e1 = combo(f, -dot(f, L), L), n = Math.sqrt(dot(e1, e1));
    e1 = [e1[0] / n, e1[1] / n, e1[2] / n];
    var e2 = cross(L, e1);
    return L.concat(S1, S2, e1, e2);
  }
  // ---- PHYSICS END ----

  // ---- drawing constants, slide px at scale 1 ----
  var ORB = 210;            // orbital radius
  var FA = 275;             // half-size of the square
  var BH = 40;              // black-hole radius
  var SPIN_PX = 180;        // arrow length per unit of spin, so |S_i| = SPIN gives about 146
  var COL = "#d00000";      // spin arrows
  var SHAFT = 7, HEAD = 26; // arrow shaft width and head length at scale 1
  var S_CMP = 0.84;         // scale of the two panels of the first slide
  var S_DET = 0.55;         // scale of the binary on the second slide
  var PANELS = [480, 1440], CY_CMP = 600;
  var CX_DET = 960, CY_DET = 590;
  var AXIS = 280, AXIS_W = 2, AXIS_HEAD = 18;   // axis half-length, width, arrow head, slide px
  var AXIS_COL = "#b8b8b8";
  var PLATE_FILL = "rgba(233,233,233,0.82)", PLATE_EDGE = "#c4c4c4";   // nearly opaque, to hide what is behind
  var LABEL_OFF = 34;       // axis label from the arrow tip, sideways, slide px
  var EARTH_D = 410, EARTH_R = 52;    // Earth's distance from the binary and radius, slide px
  var LON0 = 65;            // Earth's longitude offset, deg: Europe and Africa face the first view
  var OCEAN = [38, 104, 178], LAND = [86, 158, 76];
  var DOT_W = 3.5, DOT_GAP = 11;      // dotted orbit: dot size and spacing, slide px

  // ---- orthographic camera: right R, up U, towards the camera V (R x U = V) ----
  // As on the sphere slide: yaw about z, then pitch about the screen's right axis. Pitch 0 looks
  // down the z axis from above, -90 deg from the side with z up, -180 deg from below.
  var EL = 38 * Math.PI / 180;        // elevation above the plane perpendicular to J, first slide
  var EL_DET = 28 * Math.PI / 180, YAW_DET = -30 * Math.PI / 180;   // first view, second slide
  var DRAG = 0.008;                   // rad per slide px
  var R, U, V;
  function setCamera(yaw, pitch) {
    var cy = Math.cos(yaw), sy = Math.sin(yaw), cp = Math.cos(pitch), sp = Math.sin(pitch);
    R = [cy, -sy, 0];
    U = [sy * cp, cy * cp, -sp];
    V = [sy * sp, cy * sp, cp];
  }

  function proj(p, cx, cy, s) { return [cx + s * dot(p, R), cy - s * dot(p, U)]; }

  function square(ctx, e1, e2, cx, cy, s) {
    ctx.beginPath();
    [[-1, -1], [1, -1], [1, 1], [-1, 1]].forEach(function (q, i) {
      var p = proj([FA * (q[0] * e1[0] + q[1] * e2[0]), FA * (q[0] * e1[1] + q[1] * e2[1]),
                    FA * (q[0] * e1[2] + q[1] * e2[2])], cx, cy, s);
      if (i === 0) ctx.moveTo(p[0], p[1]); else ctx.lineTo(p[0], p[1]);
    });
    ctx.closePath();
    ctx.fillStyle = PLATE_FILL;
    ctx.fill();
    ctx.strokeStyle = PLATE_EDGE;
    ctx.lineWidth = 1.5;
    ctx.stroke();
  }

  // Arrow from the screen point a to b; the head shrinks with a foreshortened arrow.
  function arrow(ctx, a, b, s) {
    var dx = b[0] - a[0], dy = b[1] - a[1], n = Math.hypot(dx, dy);
    if (n < 1e-6) return;
    var ux = dx / n, uy = dy / n, h = Math.min(HEAD * s, 0.5 * n), w = 0.45 * h;
    var bx = b[0] - h * ux, by = b[1] - h * uy;
    ctx.strokeStyle = COL;
    ctx.fillStyle = COL;
    ctx.lineWidth = SHAFT * s;
    ctx.lineCap = "butt";
    ctx.beginPath(); ctx.moveTo(a[0], a[1]); ctx.lineTo(bx, by); ctx.stroke();
    ctx.beginPath();
    ctx.moveTo(b[0], b[1]);
    ctx.lineTo(bx - w * uy, by + w * ux);
    ctx.lineTo(bx + w * uy, by - w * ux);
    ctx.closePath();
    ctx.fill();
  }

  // The current orbit of the holes, dotted, in the plane of e1, e2.
  function orbit(ctx, e1, e2, cx, cy, s) {
    var n = 120;
    ctx.strokeStyle = "#000";
    ctx.lineWidth = DOT_W;
    ctx.lineCap = "round";
    ctx.setLineDash([0.01, DOT_GAP]);
    ctx.beginPath();
    for (var i = 0; i < n; i++) {
      var a = 2 * Math.PI * i / n, c = ORB * Math.cos(a), d = ORB * Math.sin(a);
      var p = proj([c * e1[0] + d * e2[0], c * e1[1] + d * e2[1], c * e1[2] + d * e2[2]], cx, cy, s);
      if (i === 0) ctx.moveTo(p[0], p[1]); else ctx.lineTo(p[0], p[1]);
    }
    ctx.closePath();
    ctx.stroke();
    ctx.setLineDash([]);
  }

  // Which side of the plate of normal L the camera is on: +1 or -1.
  function cameraSide(L) { return dot(L, V) >= 0 ? 1 : -1; }

  // One binary. An arrow on the far side of the plate from the camera goes behind the plate.
  // back(), front(), if given, draw what lies behind and in front of the plate.
  function binary(ctx, y, phi, cx, cy, s, back, front) {
    var L = vec(y, 0), e1 = vec(y, 9), e2 = vec(y, 12), side = cameraSide(L);
    var c = Math.cos(phi), d = Math.sin(phi);
    var pos = [c * e1[0] + d * e2[0], c * e1[1] + d * e2[1], c * e1[2] + d * e2[2]];
    var near = dot(pos, V) > 0 ? 1 : -1;     // +1: the hole at +pos is the nearer one
    function spinArrow(sg, behind) {
      var P0 = [sg * ORB * pos[0], sg * ORB * pos[1], sg * ORB * pos[2]], Si = vec(y, sg > 0 ? 3 : 6);
      if ((dot(Si, L) * side < 0) === behind) {
        arrow(ctx, proj(P0, cx, cy, s), proj(combo(P0, SPIN_PX, Si), cx, cy, s), s);
      }
    }
    if (back) back(side, L);
    [-near, near].forEach(function (sg) { spinArrow(sg, true); });
    square(ctx, e1, e2, cx, cy, s);
    orbit(ctx, e1, e2, cx, cy, s);
    [-near, near].forEach(function (sg) {
      var a = proj([sg * ORB * pos[0], sg * ORB * pos[1], sg * ORB * pos[2]], cx, cy, s);
      ctx.fillStyle = "#000";
      ctx.beginPath();
      ctx.arc(a[0], a[1], BH * s, 0, 2 * Math.PI);
      ctx.fill();
      spinArrow(sg, false);
    });
    if (front) front(side, L);
  }

  // Half of an axis, from the origin to sg times the unit vector d; arrow head on the + half.
  function axisHalf(ctx, d, sg, cx, cy) {
    var a = proj([0, 0, 0], cx, cy, 1), b = proj([sg * AXIS * d[0], sg * AXIS * d[1], sg * AXIS * d[2]], cx, cy, 1);
    var dx = b[0] - a[0], dy = b[1] - a[1], n = Math.hypot(dx, dy);
    if (n < 1e-6) return;
    var ux = dx / n, uy = dy / n, h = sg > 0 ? Math.min(AXIS_HEAD, 0.5 * n) : 0, w = 0.42 * h;
    var bx = b[0] - h * ux, by = b[1] - h * uy;
    ctx.strokeStyle = AXIS_COL;
    ctx.fillStyle = AXIS_COL;
    ctx.lineWidth = AXIS_W;
    ctx.lineCap = "butt";
    ctx.beginPath(); ctx.moveTo(a[0], a[1]); ctx.lineTo(bx, by); ctx.stroke();
    if (sg > 0) {
      ctx.beginPath();
      ctx.moveTo(b[0], b[1]);
      ctx.lineTo(bx - w * uy, by + w * ux);
      ctx.lineTo(bx + w * uy, by - w * ux);
      ctx.closePath();
      ctx.fill();
    }
  }

  // ---- Earth: a sphere with north along z, drawn pixel by pixel for the current camera. The
  // camera is orthographic, so both Earths look the same and share one image, redrawn on a drag. ----
  var maskBits = null, earthCanvas = null, earthKey = null;
  function loadMask() {
    var raw = atob(window.precess_earth_mask);
    maskBits = new Uint8Array(raw.length);
    for (var i = 0; i < raw.length; i++) maskBits[i] = raw.charCodeAt(i);
  }
  function isLand(lat, lon) {   // degrees
    var j = Math.min(179, Math.floor(90 - lat)), i = Math.floor(((lon + 180) % 360 + 360) % 360) % 360;
    var b = j * 360 + i;
    return (maskBits[b >> 3] >> (b & 7)) & 1;
  }
  function renderEarth() {
    if (!maskBits) loadMask();
    var n = EARTH_R * 4;       // the canvas has 2 device px per slide px
    if (!earthCanvas) {
      earthCanvas = document.createElement("canvas");
      earthCanvas.width = earthCanvas.height = n;
    }
    var c = earthCanvas.getContext("2d"), im = c.createImageData(n, n), d = im.data;
    var L = [-0.4, 0.5, 0.77], r = n / 2 - 1;
    for (var j = 0; j < n; j++) {
      for (var i = 0; i < n; i++) {
        var x = (i + 0.5 - n / 2) / r, y = (n / 2 - j - 0.5) / r, q = x * x + y * y, k = 4 * (j * n + i);
        var edge = (1 - Math.sqrt(q)) * r + 0.5;
        if (edge <= 0) continue;
        var z = q < 1 ? Math.sqrt(1 - q) : 0;
        var wx = x * R[0] + y * U[0] + z * V[0], wy = x * R[1] + y * U[1] + z * V[1], wz = x * R[2] + y * U[2] + z * V[2];
        var lat = Math.asin(Math.max(-1, Math.min(1, wz))) * 180 / Math.PI;
        var lon = Math.atan2(wy, wx) * 180 / Math.PI + LON0;
        var col = isLand(lat, lon) ? LAND : OCEAN;
        var sh = 0.62 + 0.45 * Math.max(0, x * L[0] + y * L[1] + z * L[2]);
        d[k] = Math.min(255, col[0] * sh);
        d[k + 1] = Math.min(255, col[1] * sh);
        d[k + 2] = Math.min(255, col[2] * sh);
        d[k + 3] = 255 * Math.min(1, edge);
      }
    }
    c.putImageData(im, 0, 0);
  }
  function earth(ctx, pos, cx, cy) {
    var p = proj(pos, cx, cy, 1);
    ctx.drawImage(earthCanvas, p[0] - EARTH_R, p[1] - EARTH_R, 2 * EARTH_R, 2 * EARTH_R);
  }

  // ---- fixed non-precessing state: everything along z ----
  var FIXED = [0, 0, 1, 0, 0, SPIN, 0, 0, SPIN,
               Math.cos(YAW), Math.sin(YAW), 0, -Math.sin(YAW), Math.cos(YAW), 0];

  var EARTHS = [[0, 0, EARTH_D], [EARTH_D, 0, 0]];
  var AXES = [{ d: [1, 0, 0], name: "x" }, { d: [0, 1, 0], name: "y" }, { d: [0, 0, 1], name: "z" }];

  // ---- one widget: simulation clock, canvas, one rAF loop ----
  // running(k): is the animation on at stage k. paint(ctx, y, phi): draw one frame.
  // onEnter(slide), if given, runs each time the slide is entered.
  function makeWidget(id, steps, startStage, paint, onEnter) {
    var ctx = null, raf = null, last = null, y = initialState(), n = 0, acc = 0, on = false;
    function clear() {
      ctx.setTransform(2, 0, 0, 2, 0, 0);
      ctx.clearRect(0, 0, 1920, 1080);
    }
    function frame() {
      clear();
      if (on) paint(ctx, y, OMEGA * n * DT);
    }
    function loop(now) {
      // A long gap (hidden tab) is not caught up: at most 0.25 s per frame.
      if (last !== null) acc += Math.min(0.25, (now - last) / 1000);
      last = now;
      while (acc >= DT) { y = rk4(y, DT); n++; acc -= DT; }
      frame();
      raf = requestAnimationFrame(loop);
    }
    function stop() { if (raf !== null) cancelAnimationFrame(raf); raf = null; last = null; }
    Deck.widget(id, {
      steps: steps,
      enter: function (slide) {
        ctx = slide.querySelector("canvas").getContext("2d");
        if (onEnter) onEnter(slide);
        stop();
      },
      leave: stop,
      step: function (slide, k, dir) {
        on = k >= startStage;
        // Arriving forwards starts from t = 0; arriving backwards keeps the state left behind.
        if (dir >= 0 && k === startStage) { y = initialState(); n = 0; acc = 0; last = null; }
        frame();
        if (on && raf === null) raf = requestAnimationFrame(loop);
        if (!on) stop();
      }
    });
  }

  makeWidget("precess-compare", 1, 1, function (ctx, y, phi) {
    setCamera(0, EL - Math.PI / 2);
    binary(ctx, FIXED, phi, PANELS[0], CY_CMP, S_CMP);
    binary(ctx, y, phi, PANELS[1], CY_CMP, S_CMP);
  });

  // ---- second slide: camera, drag, axis labels ----
  var camYaw = YAW_DET, camPitch = EL_DET - Math.PI / 2, wired = false, drag = null, labels = null;

  function wireDrag(canvas) {
    function end() { drag = null; canvas.classList.remove("precess-grabbing"); }
    canvas.addEventListener("pointerdown", function (e) {
      canvas.setPointerCapture(e.pointerId);
      drag = { x: e.clientX, y: e.clientY, k: canvas.getBoundingClientRect().width / 1920 };
      canvas.classList.add("precess-grabbing");
    });
    canvas.addEventListener("pointermove", function (e) {
      if (!drag) return;
      camYaw += (e.clientX - drag.x) / drag.k * DRAG;
      camPitch = Math.max(-Math.PI, Math.min(0, camPitch + (e.clientY - drag.y) / drag.k * DRAG));
      drag.x = e.clientX; drag.y = e.clientY;
    });
    canvas.addEventListener("pointerup", end);
    canvas.addEventListener("pointercancel", end);
  }

  // Beside the arrow tip, on the side away from the axis's own line: to the right of a
  // vertical axis, above a horizontal one.
  function placeLabel(el, d, cx, cy) {
    var o = proj([0, 0, 0], cx, cy, 1), t = proj([AXIS * d[0], AXIS * d[1], AXIS * d[2]], cx, cy, 1);
    var dx = t[0] - o[0], dy = t[1] - o[1], n = Math.hypot(dx, dy) || 1;
    var px = dy / n, py = -dx / n;
    if (py > 0.2 || (py > -0.2 && px < 0)) { px = -px; py = -py; }
    el.style.left = (t[0] + LABEL_OFF * px) + "px";
    el.style.top = (t[1] + LABEL_OFF * py) + "px";
  }

  makeWidget("precess-detectors", 0, 0, function (ctx, y, phi) {
    setCamera(camYaw, camPitch);
    var key = camYaw + "," + camPitch;
    if (key !== earthKey) { renderEarth(); earthKey = key; }
    var cx = CX_DET, cy = CY_DET;
    // Each half axis and each Earth is on one side of the plate; each side is drawn far to near.
    function items(onCameraSide, side, L) {
      var list = [];
      AXES.forEach(function (a) {
        [1, -1].forEach(function (sg) {
          if ((sg * dot(a.d, L) * side >= 0) !== onCameraSide) return;
          list.push({ z: 0.5 * AXIS * sg * dot(a.d, V), draw: function () { axisHalf(ctx, a.d, sg, cx, cy); } });
        });
      });
      EARTHS.forEach(function (c) {
        if ((dot(c, L) * side >= 0) !== onCameraSide) return;
        list.push({ z: dot(c, V), draw: function () { earth(ctx, c, cx, cy); } });
      });
      list.sort(function (a, b) { return a.z - b.z; }).forEach(function (o) { o.draw(); });
    }
    binary(ctx, y, phi, cx, cy, S_DET,
      function (side, L) { items(false, side, L); },
      function (side, L) { items(true, side, L); });
    AXES.forEach(function (a) { placeLabel(labels[a.name], a.d, cx, cy); });
  }, function (slide) {
    var canvas = slide.querySelector("canvas");
    if (!wired) {
      wireDrag(canvas);
      labels = {};
      AXES.forEach(function (a) { labels[a.name] = slide.querySelector(".precess-axis-" + a.name); });
      wired = true;
    }
    camYaw = YAW_DET;
    camPitch = EL_DET - Math.PI / 2;
  });
})();
