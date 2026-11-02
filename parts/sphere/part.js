// Non-precessing equal-mass binary and its quadrupole strain, in phase with the orbit.
// Stage 0: the binary and four rings of test masses (+z, -z, +x, -x), drawn facing the viewer.
// Stage 1: each ring becomes a bar along the stretch axis, length proportional to |h|.
// Stage 2: the bars move onto the sphere; 600 segments in 3D projection, draggable.
//
// Field, the one formula for everything drawn. With Q = (1/2)[cos 2Phi (xx - yy)
// + sin 2Phi (xy + yx)] projected on the tangent vectors e_theta, e_phi and normalised by 2:
//   h_+ = (Q_tt - Q_pp) = (1 + cos^2 theta)/2 cos 2(Phi - phi)
//   h_x = 2 Q_tp        = cos theta sin 2(Phi - phi)
// so |h| = 1 at the poles and at most 1/2 on the equator. The stretch axis is
// cos(psi) e_theta + sin(psi) e_phi with psi = (1/2) atan2(h_x, h_+).
(function () {
  "use strict";

  var PERIOD = 4000;        // orbital period, ms
  var CX = 960, CY = 540;   // slide centre, slide px
  var DIST = 360;           // ring-centre distance and sphere radius, slide px
  var ORB = 120;            // orbital radius of each hole, slide px
  var BH = 30;              // black-hole radius, slide px
  var R0 = 100;             // ring radius, slide px
  var NM = 24;              // test masses per ring
  var DOT = 6.5;            // test-mass radius, slide px
  var DISP = 0.45;          // display strain multiplying h on the rings (|h| = 1 at the poles)
  var BAR_HALF = 100;       // half-length of a stage-1 bar at |h| = 1, slide px
  var BAR_W = 11;           // stage-1 bar width, slide px
  var NG = 600;             // Fibonacci grid points (as in the workshop sphere)
  var WLEN = 0.14;          // segment length at |h| = 1, in units of the sphere radius
  var WSCALE = 2;           // segment width relative to the workshop formula, for a projector
  var ELEV = 22 * Math.PI / 180;       // camera elevation above the orbital plane
  var YAW0 = 0, PITCH0 = ELEV - Math.PI / 2;
  var XF_MS = 450, MORPH_MS = 1000;    // crossfade and move-to-sphere durations

  // ---- animated parameters, as in gw-ring ----
  function ease(f) { return f < 0.5 ? 2 * f * f : 1 - 2 * (1 - f) * (1 - f); }
  function P(v) { return { from: v, to: v, t0: 0, dur: 1 }; }
  function val(p, now) {
    var f = Math.max(0, Math.min(1, (now - p.t0) / p.dur));
    return p.from + (p.to - p.from) * ease(f);
  }
  function set(p, target, now, instant, dur) {
    if (instant) { p.from = p.to = target; return; }
    if (p.to === target) return;
    p.from = val(p, now); p.to = target; p.t0 = now; p.dur = dur;
  }
  var xf = P(0), morph = P(0);   // 0 rings -> 1 bars; 0 flat -> 1 sphere

  // ---- sky points: the grid, then the four ring directions at indices NG..NG+3 ----
  var NP = NG + 4;
  function F() { return new Float64Array(NP); }
  var Nx = F(), Ny = F(), Nz = F(), Tx3 = F(), Ty3 = F(), Tz3 = F(), Px3 = F(), Py3 = F();
  var FA = F(), FB = F(), C2 = F(), S2 = F();   // (1+ct^2)/2, ct, cos 2phi, sin 2phi
  function addPoint(j, ct, st, cp, sp) {
    Nx[j] = st * cp; Ny[j] = st * sp; Nz[j] = ct;
    Tx3[j] = ct * cp; Ty3[j] = ct * sp; Tz3[j] = -st;
    Px3[j] = -sp; Py3[j] = cp;
    FA[j] = (1 + ct * ct) / 2; FB[j] = ct;
    C2[j] = cp * cp - sp * sp; S2[j] = 2 * sp * cp;
  }
  var GOLDEN = (1 + Math.sqrt(5)) / 2;
  for (var g = 0; g < NG; g++) {
    var gz = 1 - (2 * g + 1) / NG, gphi = 2 * Math.PI * g / GOLDEN;
    addPoint(g, gz, Math.sqrt(1 - gz * gz), Math.cos(gphi), Math.sin(gphi));
  }
  // Ring directions. Exact values, so that h_x = 0 on the equator and the poles are exact.
  // sx, sy: position of the ring on the screen (y up). rx, ux: world vectors of the ring's
  // screen right and up. Top and bottom are both seen from above (x right, y up), so that
  // both turn with the holes. On the poles e_theta = up; on the equator e_theta = -up.
  var RINGS = [
    { j: NG,     ct:  1, st: 0, cp:  0, sp:  1, sx: 0, sy:  DIST, rx: [1, 0, 0],  ux: [0, 1, 0] },
    { j: NG + 1, ct: -1, st: 0, cp:  0, sp: -1, sx: 0, sy: -DIST, rx: [1, 0, 0],  ux: [0, 1, 0] },
    { j: NG + 2, ct:  0, st: 1, cp:  1, sp:  0, sx:  DIST, sy: 0, rx: [0, 1, 0],  ux: [0, 0, 1] },
    { j: NG + 3, ct:  0, st: 1, cp: -1, sp:  0, sx: -DIST, sy: 0, rx: [0, -1, 0], ux: [0, 0, 1] }
  ];
  RINGS.forEach(function (r) { addPoint(r.j, r.ct, r.st, r.cp, r.sp); });
  var MC = new Float64Array(NM), MS = new Float64Array(NM);
  for (var m0 = 0; m0 < NM; m0++) {
    MC[m0] = R0 * Math.cos(2 * Math.PI * m0 / NM);
    MS[m0] = R0 * Math.sin(2 * Math.PI * m0 / NM);
  }

  // ---- field at point j for the current orbital phase ----
  var c2Phi = 1, s2Phi = 0;      // cos 2Phi, sin 2Phi of the frame
  var hMag = 0, cPsi = 1, sPsi = 0;
  function evalField(j) {
    var hp = FA[j] * (c2Phi * C2[j] + s2Phi * S2[j]);
    var hx = FB[j] * (s2Phi * C2[j] - c2Phi * S2[j]);
    var psi = 0.5 * Math.atan2(hx, hp);
    hMag = Math.sqrt(hp * hp + hx * hx);
    cPsi = Math.cos(psi); sPsi = Math.sin(psi);
  }

  // ---- camera, as in the workshop sphere: v = Rx(pitch) Rz(yaw) p; depth = v_z ----
  var yaw = YAW0, pitch = PITCH0, cyw = 1, syw = 0, cpi = 1, spi = 0, projDirty = true;
  var X = F(), Y = F(), Z = F(), TX = F(), TY = F(), PX = F(), PY = F();
  var order = new Uint16Array(NG);
  for (var o0 = 0; o0 < NG; o0++) order[o0] = o0;
  function byDepth(a, b) { return Z[a] - Z[b]; }
  function reproject() {
    cyw = Math.cos(yaw); syw = Math.sin(yaw); cpi = Math.cos(pitch); spi = Math.sin(pitch);
    for (var j = 0; j < NP; j++) {
      var y1 = Nx[j] * syw + Ny[j] * cyw;
      X[j] = CX + DIST * (Nx[j] * cyw - Ny[j] * syw);
      Y[j] = CY - DIST * (y1 * cpi - Nz[j] * spi);
      Z[j] = y1 * spi + Nz[j] * cpi;
      y1 = Tx3[j] * syw + Ty3[j] * cyw;
      TX[j] = DIST * (Tx3[j] * cyw - Ty3[j] * syw);
      TY[j] = -DIST * (y1 * cpi - Tz3[j] * spi);
      y1 = Px3[j] * syw + Py3[j] * cyw;
      PX[j] = DIST * (Px3[j] * cyw - Py3[j] * syw);
      PY[j] = -DIST * y1 * cpi;
    }
    order.sort(byDepth);
    projDirty = false;
  }

  function whiskerWidth(u, z) { return WSCALE * (1.2 + 2.3 * u) * (0.7 + 0.3 * z); }

  // ---- the parts of a frame ----
  function drawHoles(ctx, Phi) {
    var cp = Math.cos(Phi), sp = Math.sin(Phi);
    var vx = ORB * (cp * cyw - sp * syw), vy = ORB * (cp * syw + sp * cyw) * cpi;
    ctx.fillStyle = "#000";
    ctx.globalAlpha = 1;
    ctx.beginPath();
    ctx.arc(CX + vx, CY - vy, BH, 0, 2 * Math.PI);
    ctx.moveTo(CX - vx + BH, CY + vy);
    ctx.arc(CX - vx, CY + vy, BH, 0, 2 * Math.PI);
    ctx.fill();
  }

  // Stretch axis of ring r in its screen plane (right, up): unit vector stored in axA, axB.
  var axA = 1, axB = 0;
  function ringAxis(r) {
    var j = r.j;
    var ax = cPsi * Tx3[j] + sPsi * Px3[j], ay = cPsi * Ty3[j] + sPsi * Py3[j], az = cPsi * Tz3[j];
    axA = ax * r.rx[0] + ay * r.rx[1] + az * r.rx[2];
    axB = ax * r.ux[0] + ay * r.ux[1] + az * r.ux[2];
  }

  // Each mass moves by dx_i = (1/2) h_ij x_j, h_ij = DISP |h| (2 a_i a_j - delta_ij), a = axis.
  function drawRings(ctx, alpha) {
    ctx.globalAlpha = alpha;
    ctx.fillStyle = "#000";
    ctx.beginPath();
    for (var i = 0; i < 4; i++) {
      var r = RINGS[i];
      evalField(r.j);
      ringAxis(r);
      var hd = DISP * hMag, hxx = hd * (axA * axA - axB * axB), hxy = 2 * hd * axA * axB;
      for (var k = 0; k < NM; k++) {
        var x = MC[k], y = MS[k];
        var sx = CX + r.sx + x + 0.5 * (hxx * x + hxy * y);
        var sy = CY - r.sy - y - 0.5 * (hxy * x - hxx * y);
        ctx.moveTo(sx + DOT, sy);
        ctx.arc(sx, sy, DOT, 0, 2 * Math.PI);
      }
    }
    ctx.fill();
    ctx.globalAlpha = 1;
  }

  // Bars at their flat places (m = 0), moving to their true places on the sphere (m = 1).
  function drawBars(ctx, alpha, m) {
    var fade = 1 - Math.max(0, Math.min(1, (m - 0.6) / 0.4));
    if (alpha * fade <= 0) return;
    ctx.globalAlpha = alpha * fade;
    ctx.strokeStyle = "#000";
    ctx.lineCap = "butt";
    var half = 0.5 * WLEN;
    for (var i = 0; i < 4; i++) {
      var r = RINGS[i], j = r.j;
      evalField(j);
      ringAxis(r);
      var ax = BAR_HALF * hMag * axA, ay = -BAR_HALF * hMag * axB;
      var bx = half * hMag * (cPsi * TX[j] + sPsi * PX[j]);
      var by = half * hMag * (cPsi * TY[j] + sPsi * PY[j]);
      if (ax * bx + ay * by < 0) { bx = -bx; by = -by; }   // the axis has no sign
      var hx = ax + (bx - ax) * m, hy = ay + (by - ay) * m;
      var px = CX + r.sx + (X[j] - CX - r.sx) * m, py = CY - r.sy + (Y[j] - CY + r.sy) * m;
      ctx.lineWidth = BAR_W + (whiskerWidth(hMag, Math.max(0, Z[j])) - BAR_W) * m;
      ctx.beginPath();
      ctx.moveTo(px - hx, py - hy);
      ctx.lineTo(px + hx, py + hy);
      ctx.stroke();
    }
    ctx.globalAlpha = 1;
  }

  // Near hemisphere only, far to near; opacity and width as in the workshop sphere.
  function drawSegments(ctx, m) {
    var half = 0.5 * WLEN;
    ctx.strokeStyle = "#000";
    ctx.fillStyle = "#000";
    ctx.lineCap = "round";
    for (var i = 0; i < NG; i++) {
      var j = order[i], z = Z[j];
      if (z < 0) continue;
      evalField(j);
      var u = hMag;
      var dx = half * u * (cPsi * TX[j] + sPsi * PX[j]);
      var dy = half * u * (cPsi * TY[j] + sPsi * PY[j]);
      ctx.globalAlpha = (0.55 + 0.45 * z) * m;
      ctx.lineWidth = whiskerWidth(u, z);
      ctx.beginPath();
      ctx.moveTo(X[j] - dx, Y[j] - dy);
      ctx.lineTo(X[j] + dx, Y[j] + dy);
      ctx.stroke();
      if (u < 0.02) {                     // a nearly-null point is still a point
        ctx.beginPath();
        ctx.arc(X[j], Y[j], 1.8, 0, 2 * Math.PI);
        ctx.fill();
      }
    }
    ctx.globalAlpha = 1;
  }

  // ---- frame ----
  var t0 = 0;
  function draw(ctx, now) {
    var Phi = 2 * Math.PI * (now - t0) / PERIOD;
    c2Phi = Math.cos(2 * Phi); s2Phi = Math.sin(2 * Phi);
    var x = val(xf, now), m = val(morph, now);
    if (projDirty) reproject();
    ctx.setTransform(2, 0, 0, 2, 0, 0);
    ctx.clearRect(0, 0, 1920, 1080);
    if (m > 0) {                          // the sphere is transparent: a faint outline only
      ctx.globalAlpha = 0.3 * m;
      ctx.strokeStyle = "#000";
      ctx.lineWidth = 2;
      ctx.beginPath();
      ctx.arc(CX, CY, DIST, 0, 2 * Math.PI);
      ctx.stroke();
    }
    drawHoles(ctx, Phi);
    if (x < 1) drawRings(ctx, 1 - x);
    if (m > 0) drawSegments(ctx, m);
    if (x > 0 && m < 1) drawBars(ctx, x, m);
  }

  var raf = null, ctx = null, canvas = null, stage = 0, drag = null;
  function loop(now) { draw(ctx, now); raf = requestAnimationFrame(loop); }
  function stop() { if (raf) cancelAnimationFrame(raf); raf = null; }
  function resetCamera() { yaw = YAW0; pitch = PITCH0; projDirty = true; }
  function endDrag() {
    drag = null;
    canvas.classList.remove("sphere-grabbing");
  }

  // ---- drag, as in the workshop sphere; deltas are converted to slide pixels ----
  function wire() {
    canvas.addEventListener("pointerdown", function (e) {
      if (stage !== 2) return;
      canvas.setPointerCapture(e.pointerId);
      drag = { x: e.clientX, y: e.clientY, k: canvas.getBoundingClientRect().width / 1920 };
      canvas.classList.add("sphere-grabbing");
    });
    canvas.addEventListener("pointermove", function (e) {
      if (!drag) return;
      var dx = (e.clientX - drag.x) / drag.k, dy = (e.clientY - drag.y) / drag.k;
      drag.x = e.clientX; drag.y = e.clientY;
      yaw += dx * 0.008;
      pitch = Math.max(-1.55, Math.min(1.55, pitch + dy * 0.008));
      projDirty = true;
    });
    canvas.addEventListener("pointerup", endDrag);
    canvas.addEventListener("pointercancel", endDrag);
  }

  Deck.widget("sphere-slide", {
    steps: 2,
    enter: function (slide) {
      if (!canvas) {
        canvas = slide.querySelector("canvas");
        ctx = canvas.getContext("2d");
        wire();
      }
      t0 = performance.now();
      set(xf, 0, t0, true); set(morph, 0, t0, true);
      resetCamera();
      stop();
      raf = requestAnimationFrame(loop);
    },
    leave: function () { stop(); if (canvas) endDrag(); },
    step: function (slide, k, dir) {
      var now = performance.now(), inst = dir < 0;
      stage = k;
      set(xf, k >= 1 ? 1 : 0, now, inst, XF_MS);
      set(morph, k >= 2 ? 1 : 0, now, inst, MORPH_MS);
      if (k < 2) { endDrag(); resetCamera(); }
      canvas.classList.toggle("sphere-grab", k === 2);
    },
    // Not used by the engine: the closed form as the slide computes it, for the node check.
    _check: {
      constants: { DIST: DIST, ORB: ORB, BH: BH, R0: R0, DOT: DOT, DISP: DISP, BAR_HALF: BAR_HALF,
                   WLEN: WLEN, WSCALE: WSCALE, NG: NG, PITCH0: PITCH0, YAW0: YAW0, PERIOD: PERIOD },
      field: function (theta, phi, Phi) {
        var save = NP - 1, ct = Math.cos(theta), st = Math.sin(theta);
        var keep = [Nx, Ny, Nz, Tx3, Ty3, Tz3, Px3, Py3, FA, FB, C2, S2].map(function (a) { return a[save]; });
        addPoint(save, ct, st, Math.cos(phi), Math.sin(phi));
        c2Phi = Math.cos(2 * Phi); s2Phi = Math.sin(2 * Phi);
        evalField(save);
        var out = { abs: hMag, c: cPsi, s: sPsi };
        [Nx, Ny, Nz, Tx3, Ty3, Tz3, Px3, Py3, FA, FB, C2, S2].forEach(function (a, i) { a[save] = keep[i]; });
        return out;
      }
    }
  });
})();
