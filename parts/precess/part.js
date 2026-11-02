// Equal-mass binary on a circular Newtonian orbit of fixed radius, with spin precession.
// Slide precess-compare: stage 0 title only; stage 1 non-precessing (left) and precessing
// (right) side by side, same camera and same orbital rate.
// Slide precess-detectors: the precessing binary alone, with one small LIGO along J (face on)
// and one in the plane perpendicular to J (edge on).
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
  var PANELS = [480, 1440], CY_CMP = 600;
  var CX_DET = 960, CY_DET = 590;
  // Small LIGOs: arm length (slide px), distance from the binary.
  var ARM = 100, MASS = 9, DIST_J = 510, DIST_PERP = 800;
  var PERP_AZ = 38 * Math.PI / 180;   // position in the plane perpendicular to J, from camera right

  // ---- orthographic camera: right R, up U, towards the camera V (R x U = V) ----
  var EL = 38 * Math.PI / 180;        // elevation above the plane perpendicular to J
  var R = [1, 0, 0];
  var U = [0, Math.sin(EL), Math.cos(EL)];
  var V = [0, -Math.cos(EL), Math.sin(EL)];

  function proj(p, cx, cy, s) { return [cx + s * dot(p, R), cy - s * dot(p, U)]; }

  function square(ctx, e1, e2, cx, cy, s) {
    ctx.beginPath();
    [[-1, -1], [1, -1], [1, 1], [-1, 1]].forEach(function (q, i) {
      var p = proj([FA * (q[0] * e1[0] + q[1] * e2[0]), FA * (q[0] * e1[1] + q[1] * e2[1]),
                    FA * (q[0] * e1[2] + q[1] * e2[2])], cx, cy, s);
      if (i === 0) ctx.moveTo(p[0], p[1]); else ctx.lineTo(p[0], p[1]);
    });
    ctx.closePath();
    ctx.fillStyle = "rgba(0,0,0,0.04)";
    ctx.fill();
    ctx.strokeStyle = "rgba(0,0,0,0.35)";
    ctx.lineWidth = 2;
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

  // One binary: square, then far hole and its arrow, then near hole and its arrow. A hole's
  // arrow is always on top of its own disc, so that it cannot be hidden.
  function binary(ctx, y, phi, cx, cy, s) {
    var e1 = vec(y, 9), e2 = vec(y, 12);
    square(ctx, e1, e2, cx, cy, s);
    var c = Math.cos(phi), d = Math.sin(phi);
    var pos = [c * e1[0] + d * e2[0], c * e1[1] + d * e2[1], c * e1[2] + d * e2[2]];
    var near = dot(pos, V) > 0 ? 1 : -1;     // +1: the hole at +pos is the nearer one
    [-near, near].forEach(function (sg) {
      var P0 = [sg * ORB * pos[0], sg * ORB * pos[1], sg * ORB * pos[2]];
      var Si = vec(y, sg > 0 ? 3 : 6);
      var a = proj(P0, cx, cy, s);
      ctx.fillStyle = "#000";
      ctx.beginPath();
      ctx.arc(a[0], a[1], BH * s, 0, 2 * Math.PI);
      ctx.fill();
      arrow(ctx, a, proj(combo(P0, SPIN_PX, Si), cx, cy, s), s);
    });
  }

  // Small L: corner and two end masses, centred on the point c, in the plane of u, v.
  function ligo(ctx, c, u, v, cx, cy) {
    var o = combo(combo(c, -ARM / 3, u), -ARM / 3, v);
    var p0 = proj(o, cx, cy, 1), pu = proj(combo(o, ARM, u), cx, cy, 1), pv = proj(combo(o, ARM, v), cx, cy, 1);
    ctx.strokeStyle = "#000";
    ctx.lineWidth = 4;
    ctx.lineCap = "butt";
    ctx.beginPath();
    ctx.moveTo(pu[0], pu[1]); ctx.lineTo(p0[0], p0[1]); ctx.lineTo(pv[0], pv[1]);
    ctx.stroke();
    ctx.fillStyle = "#000";
    [p0, pu, pv].forEach(function (p) {
      ctx.beginPath(); ctx.arc(p[0], p[1], MASS, 0, 2 * Math.PI); ctx.fill();
    });
  }

  // ---- fixed non-precessing state: everything along z ----
  var FIXED = [0, 0, 1, 0, 0, SPIN, 0, 0, SPIN,
               Math.cos(YAW), Math.sin(YAW), 0, -Math.sin(YAW), Math.cos(YAW), 0];

  // Detector positions and arm directions (the plane of each L is perpendicular to the line
  // from the binary to the detector).
  var D_PERP = [Math.cos(PERP_AZ), -Math.sin(PERP_AZ), 0];   // towards the camera side
  var LIGO_J = { c: [0, 0, DIST_J], u: [1, 0, 0], v: [0, 1, 0] };
  var LIGO_PERP = { c: [DIST_PERP * D_PERP[0], DIST_PERP * D_PERP[1], 0],
                    u: [Math.sin(PERP_AZ), Math.cos(PERP_AZ), 0], v: [0, 0, 1] };

  // ---- one widget: simulation clock, canvas, one rAF loop ----
  // running(k): is the animation on at stage k. paint(ctx, y, phi): draw one frame.
  function makeWidget(id, steps, startStage, paint) {
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
    binary(ctx, FIXED, phi, PANELS[0], CY_CMP, S_CMP);
    binary(ctx, y, phi, PANELS[1], CY_CMP, S_CMP);
  });

  makeWidget("precess-detectors", 0, 0, function (ctx, y, phi) {
    ligo(ctx, LIGO_J.c, LIGO_J.u, LIGO_J.v, CX_DET, CY_DET);
    ligo(ctx, LIGO_PERP.c, LIGO_PERP.u, LIGO_PERP.v, CX_DET, CY_DET);
    binary(ctx, y, phi, CX_DET, CY_DET, 1);
  });
})();
