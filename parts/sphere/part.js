// Non-precessing equal-mass binary and its quadrupole strain, in phase with the orbit.
// Everything is in one 3D scene with x, y, z axes; dragging turns the camera, at every stage.
// Stage 0: the binary and four rings of test masses at +z, -z, +x, -x. Each ring lies in the plane
// perpendicular to the line from the binary to it, so it faces the binary, on a nearly opaque square.
// Stage 1: each ring becomes a bar along the stretch axis in that plane, length proportional to |h|.
// Stage 2: the bars shrink to segments of a sphere of 600. The binary, its dotted orbit and its
// square are inside the sphere, seen through its surface; the surface is shaded by log |h|, white at
// the four nodes on the equator, where the wave vanishes at that instant.
//
// Slide sphere-precess: the precessing binary and its axes; at the next stage the same sphere around it.
// Clicking a viewing direction (without dragging) highlights its ring or bar; clicking again removes it.
//
// Field, the one formula for everything drawn. With e1, e2 the orbit's in-plane axes and
// Q = (1/2)[cos 2Phi (e1 e1 - e2 e2) + sin 2Phi (e1 e2 + e2 e1)], projected on e_theta, e_phi, normalised by 2:
//   h_+ = Q_tt - Q_pp,  h_x = 2 Q_tp
// For e1 = x, e2 = y this is h_+ = (1 + cos^2 theta)/2 cos 2(Phi - phi), h_x = cos theta sin 2(Phi - phi),
// so |h| = 1 at the poles and at most 1/2 on the equator. In general |h|^2 = 1 - 4 |Q n|^2 + (n . Q n)^2.
// The stretch axis is cos(psi) e_theta + sin(psi) e_phi with psi = (1/2) atan2(h_x, h_+).
// With precession this is the instantaneous quadrupole of the tilted orbit, a schematic: no retardation.
(function () {
  "use strict";

  // ---- PHYSICS BEGIN (tests extract this block; no DOM below this line until PHYSICS END) ----
  // Orbit-averaged equal-mass precession, coupling exaggerated, no spin-spin term:
  //   dL/dt = k S x L,  dS_i/dt = k L x S_i,  S = S1 + S2,  J = L + S conserved.
  // L then precesses about J at rate k|J|. World axes: z along J.
  // Units of L, S: |L| = 1. Time in seconds.
  var PERIOD_S = 4;          // orbital period, s (as in gw-ring)
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
  var K = 2 * Math.PI / (N_PREC * PERIOD_S * JMAG);
  var DT = 1 / 120;         // fixed RK4 step, s
  var OMEGA = 2 * Math.PI / PERIOD_S;
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
  var ELEV = 28 * Math.PI / 180;       // camera elevation above the orbital plane, first view
  var YAW0 = -35 * Math.PI / 180, PITCH0 = ELEV - Math.PI / 2;   // the rings at +-x are not edge on to the camera
  var PLANE = 150;          // half-size of the faint square around each ring, slide px
  var ORB_SQ = 190;         // half-size of the orbital square; its corners, at 269 px, are inside the sphere
  var AXIS = 470;           // axis half-length, slide px; outside the rings and the sphere
  var AXIS_W = 2, AXIS_HEAD = 18, AXIS_COL = "#b8b8b8", LABEL_OFF = 34;
  var PLATE_FILL = "rgba(233,233,233,0.82)", PLATE_EDGE = "#c4c4c4";   // nearly opaque, to hide what is behind
  var DOT_W = 3.5, DOT_GAP = 11;       // dotted orbit: dot size and spacing, slide px
  // Surface shading: gray level 1 - SH_MAX * s, s = 1 + log|h| / log(1 / H_MIN), clipped to 0..1.
  // |h| = 1 at the poles gives s = 1, |h| <= H_MIN gives s = 0: white. SURF is the surface's opacity,
  // so that what is inside the sphere is visible but clearly behind the surface.
  var H_MIN = 0.12, SH_MAX = 0.34, SURF = 0.66;
  var SPIN_PX = 103, SHAFT = 4, HEAD = 15, COL = "#d00000";   // spin arrows, px per unit of spin; highlight colour
  var PICK_MOVE = 6, PICK_R = 28, PICK_RING = 110;   // a click moves less than this; reach of a click, slide px
  var SH_N = 256;           // shading image, pixels per side; smooth, so scaled up to the sphere
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
  function addPoint(j, ct, st, cp, sp) {
    Nx[j] = st * cp; Ny[j] = st * sp; Nz[j] = ct;
    Tx3[j] = ct * cp; Ty3[j] = ct * sp; Tz3[j] = -st;
    Px3[j] = -sp; Py3[j] = cp;
  }
  var GOLDEN = (1 + Math.sqrt(5)) / 2;
  for (var g = 0; g < NG; g++) {
    var gz = 1 - (2 * g + 1) / NG, gphi = 2 * Math.PI * g / GOLDEN;
    addPoint(g, gz, Math.sqrt(1 - gz * gz), Math.cos(gphi), Math.sin(gphi));
  }
  // Ring directions. Exact values, so that h_x = 0 on the equator and the poles are exact. The
  // ring's plane is spanned by e_theta, e_phi at its point, which is perpendicular to the radius.
  var RINGS = [
    { j: NG,     ct:  1, st: 0, cp:  0, sp:  1 },
    { j: NG + 1, ct: -1, st: 0, cp:  0, sp: -1 },
    { j: NG + 2, ct:  0, st: 1, cp:  1, sp:  0 },
    { j: NG + 3, ct:  0, st: 1, cp: -1, sp:  0 }
  ];
  RINGS.forEach(function (r) { addPoint(r.j, r.ct, r.st, r.cp, r.sp); });
  var MC = new Float64Array(NM), MS = new Float64Array(NM);
  for (var m0 = 0; m0 < NM; m0++) {
    MC[m0] = Math.cos(2 * Math.PI * m0 / NM);
    MS[m0] = Math.sin(2 * Math.PI * m0 / NM);
  }

  // ---- field at point j for the current orbital phase ----
  // The orbit now: in-plane axes, normal, and the quadrupole of the phase.
  var E1 = [1, 0, 0], E2 = [0, 1, 0], LN = [0, 0, 1];
  var Qxx = 0, Qxy = 0, Qxz = 0, Qyy = 0, Qyz = 0, Qzz = 0;
  function setOrbit(e1, e2, Phi) {
    E1 = e1; E2 = e2; LN = cross(e1, e2);
    var c = 0.5 * Math.cos(2 * Phi), d = 0.5 * Math.sin(2 * Phi);
    Qxx = c * (e1[0] * e1[0] - e2[0] * e2[0]) + 2 * d * e1[0] * e2[0];
    Qyy = c * (e1[1] * e1[1] - e2[1] * e2[1]) + 2 * d * e1[1] * e2[1];
    Qzz = c * (e1[2] * e1[2] - e2[2] * e2[2]) + 2 * d * e1[2] * e2[2];
    Qxy = c * (e1[0] * e1[1] - e2[0] * e2[1]) + d * (e1[0] * e2[1] + e1[1] * e2[0]);
    Qxz = c * (e1[0] * e1[2] - e2[0] * e2[2]) + d * (e1[0] * e2[2] + e1[2] * e2[0]);
    Qyz = c * (e1[1] * e1[2] - e2[1] * e2[2]) + d * (e1[1] * e2[2] + e1[2] * e2[1]);
  }
  var hMag = 0, cPsi = 1, sPsi = 0;
  function evalField(j) {
    var tx = Tx3[j], ty = Ty3[j], tz = Tz3[j], px = Px3[j], py = Py3[j];
    var ax = Qxx * tx + Qxy * ty + Qxz * tz, ay = Qxy * tx + Qyy * ty + Qyz * tz, az = Qxz * tx + Qyz * ty + Qzz * tz;
    var bx = Qxx * px + Qxy * py, by = Qxy * px + Qyy * py;
    var hp = tx * ax + ty * ay + tz * az - (px * bx + py * by);
    var hx = 2 * (px * ax + py * ay);
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
    shadeTables();
    projDirty = false;
  }

  // World vector -> screen, at slide px per unit s.
  function sx(x, y, s) { return CX + s * (x * cyw - y * syw); }
  function sy(x, y, z, s) { return CY - s * ((x * syw + y * cyw) * cpi - z * spi); }
  function depth(x, y, z) { return (x * syw + y * cyw) * spi + z * cpi; }

  // ---- surface shading: an image of the near hemisphere, on the sphere's disc ----
  // For each pixel the field coefficients of the surface point under it. The far hemisphere is not
  // shaded: seen through a node of the near one it would make that node gray.
  var NPX = SH_N * SH_N;
  function G() { return new Float32Array(NPX); }
  var shX = G(), shY = G(), shZ = G();
  var shCanvas = null, shCtx = null, shImg = null;
  function shadeTables() {
    for (var j = 0; j < SH_N; j++) {
      for (var i = 0; i < SH_N; i++) {
        var u = (i + 0.5) / SH_N * 2 - 1, v = 1 - (j + 0.5) / SH_N * 2, q = u * u + v * v;
        if (q > 1) { var f = 1 / Math.sqrt(q); u *= f; v *= f; q = 1; }   // the rim, for smooth edges
        var w = Math.sqrt(1 - q);
        // Inverse of the camera: screen right u, up v, towards the camera w.
        var y1 = v * cpi + w * spi, z = -v * spi + w * cpi;
        var x = u * cyw + y1 * syw, y = -u * syw + y1 * cyw;
        var k = j * SH_N + i;
        shX[k] = x; shY[k] = y; shZ[k] = z;
      }
    }
  }
  function shadeLevel(h2) {   // gray level 0..1 from |h|^2
    if (h2 <= 0) return 1;
    var s = 1 + 0.5 * Math.log(h2) / Math.log(1 / H_MIN);
    return 1 - SH_MAX * (s < 0 ? 0 : s > 1 ? 1 : s);
  }
  function drawShade(ctx, m) {
    if (!shCanvas) {
      shCanvas = document.createElement("canvas");
      shCanvas.width = shCanvas.height = SH_N;
      shCtx = shCanvas.getContext("2d");
      shImg = shCtx.createImageData(SH_N, SH_N);
    }
    var d = shImg.data;
    for (var k = 0; k < NPX; k++) {
      var x = shX[k], y = shY[k], z = shZ[k];
      var a = Qxx * x + Qxy * y + Qxz * z, b = Qxy * x + Qyy * y + Qyz * z, c = Qxz * x + Qyz * y + Qzz * z;
      var nQn = x * a + y * b + z * c;
      var h2 = 1 - 4 * (a * a + b * b + c * c) + nQn * nQn;
      d[4 * k] = d[4 * k + 1] = d[4 * k + 2] = 255 * shadeLevel(h2);
      d[4 * k + 3] = 255;
    }
    shCtx.putImageData(shImg, 0, 0);
    ctx.save();
    ctx.globalAlpha = SURF * m;
    ctx.beginPath();
    ctx.arc(CX, CY, DIST, 0, 2 * Math.PI);
    ctx.clip();
    ctx.drawImage(shCanvas, CX - DIST, CY - DIST, 2 * DIST, 2 * DIST);
    ctx.restore();
  }

  function whiskerWidth(u, z) { return WSCALE * (1.2 + 2.3 * u) * (0.7 + 0.3 * z); }

  // ---- the parts of a frame ----
  // A plate: nearly opaque square.
  function plate(ctx, pts, alpha) {
    ctx.globalAlpha = alpha;
    ctx.beginPath();
    pts.forEach(function (p, i) { if (i === 0) ctx.moveTo(p[0], p[1]); else ctx.lineTo(p[0], p[1]); });
    ctx.closePath();
    ctx.fillStyle = PLATE_FILL;
    ctx.fill();
    ctx.strokeStyle = PLATE_EDGE;
    ctx.lineWidth = 1.5;
    ctx.stroke();
    ctx.globalAlpha = 1;
  }
  function inPlane(a, b) {   // the point a e1 + b e2
    return [a * E1[0] + b * E2[0], a * E1[1] + b * E2[1], a * E1[2] + b * E2[2]];
  }
  function screen(p) { return [sx(p[0], p[1], 1), sy(p[0], p[1], p[2], 1)]; }
  function drawOrbitalPlate(ctx) {
    plate(ctx, [[-1, -1], [1, -1], [1, 1], [-1, 1]].map(function (q) {
      return screen(inPlane(ORB_SQ * q[0], ORB_SQ * q[1]));
    }), 1);
  }
  // The dotted orbit and the holes.
  function drawOrbit(ctx, Phi) {
    ctx.globalAlpha = 1;
    ctx.strokeStyle = "#000";
    ctx.lineWidth = DOT_W;
    ctx.lineCap = "round";
    ctx.setLineDash([0.01, DOT_GAP]);
    ctx.beginPath();
    for (var i = 0; i < 120; i++) {
      var a = 2 * Math.PI * i / 120, p = screen(inPlane(ORB * Math.cos(a), ORB * Math.sin(a)));
      if (i === 0) ctx.moveTo(p[0], p[1]); else ctx.lineTo(p[0], p[1]);
    }
    ctx.closePath();
    ctx.stroke();
    ctx.setLineDash([]);

    var h = screen(inPlane(ORB * Math.cos(Phi), ORB * Math.sin(Phi)));
    ctx.fillStyle = "#000";
    ctx.beginPath();
    ctx.arc(h[0], h[1], BH, 0, 2 * Math.PI);
    ctx.moveTo(2 * CX - h[0] + BH, 2 * CY - h[1]);
    ctx.arc(2 * CX - h[0], 2 * CY - h[1], BH, 0, 2 * Math.PI);
    ctx.fill();
  }

  // Spin arrows of the precessing binary, those behind the plate or those in front of it.
  function drawSpins(ctx, Phi, y, behind) {
    var side = depth(LN[0], LN[1], LN[2]) >= 0 ? 1 : -1;
    [1, -1].forEach(function (sg) {
      var S = vec(y, sg > 0 ? 3 : 6);
      if ((dot(S, LN) * side < 0) !== behind) return;
      var P0 = inPlane(sg * ORB * Math.cos(Phi), sg * ORB * Math.sin(Phi));
      var a = screen(P0), b = screen(combo(P0, SPIN_PX, S));
      var dx = b[0] - a[0], dy = b[1] - a[1], n = Math.hypot(dx, dy);
      if (n < 1e-6) return;
      var ux = dx / n, uy = dy / n, h = Math.min(HEAD, 0.5 * n), w = 0.45 * h, ex = b[0] - h * ux, ey = b[1] - h * uy;
      ctx.globalAlpha = 1;
      ctx.strokeStyle = COL;
      ctx.fillStyle = COL;
      ctx.lineWidth = SHAFT;
      ctx.lineCap = "butt";
      ctx.beginPath(); ctx.moveTo(a[0], a[1]); ctx.lineTo(ex, ey); ctx.stroke();
      ctx.beginPath();
      ctx.moveTo(b[0], b[1]);
      ctx.lineTo(ex - w * uy, ey + w * ux);
      ctx.lineTo(ex + w * uy, ey - w * ux);
      ctx.closePath();
      ctx.fill();
    });
  }

  // Half of an axis, from radius r0 to r1 along sg times the unit vector d; arrow head at r = AXIS
  // on the + half.
  function axisPart(ctx, d, sg, r0, r1) {
    var ax = sx(sg * r0 * d[0], sg * r0 * d[1], 1), ay = sy(sg * r0 * d[0], sg * r0 * d[1], sg * r0 * d[2], 1);
    var bx = sx(sg * r1 * d[0], sg * r1 * d[1], 1), by = sy(sg * r1 * d[0], sg * r1 * d[1], sg * r1 * d[2], 1);
    var dx = bx - ax, dy = by - ay, n = Math.hypot(dx, dy);
    if (n < 1e-6) return;
    var ux = dx / n, uy = dy / n, head = sg > 0 && r1 === AXIS;
    var h = head ? Math.min(AXIS_HEAD, 0.5 * n) : 0, w = 0.42 * h, ex = bx - h * ux, ey = by - h * uy;
    ctx.globalAlpha = 1;
    ctx.strokeStyle = AXIS_COL;
    ctx.fillStyle = AXIS_COL;
    ctx.lineWidth = AXIS_W;
    ctx.lineCap = "butt";
    ctx.beginPath(); ctx.moveTo(ax, ay); ctx.lineTo(ex, ey); ctx.stroke();
    if (head) {
      ctx.beginPath();
      ctx.moveTo(bx, by);
      ctx.lineTo(ex - w * uy, ey + w * ux);
      ctx.lineTo(ex + w * uy, ey - w * ux);
      ctx.closePath();
      ctx.fill();
    }
  }
  var AXES = [{ d: [1, 0, 0], name: "x" }, { d: [0, 1, 0], name: "y" }, { d: [0, 0, 1], name: "z" }];
  // The halves of the axes between radii r0 and r1 that pass test(sign along the axis, axis direction).
  function drawAxes(ctx, r0, r1, test) {
    AXES.forEach(function (a) {
      [1, -1].forEach(function (sg) { if (test(sg, a.d)) axisPart(ctx, a.d, sg, r0, r1); });
    });
  }
  function isNear(sg, d) { return sg * depth(d[0], d[1], d[2]) >= 0; }
  // Which side of the orbital plate a half axis is on, against the camera's: negative is behind it.
  function plateSide(sg, d) {
    return sg * dot(d, LN) * (depth(LN[0], LN[1], LN[2]) >= 0 ? 1 : -1);
  }
  // Beside the arrow tip: to the right of a vertical axis, above a horizontal one.
  var labels = null;
  function placeLabels() {
    AXES.forEach(function (a) {
      var tx = sx(AXIS * a.d[0], AXIS * a.d[1], 1), ty = sy(AXIS * a.d[0], AXIS * a.d[1], AXIS * a.d[2], 1);
      var dx = tx - CX, dy = ty - CY, n = Math.hypot(dx, dy) || 1, px = dy / n, py = -dx / n;
      if (py > 0.2 || (py > -0.2 && px < 0)) { px = -px; py = -py; }
      labels[a.name].style.left = (tx + LABEL_OFF * px) + "px";
      labels[a.name].style.top = (ty + LABEL_OFF * py) + "px";
    });
  }

  // Square of ring r, in its plane.
  function drawPlane(ctx, r, alpha) {
    var j = r.j, c = DIST;
    plate(ctx, [[-1, -1], [1, -1], [1, 1], [-1, 1]].map(function (q) {
      var x = c * Nx[j] + PLANE * (q[0] * Tx3[j] + q[1] * Px3[j]);
      var y = c * Ny[j] + PLANE * (q[0] * Ty3[j] + q[1] * Py3[j]);
      var z = c * Nz[j] + PLANE * q[0] * Tz3[j];
      return [sx(x, y, 1), sy(x, y, z, 1)];
    }), alpha);
  }

  // In the plane of the ring a mass at p = a e_theta + b e_phi moves by dx_i = (1/2) h_ij p_j, with
  // h_ij = DISP |h| (2 n_i n_j - delta_ij) and n the stretch axis, so in the same plane.
  function drawRing(ctx, r, alpha) {
    var j = r.j;
    evalField(j);
    var hd = 0.5 * DISP * hMag;
    ctx.globalAlpha = alpha;
    ctx.fillStyle = hl[j] ? COL : "#000";
    ctx.beginPath();
    for (var k = 0; k < NM; k++) {
      var a = R0 * MC[k], b = R0 * MS[k], ap = 2 * (a * cPsi + b * sPsi);
      a += hd * (ap * cPsi - a);
      b += hd * (ap * sPsi - b);
      var x = DIST * Nx[j] + a * Tx3[j] + b * Px3[j];
      var y = DIST * Ny[j] + a * Ty3[j] + b * Py3[j];
      var z = DIST * Nz[j] + a * Tz3[j];
      var px = sx(x, y, 1), py = sy(x, y, z, 1);
      ctx.moveTo(px + DOT, py);
      ctx.arc(px, py, DOT, 0, 2 * Math.PI);
    }
    ctx.fill();
    ctx.globalAlpha = 1;
  }

  // Bar along the stretch axis in the ring's plane, shrinking to the segment's length as m -> 1.
  function drawBar(ctx, r, alpha, m) {
    var j = r.j, fade = 1 - Math.max(0, Math.min(1, (m - 0.6) / 0.4));
    if (alpha * fade <= 0) return;
    evalField(j);
    var half = (BAR_HALF + (0.5 * WLEN * DIST - BAR_HALF) * m) * hMag;
    var ax = half * (cPsi * Tx3[j] + sPsi * Px3[j]), ay = half * (cPsi * Ty3[j] + sPsi * Py3[j]), az = half * cPsi * Tz3[j];
    var x = DIST * Nx[j], y = DIST * Ny[j], z = DIST * Nz[j];
    ctx.globalAlpha = alpha * fade;
    ctx.strokeStyle = hl[j] ? COL : "#000";
    ctx.lineCap = "butt";
    ctx.lineWidth = BAR_W + (whiskerWidth(hMag, Math.max(0, Z[j])) - BAR_W) * m;
    ctx.beginPath();
    ctx.moveTo(sx(x - ax, y - ay, 1), sy(x - ax, y - ay, z - az, 1));
    ctx.lineTo(sx(x + ax, y + ay, 1), sy(x + ax, y + ay, z + az, 1));
    ctx.stroke();
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
      if (z < 0 || hl[j]) continue;
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
    // Highlighted directions on the near side, on top: a bar in colour, and a dot that stays when |h| = 0.
    ctx.strokeStyle = COL;
    ctx.fillStyle = COL;
    Object.keys(hl).forEach(function (key) {
      var j = +key, z = Z[j];
      if (z < 0) return;
      evalField(j);
      var dx = half * hMag * (cPsi * TX[j] + sPsi * PX[j]), dy = half * hMag * (cPsi * TY[j] + sPsi * PY[j]);
      ctx.globalAlpha = m;
      ctx.lineWidth = 2 * whiskerWidth(hMag, z);
      ctx.beginPath();
      ctx.moveTo(X[j] - dx, Y[j] - dy);
      ctx.lineTo(X[j] + dx, Y[j] + dy);
      ctx.stroke();
      ctx.beginPath();
      ctx.arc(X[j], Y[j], 6, 0, 2 * Math.PI);
      ctx.fill();
    });
    ctx.globalAlpha = 1;
  }

  // ---- frame ----
  var t0 = 0, prec = false, hl = {};
  // Precessing state and its clock, as on the precess slides.
  var yS = initialState(), nS = 0, accS = 0, lastS = null;
  function advance(now) {
    if (lastS !== null) accS += Math.min(0.25, (now - lastS) / 1000);
    lastS = now;
    while (accS >= DT) { yS = rk4(yS, DT); nS++; accS -= DT; }
  }
  function draw(ctx, now) {
    var Phi;
    if (prec) {
      advance(now);
      Phi = OMEGA * nS * DT;
      setOrbit(vec(yS, 9), vec(yS, 12), Phi);
    } else {
      Phi = 2 * Math.PI * (now - t0) / PERIOD;
      setOrbit([1, 0, 0], [0, 1, 0], Phi);
    }
    var x = val(xf, now), m = val(morph, now);
    if (projDirty) reproject();
    ctx.setTransform(2, 0, 0, 2, 0, 0);
    ctx.clearRect(0, 0, 1920, 1080);

    // The rings behind or in front of the binary, far to near. These plates do not meet each other
    // or the orbital plate, so this order is correct from every direction.
    function rings(near) {
      if (prec) return;
      RINGS.filter(function (r) { return (Z[r.j] >= 0) === near; })
        .sort(function (a, b) { return Z[a.j] - Z[b.j]; })
        .forEach(function (r) {
          if (m < 1) drawPlane(ctx, r, 1 - m);
          if (x < 1) drawRing(ctx, r, 1 - x);
          if (x > 0 && m < 1) drawBar(ctx, r, x, m);
        });
    }

    // Back to front, so that each plate and the surface hide what is behind them. An axis half
    // behind the orbital plate goes before it; the others go after it, and over the holes if near.
    var EPS_SIDE = 1e-9;
    function behind(sg, d) { return plateSide(sg, d) < -EPS_SIDE; }
    drawAxes(ctx, DIST, AXIS, function (sg, d) { return !isNear(sg, d); });
    rings(false);
    drawAxes(ctx, 0, DIST, behind);
    if (prec) drawSpins(ctx, Phi, yS, true);
    drawOrbitalPlate(ctx);
    drawAxes(ctx, 0, DIST, function (sg, d) { return !behind(sg, d) && !isNear(sg, d); });
    drawOrbit(ctx, Phi);
    if (prec) drawSpins(ctx, Phi, yS, false);
    drawAxes(ctx, 0, DIST, function (sg, d) { return !behind(sg, d) && isNear(sg, d); });
    if (m > 0) {
      drawShade(ctx, m);
      ctx.globalAlpha = 0.35 * m;         // the sphere's outline
      ctx.strokeStyle = "#000";
      ctx.lineWidth = 2;
      ctx.beginPath();
      ctx.arc(CX, CY, DIST, 0, 2 * Math.PI);
      ctx.stroke();
      ctx.globalAlpha = 1;
      drawSegments(ctx, m);
    }
    rings(true);
    drawAxes(ctx, DIST, AXIS, isNear);
    placeLabels();
  }

  var raf = null, ctx = null, canvas = null, stage = 0, drag = null, wired = [];
  function loop(now) { draw(ctx, now); raf = requestAnimationFrame(loop); }
  function stop() { if (raf) cancelAnimationFrame(raf); raf = null; }
  function resetCamera() { yaw = YAW0; pitch = PITCH0; projDirty = true; }
  function endDrag() {
    drag = null;
    if (canvas) canvas.classList.remove("sphere-grabbing");
  }

  // A click: the nearest direction on the sphere, on its near side, or the nearest ring or bar.
  function pick(px, py) {
    var best = -1, bestD = Infinity, j, d;
    if ((prec && stage >= 1) || (!prec && stage >= 2)) {
      for (j = 0; j < NG; j++) {
        if (Z[j] < 0) continue;
        d = Math.hypot(X[j] - px, Y[j] - py);
        if (d < bestD && d < PICK_R) { bestD = d; best = j; }
      }
    } else if (!prec) {
      RINGS.forEach(function (r) {
        d = Math.hypot(X[r.j] - px, Y[r.j] - py);
        if (d < bestD && d < PICK_RING) { bestD = d; best = r.j; }
      });
    }
    if (best < 0) return;
    if (hl[best]) delete hl[best]; else hl[best] = true;
  }

  // ---- drag, as in the workshop sphere; deltas are converted to slide pixels ----
  function wire(c) {
    c.addEventListener("pointerdown", function (e) {
      c.setPointerCapture(e.pointerId);
      drag = { x: e.clientX, y: e.clientY, k: c.getBoundingClientRect().width / 1920, moved: 0 };
      c.classList.add("sphere-grabbing");
    });
    c.addEventListener("pointermove", function (e) {
      if (!drag) return;
      var dx = (e.clientX - drag.x) / drag.k, dy = (e.clientY - drag.y) / drag.k;
      drag.x = e.clientX; drag.y = e.clientY;
      drag.moved += Math.hypot(dx, dy);
      if (drag.moved < PICK_MOVE) return;
      yaw += dx * 0.008;
      // Pitch 0 looks down the z axis from above, -pi/2 from the side, -pi from below.
      pitch = Math.max(-Math.PI, Math.min(0, pitch + dy * 0.008));
      projDirty = true;
    });
    c.addEventListener("pointerup", function (e) {
      if (drag && drag.moved < PICK_MOVE) {
        var r = c.getBoundingClientRect();
        if (projDirty) reproject();
        pick((e.clientX - r.left) / drag.k, (e.clientY - r.top) / drag.k);
      }
      endDrag();
    });
    c.addEventListener("pointercancel", endDrag);
  }

  function enterScene(slide, isPrec) {
    canvas = slide.querySelector("canvas");
    ctx = canvas.getContext("2d");
    if (wired.indexOf(canvas) < 0) { wire(canvas); wired.push(canvas); }
    labels = {};
    AXES.forEach(function (a) { labels[a.name] = slide.querySelector(".sphere-axis-" + a.name); });
    prec = isPrec;
    hl = {};
    t0 = performance.now();
    yS = initialState(); nS = 0; accS = 0; lastS = null;
    set(xf, isPrec ? 1 : 0, t0, true); set(morph, 0, t0, true);
    resetCamera();
    canvas.classList.add("sphere-grab");
    stop();
    raf = requestAnimationFrame(loop);
  }

  Deck.widget("sphere-slide", {
    steps: 2,
    enter: function (slide) { enterScene(slide, false); },
    leave: function () { stop(); endDrag(); },
    step: function (slide, k, dir) {
      var now = performance.now(), inst = dir < 0;
      stage = k;
      set(xf, k >= 1 ? 1 : 0, now, inst, XF_MS);
      set(morph, k >= 2 ? 1 : 0, now, inst, MORPH_MS);
    },
    // Not used by the engine: the closed form as the slide computes it, and the state, for the checks.
    _check: {
      camera: function () { return { yaw: yaw, pitch: pitch }; },
      setCamera: function (y, p) { yaw = y; pitch = p; projDirty = true; },
      shadeLevel: shadeLevel,
      highlighted: function () { return Object.keys(hl).map(Number); },
      screenOf: function (j) { if (projDirty) reproject(); return [X[j], Y[j], Z[j]]; },
      constants: { DIST: DIST, ORB: ORB, BH: BH, R0: R0, DOT: DOT, DISP: DISP, BAR_HALF: BAR_HALF,
                   WLEN: WLEN, WSCALE: WSCALE, NG: NG, PITCH0: PITCH0, YAW0: YAW0, PERIOD: PERIOD },
      field: function (theta, phi, Phi) {
        var save = NP - 1, ct = Math.cos(theta), st = Math.sin(theta);
        var keep = [Nx, Ny, Nz, Tx3, Ty3, Tz3, Px3, Py3].map(function (a) { return a[save]; });
        addPoint(save, ct, st, Math.cos(phi), Math.sin(phi));
        setOrbit([1, 0, 0], [0, 1, 0], Phi);
        evalField(save);
        var out = { abs: hMag, c: cPsi, s: sPsi };
        [Nx, Ny, Nz, Tx3, Ty3, Tz3, Px3, Py3].forEach(function (a, i) { a[save] = keep[i]; });
        return out;
      }
    }
  });

  Deck.widget("sphere-precess", {
    steps: 1,
    enter: function (slide) { enterScene(slide, true); },
    leave: function () { stop(); endDrag(); },
    step: function (slide, k, dir) {
      stage = k;
      set(morph, k >= 1 ? 1 : 0, performance.now(), dir < 0, MORPH_MS);
    },
    _check: {
      state: function () { return { y: yS, Phi: OMEGA * nS * DT }; },
      highlighted: function () { return Object.keys(hl).map(Number); }
    }
  });
})();
