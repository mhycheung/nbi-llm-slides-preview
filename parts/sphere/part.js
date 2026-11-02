// Non-precessing equal-mass binary and its quadrupole strain, in phase with the orbit.
// Everything is in one 3D scene with x, y, z axes; dragging turns the camera, at every stage.
// Stage 0: the binary and four rings of test masses at +z, -z, +x, -x. Each ring lies in the plane
// perpendicular to the line from the binary to it, so it faces the binary.
// Stage 1: each ring becomes a bar along the stretch axis in that plane, length proportional to |h|.
// Stage 2: the bars shrink to segments of a sphere of 600. The binary, its dotted orbit and its
// square are inside the sphere; the surface is shaded by log |h|, white at the four nodes on the
// equator, where the wave vanishes at that instant.
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
  var ELEV = 28 * Math.PI / 180;       // camera elevation above the orbital plane, first view
  var YAW0 = -35 * Math.PI / 180, PITCH0 = ELEV - Math.PI / 2;   // the rings at +-x are not edge on to the camera
  var PLANE = 150;          // half-size of the faint square around each ring, slide px
  var ORB_SQ = 190;         // half-size of the orbital square; its corners, at 269 px, are inside the sphere
  var AXIS = 470;           // axis half-length, slide px; outside the rings and the sphere
  var AXIS_W = 3, AXIS_HEAD = 24, AXIS_COL = "#444", LABEL_OFF = 36;
  var DOT_W = 3.5, DOT_GAP = 11;       // dotted orbit: dot size and spacing, slide px
  // Surface shading: black with opacity SH_MAX * s, s = 1 + log|h| / log(1 / H_MIN), clipped to 0..1.
  // |h| = 1 at the poles gives s = 1, |h| <= H_MIN gives s = 0: white, on the white slide.
  var H_MIN = 0.05, SH_MAX = 0.32;
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
  var shA = G(), shB = G(), shC = G(), shS = G();
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
        var k = j * SH_N + i, rho2 = x * x + y * y;
        shA[k] = (1 + z * z) / 2;
        shB[k] = z;
        // cos 2phi, sin 2phi; at the pole itself phi has no meaning and |h| does not depend on it.
        shC[k] = rho2 > 1e-12 ? (x * x - y * y) / rho2 : 1;
        shS[k] = rho2 > 1e-12 ? 2 * x * y / rho2 : 0;
      }
    }
  }
  function shadeAlpha(h2) {   // opacity 0..1 from |h|^2
    if (h2 <= 0) return 0;
    var s = 1 + 0.5 * Math.log(h2) / Math.log(1 / H_MIN);
    return SH_MAX * (s < 0 ? 0 : s > 1 ? 1 : s);
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
      var hp = shA[k] * (c2Phi * shC[k] + s2Phi * shS[k]), hx = shB[k] * (s2Phi * shC[k] - c2Phi * shS[k]);
      d[4 * k + 3] = 255 * shadeAlpha(hp * hp + hx * hx);
    }
    shCtx.putImageData(shImg, 0, 0);
    ctx.save();
    ctx.globalAlpha = m;
    ctx.beginPath();
    ctx.arc(CX, CY, DIST, 0, 2 * Math.PI);
    ctx.clip();
    ctx.drawImage(shCanvas, CX - DIST, CY - DIST, 2 * DIST, 2 * DIST);
    ctx.restore();
  }

  function whiskerWidth(u, z) { return WSCALE * (1.2 + 2.3 * u) * (0.7 + 0.3 * z); }

  // ---- the parts of a frame ----
  // The orbital plane: a faint square, the dotted orbit, and the holes, far hole first.
  function drawBinary(ctx, Phi) {
    ctx.globalAlpha = 1;
    ctx.beginPath();
    [[-1, -1], [1, -1], [1, 1], [-1, 1]].forEach(function (q, i) {
      var px = sx(ORB_SQ * q[0], ORB_SQ * q[1], 1), py = sy(ORB_SQ * q[0], ORB_SQ * q[1], 0, 1);
      if (i === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
    });
    ctx.closePath();
    ctx.fillStyle = "rgba(0,0,0,0.05)";
    ctx.fill();
    ctx.strokeStyle = "rgba(0,0,0,0.4)";
    ctx.lineWidth = 2;
    ctx.stroke();

    ctx.strokeStyle = "#000";
    ctx.lineWidth = DOT_W;
    ctx.lineCap = "round";
    ctx.setLineDash([0.01, DOT_GAP]);
    ctx.beginPath();
    for (var i = 0; i < 120; i++) {
      var a = 2 * Math.PI * i / 120, x = ORB * Math.cos(a), y = ORB * Math.sin(a);
      if (i === 0) ctx.moveTo(sx(x, y, 1), sy(x, y, 0, 1)); else ctx.lineTo(sx(x, y, 1), sy(x, y, 0, 1));
    }
    ctx.closePath();
    ctx.stroke();
    ctx.setLineDash([]);

    var cp = ORB * Math.cos(Phi), sp = ORB * Math.sin(Phi);
    ctx.fillStyle = "#000";
    ctx.beginPath();
    ctx.arc(sx(cp, sp, 1), sy(cp, sp, 0, 1), BH, 0, 2 * Math.PI);
    ctx.moveTo(sx(-cp, -sp, 1) + BH, sy(-cp, -sp, 0, 1));
    ctx.arc(sx(-cp, -sp, 1), sy(-cp, -sp, 0, 1), BH, 0, 2 * Math.PI);
    ctx.fill();
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
  // The halves of the axes towards (near) or away from the camera, between radii r0 and r1.
  function drawAxes(ctx, near, r0, r1) {
    AXES.forEach(function (a) {
      var sg = depth(a.d[0], a.d[1], a.d[2]) >= 0 ? 1 : -1;
      axisPart(ctx, a.d, near ? sg : -sg, r0, r1);
    });
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

  // Unit stretch axis at ring point j, for the field last evaluated, is cPsi e_theta + sPsi e_phi.
  // Faint square of ring r, in its plane.
  function drawPlane(ctx, r, alpha) {
    var j = r.j, c = DIST;
    ctx.globalAlpha = alpha;
    ctx.beginPath();
    [[-1, -1], [1, -1], [1, 1], [-1, 1]].forEach(function (q, i) {
      var x = c * Nx[j] + PLANE * (q[0] * Tx3[j] + q[1] * Px3[j]);
      var y = c * Ny[j] + PLANE * (q[0] * Ty3[j] + q[1] * Py3[j]);
      var z = c * Nz[j] + PLANE * q[0] * Tz3[j];
      if (i === 0) ctx.moveTo(sx(x, y, 1), sy(x, y, z, 1)); else ctx.lineTo(sx(x, y, 1), sy(x, y, z, 1));
    });
    ctx.closePath();
    ctx.fillStyle = "rgba(0,0,0,0.035)";
    ctx.fill();
    ctx.globalAlpha = 1;
  }

  // In the plane of the ring a mass at p = a e_theta + b e_phi moves by dx_i = (1/2) h_ij p_j, with
  // h_ij = DISP |h| (2 n_i n_j - delta_ij) and n the stretch axis, so in the same plane.
  function drawRing(ctx, r, alpha) {
    var j = r.j;
    evalField(j);
    var hd = 0.5 * DISP * hMag;
    ctx.globalAlpha = alpha;
    ctx.fillStyle = "#000";
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
    ctx.strokeStyle = "#000";
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

    // One ring's objects, for those behind or in front of the binary.
    function rings(near) {
      RINGS.forEach(function (r) {
        if ((Z[r.j] >= 0) !== near) return;
        if (m < 1) drawPlane(ctx, r, 1 - m);
        if (x < 1) drawRing(ctx, r, 1 - x);
        if (x > 0 && m < 1) drawBar(ctx, r, x, m);
      });
    }

    // Far to near: far axes, the inside, near axes, near surface, segments. Nothing inside the
    // sphere is drawn over the near surface or its segments.
    drawAxes(ctx, false, 0, AXIS);
    rings(false);
    drawBinary(ctx, Phi);
    drawAxes(ctx, true, 0, m > 0 ? DIST : AXIS);
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
      drawAxes(ctx, true, DIST, AXIS);
    }
    rings(true);
    placeLabels();
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
      canvas.setPointerCapture(e.pointerId);
      drag = { x: e.clientX, y: e.clientY, k: canvas.getBoundingClientRect().width / 1920 };
      canvas.classList.add("sphere-grabbing");
    });
    canvas.addEventListener("pointermove", function (e) {
      if (!drag) return;
      var dx = (e.clientX - drag.x) / drag.k, dy = (e.clientY - drag.y) / drag.k;
      drag.x = e.clientX; drag.y = e.clientY;
      yaw += dx * 0.008;
      // Pitch 0 looks down the z axis from above, -pi/2 from the side, -pi from below.
      pitch = Math.max(-Math.PI, Math.min(0, pitch + dy * 0.008));
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
        labels = {};
        AXES.forEach(function (a) { labels[a.name] = slide.querySelector(".sphere-axis-" + a.name); });
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
      canvas.classList.add("sphere-grab");
    },
    // Not used by the engine: the closed form as the slide computes it, for the node check.
    _check: {
      camera: function () { return { yaw: yaw, pitch: pitch }; },
      setCamera: function (y, p) { yaw = y; pitch = p; projDirty = true; },
      shadeAlpha: shadeAlpha,
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
