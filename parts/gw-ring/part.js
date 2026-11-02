// Ring of test masses deformed by the gravitational wave of a binary black hole.
// Stage 0: circular ring. Stage 1: two black holes appear left and right.
// Stage 2: ring stretched towards the holes, squeezed top and bottom.
// Stage 3: the holes orbit (speeding up to the final rate); the ring follows as for a
//          circularly polarised wave.
// Stage 4: the view turns to 3D: the orbital plane moves to the left, seen from the side,
//          and its ring fades out.
// Stage 5: a ring in a plane in front of the binary, along the orbital axis, with a delay.
// Stage 6: a second ring, further out, with twice the delay.
// Stage 7: five more rings, appearing one after the other.
// Stage 8: the orbital plane tilts to contain the propagation axis (edge-on binary);
//          new rings appear one after the other, now plus polarised.
(function () {
  "use strict";

  // World coordinates in slide pixels: X right, Y up in the orbital plane of stages 0-7,
  // Z along the orbital axis (the propagation direction of the waves drawn).
  var R0 = 300;             // ring radius
  var N = 24;               // number of test masses per ring
  var DOT = 11;             // test-mass radius, screen px
  var BH = 45;              // black-hole radius, screen px
  var D2 = 450;             // orbital radius in the 2D view
  var D3 = 220;             // orbital radius in the 3D view (clear of the first ring plane)
  var EPS = 0.26;           // strain amplitude h shown on the rings
  var PERIOD = 4000;        // final orbital period, ms
  var SLOW = 0.35;          // initial orbital rate as a fraction of the final rate
  var RAMP = 1500;          // time to reach the final rate, ms
  var NW = 7;               // number of ring planes in front of the binary
  var DZ = 400;             // spacing of the ring planes along Z
  var DELAY = PERIOD / 14;  // delay per ring plane, ms (7 planes = one wave period)
  var FB = 70;              // margin of the orbital-plane frame beyond the orbital radius
  var FW = 360;             // half-size of a ring-plane frame
  // 3D view: orthographic camera, yaw ALPHA from the Z axis (sin = 0.4: planes seen from
  // the side), pitch BETA, scale S3, image of the origin at (CX3, CY3).
  var ALPHA = Math.asin(0.4), BETA = 12 * Math.PI / 180, S3 = 0.55, CX3 = 330, CY3 = 600;

  // ---- animated parameters: value eases from `from` to `to` over [t0, t0 + dur] ----
  function ease(f) { return f < 0.5 ? 2 * f * f : 1 - 2 * (1 - f) * (1 - f); }
  function P(v) { return { from: v, to: v, t0: 0, dur: 1 }; }
  function val(p, now) {
    var f = Math.max(0, Math.min(1, (now - p.t0) / p.dur));
    return p.from + (p.to - p.from) * ease(f);
  }
  function set(p, target, now, instant, delay, dur) {
    if (instant) { p.from = p.to = target; return; }
    if (p.to === target) return;
    p.from = val(p, now); p.to = target; p.t0 = now + delay; p.dur = dur;
  }

  var par = { bh: P(0), h: P(0), ring: P(1), cam: P(0), frame: P(0), tilt: P(0) };
  var circ = [], plus = [];
  for (var j = 0; j < NW; j++) { circ.push(P(0)); plus.push(P(0)); }

  // ---- orbit: phase of the binary at time s (ms) since the orbit started ----
  // Rate w(s) = wf - (wf - wi)(1 - s/RAMP)^2 for s < RAMP, wf after; phase = its integral.
  var WF = 2 * Math.PI / PERIOD, WI = SLOW * WF;
  var orbitStart = null;    // performance.now() at the start of the orbit; null: no orbit
  // phase(s) is shared by both slides; s <= 0 means the orbit has not started.
  function phase(s) {
    if (s <= 0) return 0;
    var u = Math.min(s, RAMP) / RAMP;
    return WF * s - (WF - WI) * (RAMP / 3) * (1 - Math.pow(1 - u, 3));
  }

  // ---- projection ----
  var cam = {};
  function setCamera(c) {
    var a = Math.PI / 2 + (ALPHA - Math.PI / 2) * c, b = BETA * c;
    var B0 = [-Math.cos(a), 0, Math.sin(a)];
    cam.R = [Math.sin(a), 0, Math.cos(a)];
    cam.U = [Math.sin(b) * B0[0], Math.cos(b), Math.sin(b) * B0[2]];
    cam.s = 1 + (S3 - 1) * c;
    cam.cx = 960 + (CX3 - 960) * c;
    cam.cy = 540 + (CY3 - 540) * c;
  }
  function proj(x, y, z) {
    var R = cam.R, U = cam.U;
    return [cam.cx + cam.s * (x * R[0] + y * R[1] + z * R[2]),
            cam.cy - cam.s * (x * U[0] + y * U[1] + z * U[2])];
  }

  // Square frame of half-size a in the plane through o spanned by e1, e2.
  function frame(ctx, o, e1, e2, a, alpha) {
    if (alpha <= 0) return;
    ctx.beginPath();
    [[-1, -1], [1, -1], [1, 1], [-1, 1]].forEach(function (q, i) {
      var p = proj(o[0] + a * (q[0] * e1[0] + q[1] * e2[0]),
                   o[1] + a * (q[0] * e1[1] + q[1] * e2[1]),
                   o[2] + a * (q[0] * e1[2] + q[1] * e2[2]));
      if (i === 0) ctx.moveTo(p[0], p[1]); else ctx.lineTo(p[0], p[1]);
    });
    ctx.closePath();
    ctx.fillStyle = "rgba(0,0,0," + 0.04 * alpha + ")";
    ctx.fill();
    ctx.strokeStyle = "rgba(0,0,0," + 0.35 * alpha + ")";
    ctx.lineWidth = 2;
    ctx.stroke();
  }

  // Ring of test masses in the plane Z = z. Each mass moves by dx_i = (1/2) h_ij x_j with
  // h_+ = h cos 2phi and h_x = h sin 2phi (circular) or h_x = 0 (plus).
  function ring(ctx, z, h, phi, circular, alpha, r) {
    if (alpha <= 0) return;
    var hp = h * Math.cos(2 * phi), hx = circular ? h * Math.sin(2 * phi) : 0;
    ctx.globalAlpha = alpha;
    ctx.fillStyle = "#000";
    for (var k = 0; k < N; k++) {
      var th = 2 * Math.PI * k / N;
      var x = R0 * Math.cos(th), y = R0 * Math.sin(th);
      var p = proj(x + 0.5 * (hp * x + hx * y), y + 0.5 * (hx * x - hp * y), z);
      ctx.beginPath();
      ctx.arc(p[0], p[1], r, 0, 2 * Math.PI);
      ctx.fill();
    }
    ctx.globalAlpha = 1;
  }

  // Faint lines through the centre of the plane through o along e1 and e2.
  function planeAxes(ctx, o, e1, e2, a, alpha) {
    if (alpha <= 0) return;
    ctx.strokeStyle = "rgba(0,0,0," + 0.2 * alpha + ")";
    ctx.lineWidth = 1.5;
    [e1, e2].forEach(function (e) {
      var p = proj(o[0] - a * e[0], o[1] - a * e[1], o[2] - a * e[2]);
      var q = proj(o[0] + a * e[0], o[1] + a * e[1], o[2] + a * e[2]);
      ctx.beginPath(); ctx.moveTo(p[0], p[1]); ctx.lineTo(q[0], q[1]); ctx.stroke();
    });
  }

  // Small faint x, y, z triad in the lower left, turned with the camera.
  var TRIAD = [170, 930], TL = 110;
  function triad(ctx, alpha) {
    if (alpha <= 0) return;
    var col = "rgba(0,0,0," + 0.4 * alpha + ")";
    ctx.strokeStyle = col; ctx.fillStyle = col; ctx.lineWidth = 2;
    ctx.font = "italic 30px KaTeX_Math, 'Times New Roman', serif";
    ctx.textAlign = "center"; ctx.textBaseline = "middle";
    [["x", [1, 0, 0]], ["y", [0, 1, 0]], ["z", [0, 0, 1]]].forEach(function (ax) {
      var v = ax[1];
      var dx = v[0] * cam.R[0] + v[1] * cam.R[1] + v[2] * cam.R[2];
      var dy = -(v[0] * cam.U[0] + v[1] * cam.U[1] + v[2] * cam.U[2]);
      var n = Math.hypot(dx, dy);
      var x1 = TRIAD[0] + TL * dx, y1 = TRIAD[1] + TL * dy;
      ctx.beginPath(); ctx.moveTo(TRIAD[0], TRIAD[1]); ctx.lineTo(x1, y1); ctx.stroke();
      if (n > 0.05) {
        var ux = dx / n, uy = dy / n;
        ctx.beginPath();
        ctx.moveTo(x1, y1);
        ctx.lineTo(x1 - 12 * ux - 6 * uy, y1 - 12 * uy + 6 * ux);
        ctx.lineTo(x1 - 12 * ux + 6 * uy, y1 - 12 * uy - 6 * ux);
        ctx.closePath(); ctx.fill();
        ctx.fillText(ax[0], x1 + 24 * ux, y1 + 24 * uy);
      }
    });
  }

  function draw(ctx, now) {
    ctx.setTransform(2, 0, 0, 2, 0, 0);
    ctx.clearRect(0, 0, 1920, 1080);
    var c = val(par.cam, now);
    setCamera(c);
    var s = orbitStart === null ? 0 : now - orbitStart;
    var h = val(par.h, now), size = 1 - 0.2 * c;
    var t = val(par.tilt, now) * Math.PI / 2;
    var e1 = [1, 0, 0], e2 = [0, Math.cos(t), Math.sin(t)];   // orbital plane

    var D = D2 + (D3 - D2) * c;
    frame(ctx, [0, 0, 0], e1, e2, D + FB, val(par.frame, now));
    planeAxes(ctx, [0, 0, 0], e1, e2, D + FB, val(par.frame, now));
    triad(ctx, val(par.frame, now));
    for (var j = 0; j < NW; j++) {
      var z = (j + 1) * DZ, ret = phase(s - (j + 1) * DELAY);
      var ac = val(circ[j], now), ap = val(plus[j], now);
      frame(ctx, [0, 0, z], [1, 0, 0], [0, 1, 0], FW, Math.max(ac, ap));
      planeAxes(ctx, [0, 0, z], [1, 0, 0], [0, 1, 0], FW, Math.max(ac, ap));
      ring(ctx, z, h, ret, true, ac, DOT * size);
      ring(ctx, z, h, ret, false, ap, DOT * size);
    }
    ring(ctx, 0, h, phase(s), true, val(par.ring, now), DOT * size);

    var a = val(par.bh, now);
    if (a > 0) {
      var phi = phase(s);
      ctx.globalAlpha = a;
      ctx.fillStyle = "#000";
      for (var i = 0; i < 2; i++) {
        var sg = i === 0 ? 1 : -1, cp = sg * D * Math.cos(phi), sp = sg * D * Math.sin(phi);
        var p = proj(cp * e1[0] + sp * e2[0], cp * e1[1] + sp * e2[1], cp * e1[2] + sp * e2[2]);
        ctx.beginPath();
        ctx.arc(p[0], p[1], BH * size, 0, 2 * Math.PI);
        ctx.fill();
      }
      ctx.globalAlpha = 1;
    }
  }

  var raf = null, ctx = null;
  function loop(now) { draw(ctx, now); raf = requestAnimationFrame(loop); }
  function stop() { if (raf) cancelAnimationFrame(raf); raf = null; }

  Deck.widget("gw-ring-binary", {
    steps: 8,
    enter: function (slide) {
      ctx = slide.querySelector("canvas").getContext("2d");
      orbitStart = null;
      stop();
      raf = requestAnimationFrame(loop);
    },
    leave: function () { stop(); orbitStart = null; },
    step: function (slide, k, dir) {
      var now = performance.now(), inst = dir < 0;
      if (k < 3) orbitStart = null;
      // The orbit starts slowly on the way forward into stage 3; otherwise at full rate.
      else if (orbitStart === null) orbitStart = k === 3 && !inst ? now : now - 100 * PERIOD;
      set(par.bh, k >= 1 ? 1 : 0, now, inst, 0, 700);
      set(par.h, k >= 2 ? EPS : 0, now, inst, 0, 1200);
      set(par.ring, k <= 3 ? 1 : 0, now, inst, 0, 600);
      set(par.cam, k >= 4 ? 1 : 0, now, inst, 400, 2000);
      set(par.frame, k >= 4 ? 1 : 0, now, inst, 400, 1500);
      set(par.tilt, k >= 8 ? 1 : 0, now, inst, 500, 1500);
      for (var j = 0; j < NW; j++) {
        var on = (k === 5 && j === 0) || (k === 6 && j <= 1) || k === 7;
        set(circ[j], on ? 1 : 0, now, inst, k === 7 ? Math.max(0, j - 2) * 400 : 0, on ? 600 : 400);
        set(plus[j], k >= 8 ? 1 : 0, now, inst, 2000 + j * 350, 500);
      }
    }
  });
})();
