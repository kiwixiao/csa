"""Prescribed flow / mask-pressure profiles vs what the CFD export actually did (plain python3).

Inputs: the two STAR profile tables (time, value) and the per-step boundary series written by
cfd_boundary_series.py. The script
  * drops exact duplicate export steps (same flow, cap area and mask pressure as an earlier step),
  * rebuilds a uniform time axis (dt per step) if the export's time labels are not uniform,
  * fits density from  mdot_table(t) = rho * Q_cfd(t),
  * plots flow (STAR sign: + = into the domain at the bronchi = exhale) and mask pressure.

Usage: python3 plot_profiles_vs_cfd.py --subject ENT004 --flow F.csv --pressure P.csv --cfd series.csv
                                       --out fig.png [--dt 0.025] [--sim-end 2.183] [--note "..."]
"""
import argparse
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

INK, INK2, GRID, SURF = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
S1, S2 = '#2a78d6', '#eb6834'          # categorical slot 1 = prescribed table, slot 2 = CFD result

ap = argparse.ArgumentParser()
ap.add_argument('--subject', required=True); ap.add_argument('--flow', required=True); ap.add_argument('--pressure', required=True)
ap.add_argument('--cfd', required=True); ap.add_argument('--out', required=True)
ap.add_argument('--dt', type=float, default=0.025, help='export interval (s)')
ap.add_argument('--sim-end', type=float, default=None, help='physical time the run actually stopped at (s)')
ap.add_argument('--note', default='', help='extra line under the figure title')
ap.add_argument('--steps', default=None, help='use only these export steps, "a:b" (b exclusive) - e.g. when the export holds two runs')
a = ap.parse_args()

F = np.loadtxt(a.flow, delimiter=',', skiprows=1); P = np.loadtxt(a.pressure, delimiter=',', skiprows=1)
R = list(csv.DictReader(open(a.cfd)))
if a.steps:
    lo, hi = (int(x) for x in a.steps.split(':'))
    R = [r for r in R if lo <= int(r['step']) < hi]
    print('using export steps %d..%d (%d rows)' % (lo, hi - 1, len(R)))
col = lambda k: np.array([float(r[k]) for r in R])
step, lab, qb, ab, pm = col('step').astype(int), col('label_t_s'), col('Q_bronchi_out_Lmin'), col('A_bronchi_mm2'), col('p_mask_Pa')

# --- duplicates + time axis
# A restarted run re-computes physical times it already saved. Those steps share the IMPOSED outlet
# flow and the morphed cap area exactly (both set by the tables), while interior fields differ slightly
# (LES). Key on the imposed quantities, and keep the LATER copy: it belongs to the run that carried on.
last = {}
for i in range(len(step)):
    last[(round(qb[i], 5), round(ab[i], 5))] = i
keep = np.array(sorted(last.values()))
dups = [(step[i], step[last[(round(qb[i], 5), round(ab[i], 5))]]) for i in range(len(step)) if i not in set(keep)]
uniform = np.allclose(np.diff(lab), a.dt, atol=1e-6) and not dups
t = lab[keep] if uniform else a.dt * (np.arange(len(keep)) + 1)
if dups:
    print('dropped %d earlier copies of re-computed steps (kept the later run), e.g. %s'
          % (len(dups), ', '.join('%d->%d' % d for d in dups[:4])))
if not uniform:
    print('time labels non-uniform -> rebuilt as %.3f*(n+1): %d unique steps, %.3f..%.3f s' % (a.dt, len(t), t[0], t[-1]))
Q = -qb[keep]                                           # STAR sign: + = into the domain at the bronchi = exhale
p_cfd = pm[keep]

m_at = np.interp(t, F[:, 0], F[:, 1])
rho = float(np.sum(m_at * Q / 6e4) / np.sum((Q / 6e4) ** 2))
err = np.abs(m_at / rho * 6e4 - Q).max()
perr = np.abs(np.interp(t, P[:, 0], P[:, 1]) - p_cfd).max()
print('fitted density %.4f kg/m3; max |flow table - CFD| %.3f L/min; max |pressure table - CFD| %.2f Pa' % (rho, err, perr))
q = F[:, 1] / rho * 6e4

# --- phases of the TABLE as STAR applies it
tt = np.linspace(F[0, 0], F[-1, 0], 40001); qq = np.interp(tt, F[:, 0], q); dtt = tt[1] - tt[0]
sg = np.sign(np.where(np.abs(qq) < 1e-9, 0, qq))
nz = np.flatnonzero(sg); flips = [tt[nz[j]] for j in range(1, len(nz)) if sg[nz[j]] != sg[nz[j - 1]]]
edges = [tt[0]] + flips + [tt[-1]]
t_end = edges[-1]
t_last = a.sim_end if a.sim_end is not None else t[-1]
gap = t_end - t_last > 2 * a.dt

plt.rcParams.update({'font.size': 11, 'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2,
                     'axes.titleweight': 'bold', 'axes.titlesize': 13, 'axes.titlelocation': 'left'})
fig, (a1, a2) = plt.subplots(2, 1, figsize=(12, 9.4), sharex=True, gridspec_kw=dict(height_ratios=[1.25, 1], hspace=0.34))
fig.patch.set_facecolor(SURF)
for ax in (a1, a2):
    ax.set_facecolor(SURF); ax.grid(color=GRID, lw=0.8); ax.set_axisbelow(True)
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
    ax.axhline(0, color=INK2, lw=0.9)
    if gap:
        ax.axvspan(t_last, t_end, facecolor='none', edgecolor='#b9b8b2', hatch='///', lw=0)
        ax.axvline(t_last, color=INK, lw=1.2, ls='--')

# --- flow panel
a1.plot(F[:, 0], q, color=S1, lw=2, label='Prescribed flow at BronchiOut (table in the .sim)')
a1.plot(t, Q, 'o', ms=4.5, color=S2, mec=SURF, mew=0.8, label='CFD export: flow through BronchiOut, every %.0f ms' % (a.dt * 1e3))
qmax = np.abs(q).max()
for lo, hi in zip(edges[:-1], edges[1:]):
    m = (tt >= lo) & (tt <= hi); seg = qq[m]
    if not m.any() or np.abs(seg).max() < 0.05 * qmax: continue
    exh = seg[np.argmax(np.abs(seg))] > 0
    vol = np.abs(seg).sum() * dtt / 6e4 * 1e6
    sim = np.abs(seg[tt[m] <= t_last]).sum() * dtt / 6e4 * 1e6
    txt = '%s  (%.2f - %.2f s)\n%.0f mL %s' % ('EXHALE' if exh else 'INHALE', lo, hi, vol, 'out' if exh else 'in')
    if gap and hi > t_last: txt += ' (only %.0f mL simulated)' % sim
    a1.text((lo + hi) / 2, -0.3 * qmax if exh else 0.3 * qmax, txt, ha='center', va='center', color=INK2, fontweight='bold')
    k = np.argmax(np.abs(seg)); tk = tt[m][k]
    a1.annotate('peak %s %.1f L/min @ %.2f s' % ('exhale' if exh else 'inhale', seg[k], tk), (tk, seg[k]),
                xytext=(0, 12 if exh else -14), textcoords='offset points', ha='center', va='bottom' if exh else 'top', color=INK)
if gap:
    a1.text(t_last + 0.02, -1.05 * qmax, 'run stopped at %.3f s\n-> rest of the breath\n    never simulated' % t_last, color=INK, fontsize=10, va='bottom')
a1.set_ylabel('flow (L/min)\n+ = exhale (into domain at bronchi)')
a1.set_ylim(-1.32 * qmax, 1.32 * qmax)
a1.set_title('%s prescribed flow vs the CFD   (density fitted from the CFD: %.3f kg/m$^3$, max error %.2f L/min)' % (a.subject, rho, err), color=INK, pad=30)
a1.legend(loc='lower left', bbox_to_anchor=(0, 1.0), ncol=2, frameon=False, labelcolor=INK, fontsize=10, borderaxespad=0.1)

# --- pressure panel
a2.plot(P[:, 0], P[:, 1], color=S1, lw=2, label='Prescribed mask pressure, as STAR applied it (table %.2f - %.2f s)' % (P[0, 0], P[-1, 0]))
off = P[0, 0] - F[0, 0]
if abs(off) > 1e-6:
    a2.plot(P[:, 0] - off, P[:, 1], color=INK2, lw=1.4, ls=(0, (4, 3)), label='Same table with its time axis zeroed (%+.2f s)' % -off)
a2.plot(t, p_cfd, 'o', ms=4.5, color=S2, mec=SURF, mew=0.8, label='CFD export: mean static pressure on MaskInlet')
i_in, i_ex = np.argmin(q), np.argmax(q)
p_in, p_ex = np.interp(F[i_in, 0], P[:, 0], P[:, 1]), np.interp(F[i_ex, 0], P[:, 0], P[:, 1])
ok = p_in < 0 < p_ex
verdict = ('correct sign: %+.0f Pa at peak inhale, %+.0f Pa at peak exhale' if ok else
           'SIGN REVERSED vs the flow: %+.0f Pa at peak inhale, %+.0f Pa at peak exhale') % (p_in, p_ex)
if abs(off) > 1e-6: verdict += ', time axis offset %+.2f s' % off
pr = max(np.abs(P[:, 1]).max(), np.abs(p_cfd).max())
a2.set_ylim(-1.2 * pr, 1.2 * pr)
a2.set_ylabel('mask pressure (Pa, gauge)\n- = below ambient')
a2.set_xlabel('time (s)'); a2.set_xlim(0, max(F[-1, 0], P[-1, 0]) + 0.05)
a2.set_title('%s mask pressure: %s' % (a.subject, verdict), color=INK if ok else '#b3261e', pad=48)
a2.legend(loc='lower left', bbox_to_anchor=(0, 1.0), ncol=2, frameon=False, labelcolor=INK, fontsize=10, borderaxespad=0.1)
if a.note:
    fig.suptitle(a.note, x=0.06, ha='left', y=0.995, fontsize=12, color='#b3261e', fontweight='bold')
fig.savefig(a.out, dpi=130, bbox_inches='tight', facecolor=SURF)
print('wrote', a.out)
