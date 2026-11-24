#!/usr/bin/env python3
"""Work of breathing by anatomical region, inhale vs exhale.

For each CFD time step, a region's dissipated power = (total-pressure drop across it) x (flow
measured through it). TP is the mass-flow-weighted plane average (area average where |Q| is too
small), interpolated at the region ends; the flow is the mean measured flow of the region's own
planes, so wall motion is accounted for.
  * shared regions (pharynx, larynx, trachea) are bounded by distances along the --main-path
    (mm from its nostril opening); they lie on the shared airway, so one path covers them.
  * the NOSE is two parallel channels: its power is the sum over both nostrils of
    (drop from that nostril's opening to the top of the shared airway) x (that nostril's flow).
Energy per region and phase = sum over steps of power x dt (inhale: prescribed mass flow < 0,
STAR convention; exhale: > 0). Writes a CSV and a figure: bars (regions x inhale/exhale, mJ) and
a region x time map of dissipated power (when the load occurs).
Usage: wob_regions.py --pe-dir DIR --prof flow.csv --rho R --main-path LeftNose
                      --regions "pharynx:nose_end:150,larynx:150:195,trachea:195:end" --subject S --out-dir D
"""
import argparse, csv, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import csa_paths as CP

ap = argparse.ArgumentParser()
ap.add_argument('--pe-dir', required=True); ap.add_argument('--prof', required=True); ap.add_argument('--rho', type=float, required=True)
ap.add_argument('--main-path', default='LeftNose'); ap.add_argument('--other-path', default='RightNose')
ap.add_argument('--regions', required=True, help='name:start:end,... in mm from the main nostril; "nose_end" = top of the shared airway, "end" = outlet')
ap.add_argument('--subject', default=''); ap.add_argument('--out-dir', required=True); ap.add_argument('--qmin', type=float, default=0.5)
ap.add_argument('--bin-mm', type=float, default=0, help='also write the WOB per fixed-length segment along the main path (0 = off)')
ap.add_argument('--peak-time', type=float, default=None, help='also write the local loss (mW/mm) along the main path at this time')
ap.add_argument('--series', action='store_true', help='also write the local loss for every time step (loss_per_mm/) and the '
                'segment power vs time (for the report video)')
a = ap.parse_args()
os.makedirs(a.out_dir, exist_ok=True)

def path_data(path):
    times, D = CP.load_pe(a.pe_dir, path, qmin=a.qmin)
    out = {}
    for t in times:
        d = D[t]; o = np.argsort(d['s']); s = d['s'][o]
        tp = np.where(np.isfinite(d['tp_m'][o]), d['tp_m'][o], d['tp_a'][o])
        out[round(t, 4)] = dict(s=s, tp=tp, q=-d['Q'][o], reg=d['region'][o])        # q > 0 = toward the trachea (inhale)
    return out

M = path_data(a.main_path); O = path_data(a.other_path)
times = sorted(set(M) & set(O))
F = np.loadtxt(a.prof, delimiter=',', skiprows=1)
mdot = np.interp(times, F[:, 0], F[:, 1]); inh = mdot < 0
dt = float(np.median(np.diff(times)))

def seg_power(d, lo, hi):
    """(TP(lo) - TP(hi)) x mean measured flow of the planes in [lo, hi] (interpolated if none)."""
    k = (d['s'] >= lo) & (d['s'] <= hi)
    q = d['q'][k].mean() if k.any() else np.interp(0.5 * (lo + hi), d['s'], d['q'])
    return (np.interp(lo, d['s'], d['tp']) - np.interp(hi, d['s'], d['tp'])) * q

def nose_power(d):
    """(TP at the opening - TP at the top of the shared airway) x the nostril's mean flow."""
    nose = d['reg'] != 'DescendingAirway'
    s_top = d['s'][~nose].min()                         # first shared-airway plane
    q = d['q'][nose].mean()
    return (d['tp'][0] - np.interp(s_top, d['s'], d['tp'])) * q, s_top

spec = []
for item in a.regions.split(','):
    nm, s0, s1 = item.split(':'); spec.append((nm, s0, s1))
names = ['nose'] + [s[0] for s in spec]
P = np.zeros((len(names), len(times)))
s_top_main = []
for j, t in enumerate(times):
    d = M[t]; pn_m, st = nose_power(d); pn_o, _ = nose_power(O[t]); s_top_main.append(st)
    P[0, j] = pn_m + pn_o
    for i, (nm, s0, s1) in enumerate(spec, start=1):
        lo = st if s0 == 'nose_end' else float(s0); hi = d['s'].max() if s1 == 'end' else float(s1)
        P[i, j] = seg_power(d, lo, hi)
P *= 1e-3 / 60.0 * 1e3          # Pa x L/min -> W -> mW  (1 L/min = 1/60000 m3/s)
E_in = P[:, inh].sum(1) * dt; E_ex = P[:, ~inh].sum(1) * dt        # mW x s = mJ
T_in, T_ex = inh.sum() * dt, (~inh).sum() * dt

# cross-check against the whole-path calculation (drop opening -> outlet x local flow, both nostrils)
tot = np.array([((M[t]['tp'][0] - M[t]['tp'][-1]) * M[t]['q'][M[t]['reg'] == 'DescendingAirway'].mean()) for t in times]) * 1e-3 / 60 * 1e3
print('nose ends (top of shared airway) at %.1f-%.1f mm from the %s opening' % (min(s_top_main), max(s_top_main), a.main_path))
print('%-10s %10s %10s %10s %12s %12s' % ('region', 'inhale mJ', 'exhale mJ', 'total mJ', 'inh mean mW', 'exh mean mW'))
for i, nm in enumerate(names):
    print('%-10s %10.2f %10.2f %10.2f %12.2f %12.2f' % (nm, E_in[i], E_ex[i], E_in[i] + E_ex[i], E_in[i] / T_in, E_ex[i] / T_ex))
print('%-10s %10.2f %10.2f %10.2f   (whole-path check: %.2f mJ)' % ('sum', E_in.sum(), E_ex.sum(), E_in.sum() + E_ex.sum(), tot.sum() * dt))

with open(os.path.join(a.out_dir, '%s_wob_by_region.csv' % a.subject), 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(['region', 'bounds_mm_from_%s_opening' % a.main_path, 'inhale_mJ', 'exhale_mJ', 'total_mJ', 'inhale_mean_mW', 'exhale_mean_mW'])
    s_end = np.median([M[t]['s'].max() for t in times])
    bounds = ['0-%.0f' % np.median(s_top_main)] + ['%s-%s' % (('%.0f' % np.median(s_top_main)) if s0 == 'nose_end' else s0, ('%.0f' % s_end) if s1 == 'end' else s1) for _, s0, s1 in spec]
    for i, nm in enumerate(names):
        w.writerow([nm, bounds[i], '%.3f' % E_in[i], '%.3f' % E_ex[i], '%.3f' % (E_in[i] + E_ex[i]), '%.3f' % (E_in[i] / T_in), '%.3f' % (E_ex[i] / T_ex)])
with open(os.path.join(a.out_dir, '%s_power_by_region_vs_time.csv' % a.subject), 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(['time_s', 'phase'] + ['%s_mW' % n for n in names])
    for j, t in enumerate(times):
        w.writerow(['%.4f' % t, 'inhale' if inh[j] else 'exhale'] + ['%.4f' % P[i, j] for i in range(len(names))])

# ---- fixed-length segments along the main path ----
# nose parts (s < top of the shared airway) add both nostrils, each measured from its own opening;
# shared-airway parts come from the main path only
if a.bin_mm > 0:
    s_end = np.median([M[t]['s'].max() for t in times])
    edges = np.arange(0.0, s_end, a.bin_mm)
    edges = np.append(edges[:-1] if s_end - edges[-1] < 0.3 * a.bin_mm else edges, s_end)
    PB = np.zeros((len(edges) - 1, len(times)))
    for j, t in enumerate(times):
        for d, shared in ((M[t], True), (O[t], False)):
            st = d['s'][d['reg'] == 'DescendingAirway'].min()
            for b in range(len(edges) - 1):
                lo, hi = edges[b], edges[b + 1]
                if lo < st: PB[b, j] += seg_power(d, lo, min(hi, st))
                if shared and hi > st: PB[b, j] += seg_power(d, max(lo, st), hi)
    PB *= 1e-3 / 60.0 * 1e3
    Bi, Be = PB[:, inh].sum(1) * dt, PB[:, ~inh].sum(1) * dt
    fn = os.path.join(a.out_dir, '%s_wob_by_%gmm.csv' % (a.subject, a.bin_mm))
    with open(fn, 'w', newline='') as fh:
        w = csv.writer(fh); w.writerow(['lo_mm', 'hi_mm', 'inhale_mJ', 'exhale_mJ', 'total_mJ'])
        for b in range(len(edges) - 1):
            w.writerow(['%.0f' % edges[b], '%.0f' % edges[b + 1], '%.3f' % Bi[b], '%.3f' % Be[b], '%.3f' % (Bi[b] + Be[b])])
    print('%g mm segments: inhale %.2f + exhale %.2f = %.2f mJ  -> %s' % (a.bin_mm, Bi.sum(), Be.sum(), Bi.sum() + Be.sum(), fn))
    if a.series:
        fn = os.path.join(a.out_dir, '%s_power_by_%gmm_vs_time.csv' % (a.subject, a.bin_mm))
        with open(fn, 'w', newline='') as fh:
            w = csv.writer(fh); w.writerow(['time_s', 'dt_s', 'phase'] + ['%.0f-%.0f_mW' % (edges[b], edges[b + 1]) for b in range(len(edges) - 1)])
            for j, t in enumerate(times):
                w.writerow(['%.4f' % t, '%.4f' % dt, 'inhale' if inh[j] else 'exhale'] + ['%.4f' % v for v in PB[:, j]])

# ---- local loss along the main path at one instant (for the colour-coded centerline) ----
def loss_per_mm(d):
    from scipy.ndimage import gaussian_filter1d
    g = np.arange(0.0, d['s'].max(), 0.5)
    tpg = gaussian_filter1d(np.interp(g, d['s'], d['tp']), 5.0, mode='nearest')     # sigma 2.5 mm
    qg = gaussian_filter1d(np.interp(g, d['s'], d['q']), 5.0, mode='nearest')
    return g, -np.gradient(tpg, g) * qg / 60.0                                       # Pa/mm x L/min -> mW/mm

if a.series:
    sd = os.path.join(a.out_dir, 'loss_per_mm'); os.makedirs(sd, exist_ok=True)
    for t in times:
        g, dens = loss_per_mm(M[t])
        np.savetxt(os.path.join(sd, 't%.3f.csv' % t), np.column_stack([g, dens]), delimiter=',',
                   header='s_mm,loss_mW_per_mm', comments='', fmt='%.4f')
    print('local loss for %d time steps -> %s' % (len(times), sd))

if a.peak_time is not None:
    tk = times[int(np.argmin(np.abs(np.array(times) - a.peak_time)))]; d = M[tk]
    g, dens = loss_per_mm(d)
    fn = os.path.join(a.out_dir, '%s_loss_per_mm_t%.3f.csv' % (a.subject, tk))
    np.savetxt(fn, np.column_stack([g, dens]), delimiter=',', header='s_mm,loss_mW_per_mm', comments='', fmt='%.4f')
    print('local loss at t=%.3f s: peak %.2f mW/mm at s=%.0f mm, integral %.1f mW (path power %.1f mW) -> %s'
          % (tk, dens.max(), g[np.argmax(dens)], np.trapz(dens, g), (d['tp'][0] - d['tp'][-1]) * d['q'][d['reg'] == 'DescendingAirway'].mean() / 60.0, fn))

# ---- figure ---------------------------------------------------------------------------
INK, INK2, GRID, SURF = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
C_IN, C_EX = '#2a78d6', '#eb6834'
plt.rcParams.update({'font.size': 14, 'font.weight': 'bold', 'axes.titleweight': 'bold', 'axes.labelweight': 'bold',
                     'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK, 'ytick.color': INK})
fig = plt.figure(figsize=(17, 7.6), facecolor=SURF)
gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.25], height_ratios=[0.32, 1], wspace=0.22, hspace=0.12, left=0.06, right=0.925, top=0.86, bottom=0.12)
ax = fig.add_subplot(gs[:, 0]); ax.set_facecolor(SURF)
x = np.arange(len(names)); bw = 0.38
b1 = ax.bar(x - bw / 2 - 0.01, E_in, bw, color=C_IN, label='inhale (%.2f s)' % T_in)
b2 = ax.bar(x + bw / 2 + 0.01, E_ex, bw, color=C_EX, label='exhale (%.2f s)' % T_ex)
for bars, vals in ((b1, E_in), (b2, E_ex)):
    for bb, v in zip(bars, vals):
        ax.text(bb.get_x() + bb.get_width() / 2, v + 0.012 * max(E_ex.max(), E_in.max()), '%.1f' % v, ha='center', va='bottom', fontsize=13, color=INK)
lab = [n + '\n%s mm' % b for n, b in zip(names, bounds)]
ax.set_xticks(x); ax.set_xticklabels(lab, fontsize=13)
ax.set_ylabel('work of breathing (mJ)'); ax.grid(axis='y', color=GRID); ax.set_axisbelow(True)
for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
ax.legend(frameon=False, fontsize=13, loc='upper left')
ax.set_title('Where: resistive work by region, inhale vs exhale', loc='left', fontsize=15, color='#12305c')
# when: flow trace + region x time map
axq = fig.add_subplot(gs[0, 1]); axq.set_facecolor(SURF)
axq.plot(F[:, 0], F[:, 1] / a.rho * 6e4, color=INK2, lw=1.8); axq.axhline(0, color=GRID, lw=1)
axq.fill_between(times, -1e3, 1e3, where=inh, color=C_IN, alpha=0.10, lw=0); axq.fill_between(times, -1e3, 1e3, where=~inh, color=C_EX, alpha=0.10, lw=0)
qq = F[:, 1] / a.rho * 6e4; axq.set_ylim(qq.min() * 1.15, qq.max() * 1.15); axq.set_xlim(times[0] - dt / 2, times[-1] + dt / 2)
axq.set_ylabel('flow\n(L/min)', fontsize=12); axq.tick_params(labelbottom=False, labelsize=11)
axq.text(np.mean(np.array(times)[inh]), qq.min() * 0.6, 'INHALE', ha='center', color=C_IN, fontsize=12)
axq.text(np.mean(np.array(times)[~inh]), qq.max() * 0.6, 'EXHALE', ha='center', color=C_EX, fontsize=12)
for sp in ('top', 'right'): axq.spines[sp].set_visible(False)
axq.set_title('When: dissipated power by region over the breath', loc='left', fontsize=15, color='#12305c')
axm = fig.add_subplot(gs[1, 1], sharex=axq)
im = axm.imshow(np.clip(P, 0, None), aspect='auto', cmap='Blues', origin='upper', interpolation='nearest',
                extent=[times[0] - dt / 2, times[-1] + dt / 2, len(names) - 0.5, -0.5])
axm.set_yticks(range(len(names))); axm.set_yticklabels(names); axm.set_xlabel('time (s)')
cax = axm.inset_axes([1.012, 0.0, 0.022, 1.0]); cb = fig.colorbar(im, cax=cax); cb.set_label('dissipated power (mW)', fontsize=12)
fig.suptitle('%s: work of breathing by airway region' % a.subject, x=0.06, ha='left', y=0.97, fontsize=17, fontweight='bold', color='#12305c')
fig.savefig(os.path.join(a.out_dir, '%s_wob_by_region.png' % a.subject), dpi=130, facecolor=SURF)
print('wrote', os.path.join(a.out_dir, '%s_wob_by_region.png' % a.subject))
