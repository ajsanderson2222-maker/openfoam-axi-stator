#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import re

# ------------------------------------------------------------------
# 1. Total pressure drop
# ------------------------------------------------------------------
pt_in  = np.loadtxt("postProcessing/pt_inlet/0/surfaceFieldValue.dat")
pt_out = np.loadtxt("postProcessing/pt_outlet/0/surfaceFieldValue.dat")

t       = pt_in[:, 0]
dpt     = pt_in[:, 1] - pt_out[:, 1]
dpt_psi = dpt / 6894.76   # Pa → psi

fig, axes = plt.subplots(1, 2, figsize=(13, 5))

ax = axes[0]
ax.plot(t * 1e3, dpt_psi, lw=0.8, color='steelblue')
ax.set_xlabel("Time (ms)")
ax.set_ylabel("ΔPt (psi)")
ax.set_title("Total Pressure Drop  Pt_in − Pt_out")
ax.grid(True, alpha=0.3)
ax.axhline(dpt_psi[-1], color='red', ls='--', lw=0.8,
           label=f"Final: {dpt_psi[-1]:.2f} psi")
ax.legend()

# ------------------------------------------------------------------
# 2. Residuals from log
# ------------------------------------------------------------------
pattern = re.compile(
    r"^Time = (?P<t>[0-9eE.+\-]+)s$|"
    r"smoothSolver:.*?Solving for (?P<field>\w+),.*?Initial residual = (?P<res>[0-9eE.+\-]+)",
    re.MULTILINE
)

times, fields_data = [], {}
cur_t = None

with open("log.shockFluid") as fh:
    for line in fh:
        line = line.strip()
        m_t = re.match(r"^Time = ([0-9eE.+\-]+)s$", line)
        if m_t:
            cur_t = float(m_t.group(1))
            times.append(cur_t)
            continue
        m_r = re.match(r"smoothSolver:.*?Solving for (\w+),.*?Initial residual = ([0-9eE.+\-]+)", line)
        if m_r and cur_t is not None:
            field = m_r.group(1)
            res   = float(m_r.group(2))
            if field not in fields_data:
                fields_data[field] = ([], [])
            fields_data[field][0].append(cur_t)
            fields_data[field][1].append(res)

ax = axes[1]
colors = plt.cm.tab10.colors
for i, (field, (ft, fr)) in enumerate(fields_data.items()):
    if max(fr) > 0:
        ax.semilogy(np.array(ft) * 1e3, fr, lw=0.7,
                    color=colors[i % 10], label=field)

ax.set_xlabel("Time (ms)")
ax.set_ylabel("Initial residual")
ax.set_title("Solver Residuals")
ax.legend(fontsize=8)
ax.grid(True, which='both', alpha=0.3)

plt.tight_layout()
plt.savefig("results.png", dpi=150)
print("Saved results.png")
print(f"\nFinal ΔPt = {dpt_psi[-1]:.3f} psi  ({dpt[-1]:.0f} Pa)")
print(f"Run end time = {t[-1]*1e3:.3f} ms,  {len(t)} samples")
