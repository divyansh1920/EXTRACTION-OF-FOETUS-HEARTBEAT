# 🫀 Extraction of Foetus Heartbeat

Fetal ECG extraction from maternal abdominal signals using Adaptive Noise Cancellation (ANC) — implemented in MATLAB for both SISO and MISO configurations.

## Problem
Fetal ECG buried inside maternal ECG. Mother's heartbeat (MECG) dominates — much larger amplitude. Goal: isolate fetal signal cleanly.

## Dataset
`foetal_ecg.dat` — 9 channels, 2500 samples @ 500 Hz
- Columns 2–6: Abdominal signals (Mother + Fetus)
- Columns 7–9: Thoracic signals (Mother only)

## Approach

**Signal flow:**
Thoracic (reference noise) → Adaptive Filter → subtract from Abdominal (primary) → Fetal ECG

**Algorithms used (filter order p+1 = 12):**
| Algorithm | Step Size | Key Feature |
|-----------|-----------|-------------|
| LMS | µ = 2×10⁻⁸ | Fixed step, simple |
| NLMS | β = 0.001–0.005 | Normalized, input-power adaptive |
| LLMS | µ = 2×10⁻⁷, γ = 0.001 | Leaky — adds stability |

## Architectures

**SISO** — average all thoracic → single filter → subtract from average abdominal

**MISO** — each thoracic signal → separate filter → average outputs → subtract from average abdominal

## Results
- NLMS: best noise suppression, smallest residual maternal peaks
- LLMS: fastest convergence (LLMS > LMS > NLMS)
- MISO: more filtering stages → cleaner output than SISO

## Files
├── new_siso.m       # SISO implementation (LMS, NLMS, LLMS)

├── new_miso.m       # MISO implementation (LMS, NLMS, LLMS)

├── convm.m          # Convolution matrix utility

└── foetal_ecg.dat   # Dataset

## How to Run
```matlab
% SISO
run('new_siso.m')

% MISO
run('new_miso.m')
```

## Tech
MATLAB · Adaptive Filtering · Digital Signal Processing · Biomedical Engineering

## Authors
Ankush Badgujar & Jigyas Arora — VIT Chennai (ECE2006, Nov 2020)
