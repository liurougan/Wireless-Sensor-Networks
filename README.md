Quantized Event-Triggered Resilient Filtering for WSNs
MATLAB simulation codes for the paper "Quantized Event-Triggered Resilient Filtering for Bandwidth-Limited Wireless Sensor Networks".

Overview
This repository contains all the scripts to reproduce the experiments and figures in the paper. The proposed scheme combines:

Innovation-based event-triggered communication

Uniform quantization

Lightweight outlier rejection
for communication-efficient and attack-resilient state estimation in bandwidth-limited wireless sensor networks.

Main Files
main_experiments.m – Main script: RMSE and communication comparison (Full KF, ET-KF, QET, Periodic QKF).

exp_attack_resilient_QET.m – Sensor bias attack experiment (standard QET vs. resilient QET).

exp_delta_and_noise.m – Sensitivity analysis w.r.t. quantization step ∆ and measurement noise R.

uniform_quantizer.m – Uniform mid‑rise quantizer function.

sweep_threshold.m – Function to compute the communication‑accuracy trade‑off curve.

system_block.m – Generates the system block diagram (Fig. 1 in the paper).

How to Run
Open MATLAB.

Run main_experiments.m to generate the main performance comparisons (Figs. 2‑4).

Run exp_attack_resilient_QET.m for the attack‑resilience results (Fig. 5).

Run exp_delta_and_noise.m for the ∆ and R sensitivity plots (Fig. 6).

All figures are automatically saved as .eps and .pdf in the working directory.

License
This code is provided for academic use. Please cite the corresponding paper if you use it in your research.
