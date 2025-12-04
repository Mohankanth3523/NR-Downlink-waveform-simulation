# NR Downlink Waveform Simulation and Performance Evaluation

This repository contains a MATLAB simulation of **5G NR downlink waveform generation** and evaluation of key performance metrics.  
The project is implemented without relying on the 5G Toolbox, using custom OFDM/QAM functions.

---

## 🚀 Features
- Supports multiple numerologies (15, 30, 60, 120 kHz subcarrier spacing)
- Resource grid mapping for QPSK and QAM modulation
- OFDM modulation with cyclic prefix
- Multipath fading channel with AWGN
- Performance evaluation of:
  - Spectrum
  - Power Delay Profile (PDP)
  - Bit Error Rate (BER) vs SNR
  - Error Vector Magnitude (EVM) vs SNR
  - Throughput vs SNR

---

## ⚙️ Requirements
- MATLAB R202x (tested on R2023b)
- No 5G Toolbox required

---

## ▶️ How to Run
1. Clone the repository:
   ```bash
   git clone https://github.com/<your-username>/nr-downlink-waveform-simulation.git
