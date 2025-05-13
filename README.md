# Versatile Reservoir Computing for Heterogeneous Complex Networks

This project implements the methods described in the paper:  
**[Versatile Reservoir Computing for Heterogeneous Complex Networks](https://arxiv.org/abs/xxxxxxx)**.  
We propose a novel machine learning scheme termed **versatile reservoir computing ** that enables a single small-scale reservoir computer, trained on partial time series data, to sustain the collective dynamics of large-scale heterogeneous complex networks. The trained machine can accurately replicate the dynamics of any element in the network and preserve the network's collective behavior when substituting failed components over finite time horizons. Effectiveness is validated on three representative network models.

---

## Features  
- **Heterogeneous Network Adaptation**: Handle networks with diverse node parameters and connection topologies.  
- **Universal Substitution**: A single trained machine substitutes *any* network element while preserving collective dynamics.  
- **Multi-Model Support**: Validated on:  
  - Homogeneous phase oscillators  
  - Heterogeneous phase oscillators  
  - Chaotic Lorenz oscillators  
- **Efficient Training**: Optimized reservoir hyperparameters for robust performance across dynamical regimes.  

---

## Project Structure  
```bash
.
├── Model_A/                         # Homogeneous phase oscillators
│   ├── train/                       # Training scripts
│   │   └── opt_Kur_para.m           # Hyperparameter optimization
│   └── substitute_oscillator/       
│       └── R_predict.m              # Substitution test (Fig. 2)
│
├── Model_B/                         # Heterogeneous phase oscillators
│   ├── train/                       
│   │   └── opt_Kur_para.m           # Shared training framework
│   └── substitute_oscillator/       
│       └── R_predict_delay.m        # Delayed substitution (Fig. 3)
│
└── Model_C/                         # Chaotic Lorenz oscillators
    ├── train/
    │   └── opt_chua_coupled_classic.m  # Chaotic dynamics training
    └── substitute_oscillator/
        └── R_predict_chua_delay_2.m    # Chaotic substitution (Fig. 4)
