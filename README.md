# Project BIO580 - Binary Detectability in MR Reconstruction

Magnetic Resonance Imaging (MRI) provides detailed anatomical images, but acquiring high-quality images can be time-consuming. Acceleration techniques that skip k-space lines are commonly used to reduce acquisition time, but this leads to image artifacts and loss of diagnostic information.

This project explores **binary detectability** in undersampled MRI reconstructions using **Hotelling Observer (HO)** and **Channelized Hotelling Observer (CHO)** models. The focus is to evaluate how acceleration, noise, and lesion characteristics affect the ability to detect small structures in reconstructed images—going beyond traditional metrics such as SSIM and MSE.

---

## 🔬 Problem Setup

We simulate elliptical lesions embedded in background patches to model signal-present and signal-absent conditions under varying noise and acceleration levels.

### Hypotheses

- **Signal Present (H₁):**
g₁ = H(fₛ + f_b) + n
- **Signal Absent (H₂):**
g₂ = H(f_b) + n

Where:
- `fₛ`: signal (elliptical lesion)  
- `f_b`: background  
- `n`: additive noise  
- `H`: reconstruction operator (model)  

---

## 🧠 Hotelling Observer

The Hotelling observer template is defined as:
W_hot = [0.5 * (K₁ + K₂)]⁻¹ · (E[g₁] - E[g₂])

- `K₁`, `K₂`: Covariance matrices under hypotheses H₁ and H₂  
- `E[g₁]`, `E[g₂]`: Mean vectors of signal-present and signal-absent images

The test statistic for a given patch `g` is:

t_hot = W_hotᵀ · g



This allows us to compute the **AUC** and evaluate detectability performance across varying conditions.

---

## 📦 Project Structure

| Module                      | Description                                                                 |
|----------------------------|-----------------------------------------------------------------------------|
| `Data Generation`          | Simulates signal-present/absent patches with varying noise and acceleration |
| `HotellingObserver`        | Implements standard Hotelling observer                                      |
| `ChannelizedObserver`      | Implements CHO using frequency and spatial filters                          |
| `MetricComparison`         | Computes and compares SSIM, MSE, AUC, and SNR                               |
| `Visualization`            | Plots ROC curves, confusion matrices, and test statistic distributions      |
| `AnalysisTools`            | Functions for covariance estimation, test statistic computation, etc.       |

---

## 📊 Metrics

While traditional metrics like:

- **MSE** – Average pixel-wise error  
- **SSIM** – Structural Similarity Index
- **AUC_HO** - AUC for hotelling observer
- **AUC_CHO** - AUC for hotelling observer
  
are included, this project emphasizes **task-based** metrics like AUC derived from observer models, which are more aligned with clinical goals.

---

