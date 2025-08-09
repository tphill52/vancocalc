import streamlit as st
import math

st.title("Vancomycin Bayesian AUC Calculator")

# User inputs
weight = st.number_input("Weight (kg)", value=75.0)
height = st.number_input("Height (in)", value=72.0)
age = st.number_input("Age (years)", value=80)
sex = st.selectbox("Sex", ["Male", "Female"])
scr = st.number_input("Serum Creatinine (mg/dL)", value=0.9)
trough = st.number_input("Observed Trough (mg/L)", value=14.5)
mic = st.number_input("MIC", value=1.0)
auc_low, auc_high = 400, 600

# Calculations
bmi = weight / ((height * 0.0254) ** 2)
if sex == "Male":
    crcl = ((140 - age) * weight) / (72 * scr)
else:
    crcl = 0.85 * ((140 - age) * weight) / (72 * scr)

vd = 0.7 * weight
cl = crcl * 0.06  # L/hr
ke = cl / vd
t_half = math.log(2) / ke
loading_dose = 25 * weight

# AUC targets
def daily_dose_for_auc(target_auc):
    return target_auc * cl

auc_targets = {auc: daily_dose_for_auc(auc) for auc in [auc_low, 500, auc_high]}

# Predicted trough at 12 hours (prior Bayesian)
predicted_trough_prior = (loading_dose / vd) * math.exp(-ke * 12)

# Bayesian update for trough (weighted average approach)
prior_sigma = 3.0   # uncertainty in predicted trough (population PK uncertainty)
obs_sigma = 1.5     # uncertainty in observed lab measurement
posterior_mu = (predicted_trough_prior / prior_sigma**2 + trough / obs_sigma**2) / (1/prior_sigma**2 + 1/obs_sigma**2)
posterior_sigma = 1.0 / math.sqrt(1/prior_sigma**2 + 1/obs_sigma**2)

# Output
st.subheader("Results")
st.write(f"**BMI:** {bmi:.2f} kg/m²")
st.write(f"**CrCl:** {crcl:.1f} mL/min")
st.write(f"**Vd:** {vd:.1f} L")
st.write(f"**Clearance:** {cl:.3f} L/hr")
st.write(f"**Half-life:** {t_half:.2f} h")
st.write(f"**Loading Dose:** {loading_dose:.0f} mg")

st.write(f"**Predicted trough (prior) at 12h:** {predicted_trough_prior:.2f} mg/L")
st.write(f"**Observed trough:** {trough} mg/L")
st.write(f"**Posterior trough estimate (Bayesian update):** {posterior_mu:.2f} ± {posterior_sigma:.2f} mg/L")

df = pd.DataFrame({
    "Target AUC": list(auc_targets.keys()),
    "Total Daily Dose (mg)": [round(v, 1) for v in auc_targets.values()],
    "q12h Dose (mg)": [round(v/2, 1) for v in auc_targets.values()]
})
st.table(df)

st.write(f"**Suggested dosing interval:** q12h (based on CrCl of {crcl:.1f} mL/min)")
