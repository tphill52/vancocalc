
import streamlit as st
import math
import pandas as pd  # Importing Pandas for DataFrame functionality
import matplotlib.pyplot as plt  # Importing the Matplotlib library for plotting

# User Inputs
st.title("Vancomycin Bayesian AUC Calculator")

# Inputs for patient characteristics
weight = st.number_input("Weight (kg)", value=75.0)
height = st.number_input("Height (in)", value=72.0)
age = st.number_input("Age (years)", value=80)
sex = st.selectbox("Sex", ["Male", "Female"])
scr = st.number_input("Serum Creatinine (mg/dL)", value=0.9)
trough = st.number_input("Observed Trough (mg/L)", value=14.5)
mic = st.number_input("MIC", value=1.0)

# Additional inputs for patient classification
bmi = weight / ((height * 0.0254) ** 2)
is_icu = st.checkbox("Is this patient in ICU?")
extreme_obesity = weight >= 120 or bmi >= 40  # BMI > 40 or weight >= 120kg

# Define Bayesian models and CL and Vd calculations based on patient profile
def select_bayesian_model():
    global cl_vanco, vd
    if extreme_obesity:
        # Adane et al., 2015 (Extreme obesity model)
        st.write("Using Extreme Obesity Model (Adane et al., 2015)")
        # CL calculation for extreme obesity
        cl_vanco = (6.54 * scr * weight) / 125
        # Vd calculation for extreme obesity
        vd = 0.51 * weight
        return 3.0, 1.5  # Uncertainty values for the model
    elif is_icu and bmi < 30:
        # Roberts et al., 2011 (ICU model)
        st.write("Using ICU Model (Roberts et al., 2011)")
        # CL calculation for ICU patients
        cl_vanco = (4.58 * scr * 1.73) / 100
        # Vd calculation for ICU patients
        vd = 1.53 * weight
        return 2.5, 1.2  # Uncertainty values for the ICU model
    elif is_icu and bmi >= 30:
        # Masich et al., 2020 (Obesity + ICU model)
        st.write("Using ICU + Obesity Model (Masich et al., 2020)")
        # CL calculation for critically ill with obesity
        cl_vanco = 3.23 * (scr / 40) ** 0.69
        # Vd calculation for ICU + Obesity patients
        vd = 0.78 * weight
        return 2.7, 1.3  # Uncertainty values for obesity + ICU model
    else:
        # Buelga et al., 2005 (Default model for general hospitalized patients)
        st.write("Using Default Model (Buelga et al., 2005)")
        # CL calculation for general hospitalized patients
        cl_vanco = (scr * 60) / 1000 * 1.08
        # Vd calculation for general hospitalized patients
        vd = 0.98 * weight
        return 3.0, 1.5  # Default model uncertainty values

# Select the appropriate Bayesian model
prior_sigma, obs_sigma = select_bayesian_model()

# Predicted trough calculation (prior)
predicted_trough = (25 * weight / vd) * math.exp(-cl_vanco / vd * 12)

# Bayesian update for posterior mean (mu) and posterior sigma
posterior_mu = (predicted_trough / (prior_sigma ** 2) + trough / (obs_sigma ** 2)) / (1 / (prior_sigma ** 2) + 1 / (obs_sigma ** 2))
posterior_sigma = 1.0 / math.sqrt(1 / (prior_sigma ** 2) + 1 / (obs_sigma ** 2))

# Output
st.subheader("Results")
st.write(f"**BMI:** {bmi:.2f} kg/m²")
st.write(f"**CrCl:** {scr:.1f} mL/min")
st.write(f"**Vd:** {vd:.1f} L")
st.write(f"**Clearance (CL):** {cl_vanco:.3f} L/hr")
st.write(f"**Half-life:** {math.log(2)/cl_vanco:.2f} h")
st.write(f"**Loading Dose:** {25 * weight:.0f} mg")

st.write(f"**Predicted trough (prior) at 12h:** {predicted_trough:.2f} mg/L")
st.write(f"**Observed trough:** {trough} mg/L")
st.write(f"**Posterior trough estimate (Bayesian update):** {posterior_mu:.2f} ± {posterior_sigma:.2f} mg/L")

# Plot Vancomycin Concentration over Time
time_hours = range(0, 25)  # Time in hours (0 to 24)
concentration = [(25 * weight / vd) * math.exp(-cl_vanco / vd * t) for t in time_hours]

# Plotting the graph
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(time_hours, concentration, label="Predicted Vancomycin Concentration", color='b', linewidth=2)
ax.set_xlabel('Time (hours)', fontsize=12)
ax.set_ylabel('Vancomycin Concentration (mg/L)', fontsize=12)
ax.set_title('Predicted Vancomycin Concentration-Time Curve', fontsize=14)
ax.grid(True)
ax.legend()

# Display the plot in Streamlit
st.pyplot(fig)

# Optionally display the recommended AUC-based doses
auc_low, auc_high = 400, 600
def daily_dose_for_auc(target_auc):
    return target_auc * cl_vanco

auc_targets = {auc: daily_dose_for_auc(auc) for auc in [auc_low, 500, auc_high]}
df = pd.DataFrame({
    "Target AUC": list(auc_targets.keys()),
    "Total Daily Dose (mg)": [round(v, 1) for v in auc_targets.values()],
    "q12h Dose (mg)": [round(v/2, 1) for v in auc_targets.values()]
})
st.table(df)

st.write(f"**Suggested dosing interval:** q12h (based on CrCl of {scr:.1f} mL/min)")
