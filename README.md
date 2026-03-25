# Vapour-Compression-Refrigeration-System-3E-Analysis
3E (Energy, Exergy, Economic) performance analysis and optimization of a Vapour Compression Refrigeration System (VCRS) using R1234yf, R134a, and R744 refrigerants.

# Project Overview
This project provides a **3E (Energy, Exergy, Economic)** analysis of a Vapour Compression Refrigeration System. 
The study evaluates **R1234yf, R134a and R744** as zero-ODP alternatives to traditional refrigerants using a computational model.

**Energy:** Evaluates the **Coefficient of Performance (COP)** by analyzing refrigerating effect versus compressor.

**Exergy:** Quantifies **Exergy Destruction ($E_d$)** and **Exergetic Efficiency ($\eta_{ex}$)** to identify system irreversibilities.

**Economic:** Assesses the cost-viability of substituting refrigerants based on thermodynamic performance.

## How to Run the Model
1.  **Clone the Repository:**
    ```bash
    git clone https://github.com/Akshat-Singh-Kshatriya/VCRS-4E-Analysis.git
    cd VCRS-4E-Analysis
    ```
2.  **Install Dependencies:**
    ```bash
    pip install CoolProp
    ```
    ```bash
    pyenv
    ```
3.  **Run the Analysis:**
    ```bash
    vcrs
    ```
    ```bash
    vcrs_with_ejector_pin
    ```
     ```bash
     ejector_pin_with_plots
    ```
      ```bash
    combined_plots_varying_T_evap
    ```
    ```bash
    combined_plots_varying_T_cond
    ```


