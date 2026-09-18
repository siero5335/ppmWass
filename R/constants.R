# Internal constants

# Typical neutral losses used by Mref estimation and optional loss-projection.
# NOTE: Multiple entries can share the same nominal mass; matching is performed by exact mass within ppm.
#       When multiple candidates fall within tolerance (rare), the closest exact mass is used (see project_to_typical_losses()).
# Organized into sets so derivatization-related losses can be toggled on/off.

TYPICAL_LOSSES_CORE <- list(
  list(nominal = 15, exact = 15.0235, formula = "CH3", class = "alkyl"),
  list(nominal = 17, exact = 17.0027, formula = "OH", class = "hydroxyl"),
  list(nominal = 17, exact = 17.0265, formula = "NH3", class = "amine"),
  list(nominal = 18, exact = 18.0106, formula = "H2O", class = "hydroxyl"),
  list(nominal = 28, exact = 27.9949, formula = "CO", class = "carbonyl"),
  list(nominal = 28, exact = 28.0313, formula = "C2H4", class = "hydrocarbon"),
  list(nominal = 29, exact = 29.0027, formula = "CHO", class = "aldehyde"),
  list(nominal = 29, exact = 29.0391, formula = "C2H5", class = "alkyl"),
  list(nominal = 31, exact = 31.0184, formula = "CH3O", class = "methoxy"),
  list(nominal = 31, exact = 31.0422, formula = "CH5N", class = "amine"),
  list(nominal = 32, exact = 32.0262, formula = "CH3OH", class = "methanol"),
  list(nominal = 42, exact = 42.0106, formula = "C2H2O", class = "ketene"),
  list(nominal = 42, exact = 42.0470, formula = "C3H6", class = "hydrocarbon"),
  list(nominal = 43, exact = 43.0184, formula = "CH3CO", class = "acetyl"),
  list(nominal = 43, exact = 43.0548, formula = "C3H7", class = "alkyl"),
  list(nominal = 44, exact = 43.9898, formula = "CO2", class = "carboxyl"),
  list(nominal = 44, exact = 44.0262, formula = "C2H4O", class = "acetaldehyde"),
  list(nominal = 45, exact = 44.9977, formula = "CHO2", class = "carboxyl"),
  list(nominal = 45, exact = 45.0340, formula = "C2H5O", class = "ethoxy"),
  list(nominal = 45, exact = 45.0215, formula = "CH3NO", class = "amide"),
  list(nominal = 46, exact = 46.0055, formula = "CH2O2", class = "formic_acid"),
  list(nominal = 46, exact = 46.0419, formula = "C2H6O", class = "ethanol"),
  list(nominal = 56, exact = 55.9898, formula = "C2O2", class = "dicarbonyl"),
  list(nominal = 56, exact = 56.0626, formula = "C4H8", class = "hydrocarbon"),
  list(nominal = 60, exact = 60.0211, formula = "C2H4O2", class = "acetic_acid")
)

TYPICAL_LOSSES_EXTENDED <- list(
  list(nominal = 16, exact = 16.0313, formula = "CH4", class = "hydrocarbon"),
  list(nominal = 20, exact = 20.0062, formula = "HF", class = "halogen_acid"),
  list(nominal = 26, exact = 26.0157, formula = "C2H2", class = "hydrocarbon"),
  list(nominal = 27, exact = 27.0109, formula = "HCN", class = "nitrile"),
  list(nominal = 28, exact = 28.0061, formula = "N2", class = "nitrogen"),
  list(nominal = 30, exact = 29.9980, formula = "NO", class = "nitroso"),
  list(nominal = 30, exact = 30.0106, formula = "CH2O", class = "aldehyde"),
  list(nominal = 30, exact = 30.0470, formula = "C2H6", class = "hydrocarbon"),
  list(nominal = 34, exact = 33.9877, formula = "H2S", class = "sulfide"),
  list(nominal = 36, exact = 35.9767, formula = "HCl", class = "halogen_acid"),
  list(nominal = 44, exact = 44.0626, formula = "C3H8", class = "hydrocarbon"),
  list(nominal = 46, exact = 45.9929, formula = "NO2", class = "nitro"),
  list(nominal = 64, exact = 63.9619, formula = "SO2", class = "sulfur_oxide"),
  list(nominal = 80, exact = 79.9568, formula = "SO3", class = "sulfur_oxide"),
  list(nominal = 80, exact = 79.9262, formula = "HBr", class = "halogen_acid"),
  list(nominal = 128, exact = 127.9123, formula = "HI", class = "halogen_acid")
)

TYPICAL_LOSSES_DERIV_TMS <- list(
  list(nominal = 72, exact = 72.0395, formula = "C3H8Si", class = "deriv_TMS"),
  list(nominal = 73, exact = 73.0474, formula = "C3H9Si", class = "deriv_TMS"),
  list(nominal = 89, exact = 89.0423, formula = "C3H9OSi", class = "deriv_TMS"),
  list(nominal = 90, exact = 90.0501, formula = "C3H10OSi", class = "deriv_TMS")
)

TYPICAL_LOSSES_DERIV_TBDMS <- list(
  list(nominal = 57, exact = 57.0704, formula = "C4H9", class = "deriv_TBDMS"),
  list(nominal = 115, exact = 115.0943, formula = "C6H15Si", class = "deriv_TBDMS"),
  list(nominal = 132, exact = 132.0970, formula = "C6H16OSi", class = "deriv_TBDMS")
)

# Backward-compatible alias (core only)
TYPICAL_LOSSES <- TYPICAL_LOSSES_CORE
META_COLS <- c(
  "Alignment ID", "Average Rt(min)", "Average RI", "Quant mass",
  "Metabolite name", "Fill %", "Reference RT", "Reference RI",
  "Formula", "Ontology", "INCHIKEY", "SMILES",
  "Annotation tag (VS1.0)", "RT/RI matched", "EI-MS matched",
  "Comment", "Manually modified for quantification",
  "Manually modified for annotation", "Total score",
  "RT similarity", "RI similarity", "Total spectrum similarity",
  "Dot product", "Reverse dot product", "Fragment presence %",
  "S/N average", "Spectrum reference file name", "EI spectrum"
)

INTERNAL_COLS <- c(
  "mean", "known", "compound_id", "blank_mean", "nonzero"
)
