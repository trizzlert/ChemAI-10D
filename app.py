import streamlit as st
import cv2
import numpy as np
import pytesseract
from chempy import balance_stoichiometry
from chempy.chemistry import Substance
from collections import defaultdict
from periodictable import elements  # NEW: for atomic number to symbol

# --- Image Preprocessing ---
def preprocess_image(image_bytes):
    nparr = np.frombuffer(image_bytes, np.uint8)
    img = cv2.imdecode(nparr, cv2.IMREAD_COLOR)
    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    blur = cv2.GaussianBlur(gray, (3, 3), 0)
    _, thresh = cv2.threshold(blur, 127, 255, cv2.THRESH_BINARY)
    return thresh

# --- OCR Processing ---
def extract_equation_from_image(image_bytes):
    preprocessed = preprocess_image(image_bytes)
    config = r'--oem 3 --psm 6'
    text = pytesseract.image_to_string(preprocessed, config=config)
    return text.strip().replace("=", "->")

# --- Equation Formatter ---
def format_equation(reac, prod):
    def fmt(compounds):
        return " + ".join(f"{'' if v == 1 else v}{k}" for k, v in compounds.items())
    return f"{fmt(reac)} → {fmt(prod)}"

# --- Reaction Type Detection ---
def identify_reaction_type(reactants, products):
    reac_len = len(reactants)
    prod_len = len(products)

    if reac_len == 1 and prod_len > 1:
        return "Decomposition"
    elif reac_len > 1 and prod_len == 1:
        return "Composition (Synthesis)"
    elif reac_len == 2 and prod_len == 2:
        return "Double Displacement"
    elif reac_len == 2 and prod_len == 1 or reac_len == 1 and prod_len == 2:
        return "Single Displacement"
    else:
        return "Other"

# --- Element Symbol Utility ---
def element_symbol(e):
    if isinstance(e, int):
        try:
            return elements[e].symbol
        except:
            return str(e)
    return str(e)

# --- Redox Detection ---
def oxidation_state_changes(reactants, products):
    def get_states(compounds):
        states = defaultdict(dict)
        for compound, coeff in compounds.items():
            try:
                s = Substance.from_formula(compound)
                ox_states = s.composition
                for elem, ox in ox_states.items():
                    states[compound][elem] = ox
            except Exception:
                continue
        return states

    reac_states = get_states(reactants)
    prod_states = get_states(products)

    redox_found = False
    oxidation = reduction = False
    changed_elements = []

    for r_cmpd, r_elems in reac_states.items():
        for p_cmpd, p_elems in prod_states.items():
            for elem in r_elems:
                if elem in p_elems:
                    if r_elems[elem] != p_elems[elem]:
                        redox_found = True
                        changed_elements.append(elem)
                        if r_elems[elem] < p_elems[elem]:
                            oxidation = True
                        elif r_elems[elem] > p_elems[elem]:
                            reduction = True

    changed_symbols = ', '.join(element_symbol(e) for e in set(changed_elements))

    if not redox_found:
        return "Not a redox reaction"

    if oxidation and reduction:
        return f"Redox Reaction (elements changed: {changed_symbols})"
    elif oxidation:
        return f"Oxidation Reaction (element: {changed_symbols})"
    elif reduction:
        return f"Reduction Reaction (element: {changed_symbols})"
    else:
        return "Unknown redox type"

# --- Chemical Balancing ---
def balance_equation(equation):
    try:
        if "->" not in equation:
            return "❌ Error: Equation must contain '->'", None, None

        lhs, rhs = equation.split("->")
        reactants = {r.strip() for r in lhs.strip().split("+")}
        products = {p.strip() for p in rhs.strip().split("+")}

        reac, prod = balance_stoichiometry(reactants, products)
        balanced_eq = format_equation(reac, prod)
        return balanced_eq, reac, prod

    except Exception as e:
        return f"❌ Error: {str(e)}", None, None

# --- Streamlit UI ---
st.set_page_config(page_title="Chemical Equation Balancer", page_icon="🧪")

st.title("🧪ChemAI ")
st.markdown("Balance chemical equations and identify reaction types, including redox classification. Made by Students of class 10-D")

# --- OCR Upload Option ---
st.subheader("📷 Upload an Image")
uploaded_image = st.file_uploader("Upload a photo of the equation (e.g. from textbook or whiteboard)", type=["png", "jpg", "jpeg"])

ocr_equation = ""
if uploaded_image:
    st.image(uploaded_image, caption="Uploaded Image", use_column_width=True)
    ocr_equation = extract_equation_from_image(uploaded_image.read())
    st.write(f"🔍 **Detected Equation:** `{ocr_equation}`")

# --- Manual Input Option ---
st.subheader("✍️ Or Enter Equation Manually")
manual_equation = st.text_input("")

# --- Combine and Process ---
equation_to_balance = manual_equation or ocr_equation

# ✅ Only one button
if st.button("⚖️ Balance Equation"):
    if equation_to_balance:
        st.info(f"Unbalanced: `{equation_to_balance}`")
        balanced_result, reac, prod = balance_equation(equation_to_balance)

        if reac and prod:
            st.success(f"✅ Balanced Equation: `{balanced_result}`")

            # Detect Reaction Type
            reaction_type = identify_reaction_type(reac, prod)
            st.info(f"🔬 Reaction Type: **{reaction_type}**")

            # Redox Classification
            redox_result = oxidation_state_changes(reac, prod)
            st.info(f"🧪 Redox Analysis: **{redox_result}**")
        else:
            st.error(balanced_result)
    else:
        st.warning("Please upload an image or enter an equation manually.")
