# === BLOC 1 : CONFIGURATION, SHIFT VERTICAL ET CALCULS ===
import streamlit as st
import numpy as np
import plotly.graph_objects as go

st.set_page_config(page_title="Cinématique Relevage", layout="wide")
st.title("📐 Schéma Cinématique : Théorème des Forces Concourantes")
st.markdown("Démonstration graphique de la statique : Ligne d'action des forces avec origine Z calée au niveau du sol.")

angle_bras_deg = st.sidebar.slider("🕹️ Position : Angle du bras supérieur (1) [°]", -15.0, 45.0, 5.0, step=1.0)
L_3eme_point = st.sidebar.slider("🔗 Longueur du Troisième Point (2) [m]", 0.65, 0.95, 0.82, step=0.01)
charge_kg = st.sidebar.slider("🏋️ Masse de la charrue (au CG) [kg]", 100, 3000, 1500, step=100)

# --- VARIABLES DE RÉGLAGE DES DÉCALAGES (SHIFTS) ---
decalage_relevage_X = -0.75
shift_sol_Z = 0.85  # Décalage pour amener l'interface pneu/sol à l'origine Z = 0

# --- CONFIGURATION DES PIVOTS FIXES TRANSLATÉS (X et Z) ---
O6 = np.array([0.95 + decalage_relevage_X, -0.42 + shift_sol_Z])
O4 = np.array([1.00 + decalage_relevage_X, -0.15 + shift_sol_Z])
O1 = np.array([1.02 + decalage_relevage_X, 0.18 + shift_sol_Z])
O2 = np.array([0.96 + decalage_relevage_X, 0.32 + shift_sol_Z])

L_O1_B, L_O1_A, L_chandelle, L_O6_D, L_O6_C = 0.35, 0.22, 0.50, 0.85, 0.42
H_attelage_charrue, Dist_CG_X, Dist_CG_Z = 0.45, 0.60, 0.05

def intersect_cercles(O, r_O, B, r_B, bas=False):
    d = np.linalg.norm(B - O)
    if d > r_O + r_B or d < abs(r_O - r_B) or d == 0: return None
    a = (r_O**2 - r_B**2 + d**2) / (2 * d)
    h = np.sqrt(max(0.0, r_O**2 - a**2))
    p2 = O + a * (B - O) / d
    u_perp = np.array([-(B[1] - O[1]), (B[0] - O[0])]) / d
    sol1, sol2 = p2 + h * u_perp, p2 - h * u_perp
    if bas: return sol1 if sol1[1] < sol2[1] else sol2
    return sol1 if sol1[0] > sol2[0] else sol2

angle_rad = np.radians(angle_bras_deg)
u_bras = np.array([np.cos(angle_rad), np.sin(angle_rad)])
A, B = O1 + L_O1_A * u_bras, O1 + L_O1_B * u_bras
L_verin_actuelle = np.hypot(A[0] - O4[0], A[1] - O4[1])

C = intersect_cercles(O6, L_O6_C, B, L_chandelle, bas=True)
if C is None: st.stop()

alpha_inf = np.atan2(C[1] - O6[1], C[0] - O6[0])
u_inf = np.array([np.cos(alpha_inf), np.sin(alpha_inf)])
D = O6 + L_O6_D * u_inf
E = intersect_cercles(O2, L_3eme_point, D, H_attelage_charrue)
if E is None: st.stop()

u_tete = (E - D) / np.linalg.norm(E - D)
u_bati = np.array([u_tete[1], -u_tete[0]])
G = D + Dist_CG_X * u_bati + Dist_CG_Z * u_tete
P_fin_bati = D + 1.3 * u_bati

F_gravity = charge_kg * 9.81
M_poids_D = - F_gravity * (G[0] - D[0])
u_3eme = (O2 - E) / np.linalg.norm(O2 - E)
levier_3eme_D = (E[0] - D[0]) * u_3eme[1] - (E[1] - D[1]) * u_3eme[0]
F_3eme_norme = - M_poids_D / levier_3eme_D if levier_3eme_D != 0 else 0
F_3eme_sur_charrue = F_3eme_norme * u_3eme
R_D = F_3eme_sur_charrue - np.array([0.0, F_gravity])

M_RD_O6 = (D[0] - O6[0]) * R_D[1] - (D[1] - O6[1]) * R_D[0]
u_chandelle = (B - C) / np.linalg.norm(B - C)
levier_chandelle_O6 = (C[0] - O6[0]) * u_chandelle[1] - (C[1] - O6[1]) * u_chandelle[0]
F_chandelle_norme = - M_RD_O6 / levier_chandelle_O6 if levier_chandelle_O6 != 0 else 0
R_C_chandelle = F_chandelle_norme * u_chandelle
R_O6 = R_C_chandelle + R_D

def intersect_droites(P1, u1, P2, u2):
    A_mat = np.array([[u1[0], -u2[0]], [u1[1], -u2[1]]])
    try: return P1 + np.linalg.solve(A_mat, P2 - P1)[0] * u1
    except np.linalg.LinAlgError: return None

I_charrue = intersect_droites(G, np.array([0.0, -1.0]), E, u_3eme)
I_bras = intersect_droites(D, R_D/np.linalg.norm(R_D), C, u_chandelle)
# === BLOC 2 : TRACÉ GRAPHIQUE AVEC ESSIEUX DÉCALÉS ET SANS GRILLE ===
col1, col2, col3 = st.columns(3)
col1.metric("Longueur du Vérin (4)", f"{L_verin_actuelle:.3f} m")
col2.metric("🏋️ Effort Vérin (Poussée)", f"{abs(F_chandelle_norme * 1.25)/1000:.1f} kN")
col3.metric("🔗 Effort Traction 3ème Point", f"{F_3eme_norme/1000:.1f} kN")

fig = go.Figure()

# Le bloc de sol s'étend désormais de Z = 0 (Nouveau niveau du sol) jusqu'à Z = -0.15m (bas du graphique)
fig.add_trace(go.Scatter(x=[-4.1, 4.1, 4.1, -4.1, -4.1], y=[0.0, 0.0, -1, -0.8, 0.0], mode='lines', fill='toself', fillcolor='rgba(139, 69, 19, 0.25)', line=dict(color='rgba(0,0,0,0)'), showlegend=False))
fig.add_trace(go.Scatter(x=[-4.2, 4.2], y=[0.0, 0.0], mode='lines', line=dict(color='gray', width=2), showlegend=False))

# Décalage des centres des roues en Z (+shift_sol_Z)
th = np.linspace(0, 2*np.pi, 60)
fig.add_trace(go.Scatter(x=0.85*np.cos(th), y=shift_sol_Z + 0.85*np.sin(th), mode='lines', line=dict(color='lightgray', width=2), showlegend=False))
fig.add_trace(go.Scatter(x=-2.35 + 0.60*np.cos(th), y=(shift_sol_Z - 0.25) + 0.60*np.sin(th), mode='lines', line=dict(color='lightgray', width=2), showlegend=False))

# Décalage du châssis, cabine et toit en Z
tracteur_z_shift = np.array([0.4, 0.6, 0.6, -0.3, -0.3]) + shift_sol_Z
fig.add_trace(go.Scatter(x=[-2.3, -2.25, -0.6, -1.75, -0.81], y=tracteur_z_shift, mode='markers+lines', line=dict(color='lightgray', width=2), showlegend=False))

cab_x = np.array([-0.60, 0.60, 0.50, -0.50, -0.60]) - 0.6
cab_z = np.array([0.60, 0.60, 1.45, 1.45, 0.60]) + shift_sol_Z
fig.add_trace(go.Scatter(x=cab_x, y=cab_z, mode='lines', line=dict(color='lightgray', width=2), showlegend=False))

toit_x = np.array([-0.55, 0.65, 0.60, -0.50, -0.55]) - 0.65
toit_z = np.array([1.45, 1.45, 1.52, 1.52, 1.45]) + shift_sol_Z
fig.add_trace(go.Scatter(x=toit_x, y=toit_z, mode='lines', line=dict(color='lightgray', width=2), showlegend=False))

def l(p1, p2, couleur='darkgray', w=4):
    fig.add_trace(go.Scatter(x=[p1[0], p2[0]], y=[p1[1], p2[1]], mode='lines+markers', line=dict(color=couleur, width=w), marker=dict(size=5), showlegend=False))

l(O1, B); l(O2, E); l(O4, A, 'lightgray', 6); l(B, C); l(O6, D); l(D, E)
fig.add_trace(go.Scatter(x=[D[0], P_fin_bati[0]], y=[D[1], P_fin_bati[1]], mode='lines', line=dict(color='black', width=6), showlegend=False))
fig.add_trace(go.Scatter(x=[G[0]], y=[G[1]], mode='markers', marker=dict(size=8, color='black', symbol='x'), showlegend=False))

echelle = 0.25 / 15000
def f(orig, vec, inter):
    fin = orig + vec * echelle
    fig.add_trace(go.Scatter(x=[orig[0], fin[0]], y=[orig[1], fin[1]], mode='lines+markers', line=dict(color='red', width=4), marker=dict(symbol='arrow-up', size=8, angleref='previous'), showlegend=False))
    if inter is not None: fig.add_trace(go.Scatter(x=[orig[0], inter[0]], y=[orig[1], inter[1]], mode='lines', line=dict(color='red', width=1.5, dash='dashdot'), showlegend=False))

f(D, R_D, I_charrue); f(O6, R_O6, I_bras); f(C, R_C_chandelle, I_bras); f(E, F_3eme_sur_charrue, I_charrue)
if I_bras is not None: fig.add_trace(go.Scatter(x=[D[0], I_bras[0]], y=[D[1], I_bras[1]], mode='lines', line=dict(color='red', width=1.5, dash='dashdot'), showlegend=False))
if I_charrue is not None: fig.add_trace(go.Scatter(x=[G[0], I_charrue[0]], y=[G[1], I_charrue[1]], mode='lines', line=dict(color='red', width=1.5, dash='dashdot'), showlegend=False))
fig.add_trace(go.Scatter(x=[G[0], G[0]], y=[G[1], G[1] - 0.20], mode='lines+markers', line=dict(color='red', width=4), marker=dict(symbol='arrow-down', size=8), showlegend=False))

if I_charrue is not None: fig.add_trace(go.Scatter(x=[I_charrue[0]], y=[I_charrue[1]], mode='markers+text', marker=dict(size=10, color='red', symbol='circle-open'), text=["I"], textposition="top right", showlegend=False))
if I_bras is not None: fig.add_trace(go.Scatter(x=[I_bras[0]], y=[I_bras[1]], mode='markers+text', marker=dict(size=10, color='red', symbol='circle-open'), text=["I'"], textposition="bottom left", showlegend=False))

corps_x = [0.35, 0.85]
for cx in corps_x:
    p_base = D + cx * u_bati
    p_soc_bas = p_base - 0.22 * u_tete + 0.12 * u_inf
    fig.add_trace(go.Scatter(x=[p_base[0], p_soc_bas[0]], y=[p_base[1], p_soc_bas[1]], mode='lines', line=dict(color="black", width=3), showlegend=False))
    p_soc_haut = p_base + 0.32 * u_tete + 0.12 * u_inf
    fig.add_trace(go.Scatter(x=[p_base[0], p_soc_haut[0]], y=[p_base[1], p_soc_haut[1]], mode='lines', line=dict(color="black", width=1.5), showlegend=False))

# --- RECONFIGURATION AXES : VERROUILLAGE X [-4, +4] ET ZÉRO GRID ---
fig.update_layout(
    xaxis=dict(title="Longueur X (m)", range=[-4.0, 4.0], showgrid=False, gridcolor='rgba(0,0,0,0)'),
    yaxis=dict(title="Hauteur Z (m)", range=[-1, 2.8], showgrid=False, gridcolor='rgba(0,0,0,0)', scaleanchor="x", scaleratio=1),
    plot_bgcolor='white', height=650, showlegend=False
)
st.plotly_chart(fig, use_container_width=True)
