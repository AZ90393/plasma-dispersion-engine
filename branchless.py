
# ============================================================================
#COPYRIGHT & ALGORITHM ATTRIBUTION::
# © 2026 RMIO ENTERPRISE, LLC. All rights reserved.
# Derived from the Dawson Integral  in the Cephes Mathematical Library.
# Copyright 1984, 1987, 2000 by Stephen L. Moshier.
# Permission to use, copy, modify, and distribute this software for any
# purpose without fee is hereby granted, provided that this entire notice 
# is included in all copies of any software which is or includes a copy 
# or modification of this software.
# Optimized into a branchless NaN-free version  CMOV configuration using IEEE 754 hardware operations.
# ============================================================================"""
# ============================================================================ #
# ALGORITHM ATTRIBUTION & LINEAGE (PYTHON MODULE) 
#This module implements the rational approximations of the Dawson Integral 
# derived from the Cephes Mathematical Library backend.It replicates the
# high-precision numerical framework utilized by scipy.special.dawsn 
# without requiring heavy external dependencies.
 # ============================================================================#
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from scipy import special
# [Insert your dawsn_branchless, cephes_polevl, and cephes_p1evl functions here]


def cephes_polevl(x: float, coef: list) -> float:
    #Evaluates a standard polynomial using Horner's Method.
    ans = coef[0]
    for c in coef[1:]:
        ans = ans * x + c
    return ans

def cephes_p1evl(x: float, coef: list) -> float:
    #Evaluates a polynomial where the leading coefficient is implicitly 1.0.
    ans = x + coef[0]
    for c in coef[1:]:
        ans = ans * x + c
    return ans

# ============================================================================
# BRANCHLESS HIGH-PERFORMANCE DAWSON INTEGRAL FUNCTION
# ============================================================================

# Apply global high-res settings
#plt.rcParams['figure.dpi'] = 300
#plt.rcParams['savefig.dpi'] = 300

def dawsn_branchless(x: float,n:float) -> float: #def f(x: float, n: float) -> float:
    PI = 4 * np.arctan(1)
   # x = np.array(x, dtype=np.float128)
    e1_x=(np.exp(800*(x-3.25)))**n#SLIDER
    #e1_x=(np.exp(2*800*(x-3.25)))
    e2_x=(np.exp(800*(x+3.25)))**n#SLIDER
    #e2_x=(np.exp(2*800*(x+3.25)))
    e3_x=(np.exp(800*(x-6.25)))**n#SLIDER
    #e3_x=(np.exp(2*800*(x-6.25)))
    e4_x=(np.exp(800*(x+6.25)))**n#SLIDER
    #e4_x=(np.exp(2*800*(x+6.25)))
    R5_x=np.exp(-e1_x)*(-1*np.exp(-e2_x)+1)
    R6_x=((-1*np.exp(-e1_x)+1)+(np.exp(-e2_x)))*(np.exp(-e3_x))+((-1*np.exp(-e1_x)+1)+(-1*np.exp(-e2_x)))*(np.exp(-e4_x))
    R7_x=(-1*np.exp(-e3_x)+1)+(np.exp(-e4_x))
    R8_x=(np.exp(-e3_x))*(-1*np.exp(-e4_x)+1)
    R9_x=(2*np.arctan((1e300)*x))/(PI)
    #abs_x = abs(x)
    
    # 1. Safe domain selection to isolate polynomial inputs and prevent division-by-zero
   # x_A = abs_x if abs_x <= 3.25 else 1.0
    #x_BC = abs_x if abs_x > 3.25 else 5.0 # Avoids 0.0 down the road

   # t_A = x * x
    #t_BC = 1.0 / (x_BC * x_BC)

    # ========================================================================
    # INTERVAL A: 0 <= |x| <= 3.25 (Uses standard polevl for the denominator)
    # ========================================================================
    AN=[ 
     1.13681498971755972054E-11,
     8.49262267667473811108E-10,
     1.94434204175553054283E-8,
     9.53151741254484363489E-7,
     3.07828309874913200438E-6,
     3.52513368520288738649E-4,
     -8.50149846724410912031E-4,
     4.22618223005546594270E-2,
     -9.17480371773452345351E-2,
     9.99999999999999994612E-1,
     ]

    AD = [
            2.40372073066762605484E-11,
             1.48864681368493396752E-9,
             5.21265281010541664570E-8,
             1.27258478273186970203E-6,
             2.32490249820789513991E-5,
             3.25524741826057911661E-4,
             3.48805814657162590916E-3,
             2.79448531198828973716E-2,
             1.58874241960120565368E-1,
             5.74918629489320327824E-1,
             1.00000000000000000539E0,
    ]

    # FIX MAINTAINED: Both are standard cephes_polevl here to honor trailing 1.0
   
    t_A = abs(x)**(2*R5_x)
    F3_X =R5_x*(abs(x)**R5_x)*(cephes_polevl(t_A, AN) / cephes_polevl(t_A, AD))

    # ========================================================================
    # INTERVAL B: 3.25 < |x| <= 6.25 (Uses p1evl for the denominator)
    # ========================================================================
    BN = [
            5.08955156417900903354E-1,
           -2.44754418142697847934E-1,
            9.41512335303534411857E-2,
            -2.18711255142039025206E-2,
            3.66207612329569181322E-3,
            -4.23209114460388756528E-4,
            3.59641304793896631888E-5,
            -2.14640351719968974225E-6,
            9.10010780076391431042E-8,
            -2.40274520828250956942E-9,
            3.59233385440928410398E-11, 
    ]

    BD = [
           -6.31839869873368190192E-1,
            2.36706788228248691528E-1,
            -5.31806367003223277662E-2,
            8.48041718586295374409E-3,
            -9.47996768486665330168E-4,
            7.81025592944552338085E-5,
            -4.55875153252442634831E-6,
            1.89100358111421846170E-7,
            -4.91324691331920606875E-9,
            7.18466403235734541950E-11,
    ]
    #double t_BC = 1.0 / (std::pow(std::abs(x)+R5_x,2));
    t_BC = (1.0) / ((abs(x)+R5_x)**(2))
    #double F5_X=(1.0 / (2.0 *std::abs(x)+R5_x)) + (t_BC / (2.0 *std::abs(x)+R5_x)) * 
                   #(cephes_polevl(t_BC, BN, 10) / cephes_p1evl(t_BC, BD, 10));
    F5_X = (1.0) / (2.0 * (abs(x)+R5_x)) + ((t_BC) / (2.0 * (abs(x)+R5_x)) )* \
            (cephes_polevl(t_BC, BN) / cephes_p1evl(t_BC, BD))

    # ========================================================================
    # INTERVAL C: |x| > 6.25 (Uses p1evl for the denominator)
    # ========================================================================
    CN = [
            -5.90592860534773254987E-1,
            6.29235242724368800674E-1,
            -1.72858975380388136411E-1,
            1.64837047825189632310E-2,
            -4.86827613020462700845E-4,
    ]

    CD = [
           -2.69820057197544900361E0,
            1.73270799045947845857E0,
            -3.93708582281939493482E-1,
            3.44278924041233391079E-2,
            -9.73655226040941223894E-4,
    ]
    t1_BC = (1.0) / ((abs(x)+R8_x)**(2))
    F7_X = (1.0) / (2.0 * (abs(x)+R8_x)) +  ((t1_BC) / (2.0 * (abs(x)+R8_x)) ) * \
            (cephes_polevl(t1_BC, CN) / cephes_p1evl(t1_BC, CD))

    
    final_res_x=(R9_x)*(F3_X+(R6_x)*F5_X+(R7_x)*F7_X)
    
    return final_res_x  


def sum_error(n: float)->float:
 diff1=abs(special.dawsn(-3.25-.001)-dawsn_branchless(-3.25-.001,n))
 diff2=abs(special.dawsn(-3.25+.001)-dawsn_branchless(-3.25+.001,n))
 diff3=abs(special.dawsn(3.25-.001)-dawsn_branchless(3.25-.001,n))
 diff4=abs(special.dawsn(3.25+.001)-dawsn_branchless(3.25+.001,n))
 sum=(diff1+diff2+diff3+diff4)/4
 return sum

# 1. Generate smooth distribution samples across the domains
num=float('-inf')
# Force vector-based text rendering
plt.rcParams['svg.fonttype'] = 'none'

# STEP 1: Use a perfectly square figure canvas size (e.g., 8x8 inches)
#fig, ax = plt.subplots(figsize=(8, 8), dpi=300)

x_vals = np.linspace(-7, 7.0, 100000)
#y_vals = [dawsn_branchless(x) for x in x_vals]
initial_n = 2.0
y_vals = dawsn_branchless(x_vals, initial_n)#test
y_sum =sum_error(initial_n)

 #Create the main figure and axis for the plot
fig, ax = plt.subplots(figsize=(8, 6))#(8,6)
plt.subplots_adjust(bottom=0.25)  # Make room at the bottom for the slider
# Draw the initial line
line, = ax.plot(x_vals, y_vals, lw=1, color='blue',)
#line, = ax.plot(x_vals, y_vals, lw=1, color='blue', label=f'Relative Error around breack points = {sum_error(initial_n):.2f}')
ax.text(0.05, 0.85, f'Relative Error around breack points : {sum_error(initial_n):.2f}', 
                      transform=ax.transAxes, fontsize=11, fontweight='bold',
                      bbox=dict(facecolor='white', alpha=0.8, edgecolor='gray', boxstyle='round,pad=0.5'))

ax.set_xlabel('x')
ax.set_ylabel('y')
plt.suptitle("Branchless Dawson Code", fontsize=10, color="black", weight="bold")
ax.set_title(r'localizing cusp-like behavior at defined coordinates: (-6.25,-3.25,3.25,6.25)')
#plt.suptitle("hellow", fontsize=20, color="black", weight="bold")
ax.grid(True)
ax.legend()
# STEP 2: Make the curve incredibly thin ('lw=1.0' or 'lw=0.8') for that sharp math style
#ax.plot(x_vals, y_vals, color="#1d6fca", lw=1.0, zorder=3)

#print(dawsn_branchless(3.0,initial_n))#test
print("actual",special.dawsn(num))#relative error
plt.axvline(x=3.25,color='red',lw=1.0,linestyle='--',)
plt.axvline(x=3.25+.001,color='red',lw=1.0,linestyle='-',)
plt.axvline(x=3.25-.001,color='red',lw=1.0,linestyle='-',)
plt.axvline(x=-3.25,color='red',lw=1.0,linestyle='--',)
plt.axvline(x=-3.25+.001,color='red',lw=1.0,linestyle='-',)
plt.axvline(x=-3.25-.001,color='red',lw=1.0,linestyle='-',)
plt.axvline(x=6.25,color='red',lw=1.0,linestyle='--',)
plt.axvline(x=6.25+.001,color='red',lw=1.0,linestyle='-',)
plt.axvline(x=6.25-.001,color='red',lw=1.0,linestyle='-',)
plt.axvline(x=-6.25,color='red',lw=1.0,linestyle='--',)
plt.axvline(x=-6.25+.001,color='red',lw=1.0,linestyle='-',)
plt.axvline(x=-6.25-.001,color='red',lw=1.0,linestyle='-',)

#• plt.axvline(x=2): Draws a vertical line (goes up and down).
#• plt.axhline(y=2): Draws a horizontal line (goes left to right).
# STEP 3: FORCE THE ACCURATE 1:1 MATHEMATICAL ASPECT RATIO
# This guarantees that 1 unit on the X-axis is the exact same screen distance as 1 unit on the Y-axis
#ax.set_box_aspect(1) 
# 4. Create an axis dedicated to holding the slider geometry
# [left, bottom, width, height]
ax_slider = plt.axes([0.15, 0.1, 0.65, 0.03])

# 5. Initialize the Matplotlib Slider
n_slider = Slider(
    ax=ax_slider,
    label='increase n ',
    valmin=0.1,
    valmax=20.0,
    valinit=initial_n,
    valstep=0.1
)
# absolute error around break points
def update(val):
    current_n = n_slider.val
    # Recalculate y data based on new slider value
    line.set_ydata(dawsn_branchless(x_vals, current_n))
    # Dynamic label update
    #line.set_label(f'Relative Error around breack points = {sum_error(current_n):.2f}')#displays
    ax.text(0.05, 0.85, f'Relative Error around breack points : {sum_error(current_n):.2f}', 
                      transform=ax.transAxes, fontsize=11, fontweight='bold',
                      bbox=dict(facecolor='white', alpha=0.8, edgecolor='gray', boxstyle='round,pad=0.5'))
    
    
    ax.legend()
    # Redraw the canvas
    fig.canvas.draw_idle()
    print(dawsn_branchless(3.25+.001,current_n))
# Configure the grid limits and boundaries
ax.set_xlim([-7, 7])
ax.set_ylim([-7, 7]) # Match the limits so the viewport is balanced

# Thin out axis lines and grids to keep them subtle
ax.grid(True, which='both', color='#e5e5e5', linestyle='-', linewidth=0.5)
ax.axhline(0, color='#444444', linewidth=0.8, zorder=2)
ax.axvline(0, color='#444444', linewidth=0.8, zorder=2)
#ax.text(f'Percentage Error (MAPE): {current_n:.2f}%')
#print(dawsn_branchless(3.25+.001,current_n))
# Save as a pure mathematical vector file
plt.savefig('dawsn_perfect_desmos.svg', bbox_inches='tight', transparent=True)
#print(f"f(3.001,{current_n:.2f})=")
n_slider.on_changed(update)

plt.show()#branchless nnamdi

