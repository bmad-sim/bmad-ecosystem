from pytao import Tao
import numpy as np
import matplotlib.pyplot as plt
import os
import time
import math
import sys
import csv



def generate_spin_resonant_tunes(NU, Ki, Kf):
    # Initialize a dictionary to hold valid K values and their annotations
    valid_K_values_with_annotations = {}

    # Estimate a reasonable range for N
    N_min = int(Ki - abs(NU)) - 1
    N_max = int(Kf + abs(NU)) + 1

    for N in range(N_min, N_max):
        # Calculate potential K values and check if they meet the criteria
        K_values = [
            (NU + N, 2 * math.pi),
            (abs(NU - N), -2 * math.pi if N > NU else 2 * math.pi),
            (N + NU, 2 * math.pi),
        ]

        for K, annotation in K_values:
            if Ki <= K <= Kf and K > 0 and K not in valid_K_values_with_annotations:
                valid_K_values_with_annotations[K] = annotation

    # Sort the dictionary by key (K value) to return a sorted list of tuples
    sorted_K_values_and_annotations = sorted(valid_K_values_with_annotations.items())

    # Return the sorted list of unique K values and their annotations
    return sorted_K_values_and_annotations

def generate_strengths(inf,out,plot_file,particle,SPRINT,lo,hi,J):
    tao=Tao(f'-lat {inf} -noplot')
    #tao.cmd('set element overlay::* is_on = F')
    #tao.cmd('set element * field_master = F')
    if SPRINT:
        tao.cmd('set ele * spin_tracking_method = sprint')
    tao.cmd('set bmad_com spin_tracking_on = T')
    #tao.cmd('show -write bmad.twiss lat beginning:end -att l -att ang -att beta_a -att alpha_a -att beta_b -att alpha_b  -att phi_b/(2*pi) -att orbit_y -att k1*l -att k2*l -att h1 -att h2 -att orbit_x -att orbit_px -att orbit_py')
    data = tao.cmd('show val lat::tune.b[0]/twopi')
    Qy =float(data[0])
    Kn = generate_spin_resonant_tunes(Qy, lo, hi)

    strengths = []

    for i in range(0,len(Kn)):
        tao.cmd(f'set element 0 e_tot = {Kn[i][0]} / anomalous_moment_of({particle}) * mass_of({particle})')
        s = tao.cmd('python spin_resonance')
        d1 = float(s[6][s[6].find('dq_b_sum;REAL;F;')+len('dq_b_sum;REAL;F;'):])
        d2 = float(s[7][s[7].find('dq_b_diff;REAL;F;')+len('dq_b_diff;REAL;F;'):])
        e1 = float(s[8][s[8].find('xi_res_b_sum;REAL;F;')+len('xi_res_b_sum;REAL;F;'):])
        e2 = float(s[9][s[9].find('xi_res_b_diff;REAL;F;')+len('xi_res_b_diff;REAL;F;'):])
        if abs(d1) > abs(d2):
            strengths.append(e2*np.sqrt(J))
        else:
            strengths.append(e1*np.sqrt(J))

    lst = zip([Kn[i][0] for i in range(len(Kn))],strengths)
    open(out, 'w').close()
    with open(out,"w") as f:
        writer = csv.writer(f,delimiter='\t')
        writer.writerows(lst)

    plt.scatter([k[0] for k in Kn],strengths)
    plt.xlabel(r"$G\gamma$")
    plt.ylabel("$\epsilon$")
    plt.savefig(plot_file, dpi=1200)


def main():
    inf = sys.argv[1]
    out = sys.argv[2]
    plot_file = sys.argv[3]
    particle = sys.argv[4]
    SPRINT = True if sys.argv[5].upper() == "T" else False
    lo = float(sys.argv[6])
    hi = float(sys.argv[7])
    J = float(sys.argv[8])
    generate_strengths(inf,out,plot_file,particle,SPRINT,lo,hi,J)



if __name__ == "__main__":
    main()
