import math

def calc_alphak(loading, big_dimension, small_dimension, radius, case):
    # loading: 1: tension/compression
    #          2: bending
    #          3: torsion

    # case     1: flat bar notched
    #          2: flat bar shoulder
    #          3: round bar notched
    #          4: round bar shoulder

    # flat bar notched
    if case == 1:
       if loading == 1:  # Tension
            constant_A = 0.1
            constant_B = 0.7
            constant_c = 0.13
            constant_k = 1
            constant_l = 2
            constant_m = 1.25

    elif loading == 2:  # bending
            constant_A = 0.08
            constant_B = 2.2
            constant_c = 0.2
            constant_k = 0.66
            constant_l = 2.25
            constant_m = 1.33

    # flat bar shoulder
    elif case == 2:
        if loading == 1:  # tension
            constant_A = 0.55
            constant_B = 1.1
            constant_c = 0.2
            constant_k = 0.8
            constant_l = 2.2
            constant_m = 1.33
        elif loading == 2:  # bending
            constant_A = 0.4
            constant_B = 3.8
            constant_c = 0.2
            constant_k = 0.8
            constant_l = 2.2
            constant_m = 1.33

    # Round bar notched
    elif case == 3:
        if loading == 1:  # tension/Compression
            constant_A = 0.1
            constant_B = 1.6
            constant_c = 0.11
            constant_k = 0.55
            constant_l = 2.50
            constant_m = 1.50
        elif loading == 2:  # bending
            constant_A = 0.12
            constant_B = 4
            constant_c = 0.1
            constant_k = 0.45
            constant_l = 2.66
            constant_m = 1.2
        elif loading == 3:  # torsion
            constant_A = 0.4
            constant_B = 15
            constant_c = 0.1
            constant_k = 0.35
            constant_l = 2.75
            constant_m = 1.50

    # Round bar shoulder
    elif case == 4:
        if loading == 1:  # tension/Compression
            constant_A = 0.44
            constant_B = 2
            constant_c = 0.3
            constant_k = 0.6
            constant_l = 2.2
            constant_m = 1.6
        elif loading == 2:  # bending
            constant_A = 0.4
            constant_B = 6
            constant_c = 0.8
            constant_k = 0.4
            constant_l = 2.75
            constant_m = 1.50
        elif loading == 3:  # torsion
            constant_A = 0.4
            constant_B = 25
            constant_c = 0.2
            constant_k = 0.45
            constant_l = 2.25
            constant_m = 2

    t = (big_dimension - small_dimension) / 2
    a = small_dimension / 2
    rho = radius

    rt1 = constant_A / ((t / rho) ** constant_k)
    rt2 = constant_B*((1+a/rho)/(a/rho*math.sqrt(a/rho)))**constant_l
    rt3 = constant_c * (a / rho) / ((a / rho + t / rho) * (t / rho) ** constant_m)

    root = math.sqrt(rt1 + rt2 + rt3)
    alphak = 1 + 1 / root

    return alphak

def calk_alphak_ld1(bigd, smalld, t, r):
    bigd = 20
    smalld = bigd * d_to_D
    t = (bigd - smalld) / 2
    r = r_to_t * t
    dom_1 = np.sqrt(0.62 * r / t + 7 * r / smalld * (1 + 2 * r / smalld) ** 2)
    alphak_ld_1 = 1 + (1 / dom_1)
    return alphak_ld_1

def calk_alphak_ld2(bigd, smalld, t, r):
    bigd = 20
    smalld = bigd * d_to_D
    t = (bigd - smalld) / 2
    r = r_to_t * t
    dom_2 = np.sqrt(0.62 * r / t + 11.6 * r / smalld * (1 + 2 * r / smalld) ** 2 + 0.2 * (r / t) ** 3 * smalld / bigd)
    alphak_ld_2 = 1 + (1 / dom_2)
    return alphak_ld_2

def calk_alphak_ld3(bigd, smalld, t, r):
    bigd = 20
    smalld = bigd * d_to_D
    t = (bigd - smalld) / 2
    r = r_to_t * t
    dom_3 = np.sqrt(3.4 * r / t + 38 * r / smalld * (1 + 2 * r / smalld) ** 2 + 1.0 * (r / t) ** 2 * smalld / bigd)
    alphak_ld_3 = 1 + (1 / dom_3)
    return alphak_ld_3


#TODO:
import numpy as np
import matplotlib.pyplot as plt

# fix outer diameter
bigd = 20

d_to_D_relation = np.linspace(0.4, 0.9, 20)
r_to_t_relation = np.array(
  [0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 1, 1.5, 2, 2.5, 5, 10])
np.seterr(invalid='ignore')

#preallocate alpha_array(20,20) _ALL ARRAYS_
alpha_array = np.zeros((20, 20))

#calculate alphak using equation for round bar shoulder
# Outer loop for r_to_t
for index_ld1_clk, r_to_t in enumerate(r_to_t_relation):
    # Inner loop for d_to_D
  for j_ld1_clk, d_to_D in enumerate(d_to_D_relation):
       smalld = bigd * d_to_D
       t = (bigd - smalld) / 2
       r = r_to_t * t
       #calculate alphak for round bar shoulder case 4 loading 1 tension/compression
       alphak = calc_alphak(1, bigd, smalld, r, 4)
       alpha_array[index_ld1_clk, j_ld1_clk] = alphak
       print("alphak_calculation_ld1_ = ", alphak)
plt.plot(d_to_D_relation, np.transpose(alpha_array))
plt.xlabel('d_to_D_relation')
plt.ylabel('alphak_calculation_ld1')
plt.show()

# Outer loop for r_to_t
for index_ld2_clk, r_to_t in enumerate(r_to_t_relation):
    # Inner loop for d_to_D
  for j_ld2_clk, d_to_D in enumerate(d_to_D_relation):
       smalld = bigd * d_to_D
       t = (bigd - smalld) / 2
       r = r_to_t * t
       # calculate alphak for round bar shoulder case 4 loading 2 bending
       alphak = calc_alphak(2, bigd, smalld, r, 4)
       alpha_array[index_ld2_clk, j_ld2_clk] = alphak
       print("alphak_calculation_ld2_ = ", alphak)
plt.plot(d_to_D_relation, np.transpose(alpha_array))
plt.xlabel('d_to_D_relation')
plt.ylabel('alphak_calculation_ld2_')
plt.show()

# Outer loop for r_to_t
for index_ld3_clk, r_to_t in enumerate(r_to_t_relation):
    # Inner loop for d_to_D
  for j_ld3_clk, d_to_D in enumerate(d_to_D_relation):
       smalld = bigd * d_to_D
       t = (bigd - smalld) / 2
       r = r_to_t * t
       # calculate alphak for round bar shoulder case 4 loading 3 torsion
       alphak = calc_alphak(3, bigd, smalld, r, 4)
       alpha_array[index_ld3_clk, j_ld3_clk] = alphak
       print("alphak_calculation_ld3_ = ", alphak)
plt.plot(d_to_D_relation, np.transpose(alpha_array))
plt.xlabel('d_to_D_relation')
plt.ylabel('alphak_calculation_ld3_')
plt.show()

#calculate alphak using graph for round bar shoulder
#calculate alpha for round bar shoulder tension/compression loading using the graph equation loading = 1 , case = 4
# Outer loop for r_to_t
for index_ld1, r_to_t in enumerate(r_to_t_relation):
    # Inner loop for d_to_D
  for j_ld1, d_to_D in enumerate(d_to_D_relation):
       smalld = bigd * d_to_D
       t = (bigd - smalld) / 2
       r = r_to_t * t
       # call function calk_alphak_ld1
       alphak_ld_1 = calk_alphak_ld1(bigd, smalld, t, r)
       #alpha_array[index_ld1, j_ld1] = alphak_ld_1
       print("alphak_ld_1 = ", alphak_ld_1)
plt.plot(d_to_D_relation, np.transpose(alpha_array))
plt.xlabel('d_to_D_relation')
plt.ylabel('alphak_ld_1')
plt.show()

#calculate alpha for round bar shoulder bending loading using the graph equation loading = 2 , case = 4
# Outer loop for r_to_t
for index_ld2, r_to_t in enumerate(r_to_t_relation):
    # Inner loop for d_to_D
  for j_ld2, d_to_D in enumerate(d_to_D_relation):
       smalld = bigd * d_to_D
       t = (bigd - smalld) / 2
       r = r_to_t * t
       #call function
       alphak_ld2 = calk_alphak_ld2(bigd, smalld, t, r)
       #alpha_array[index_ld2, j_ld2] = alpha_ld_2
       print("alphak_ld_2 = ", alphak_ld2)
plt.plot(d_to_D_relation, np.transpose(alpha_array))
plt.xlabel('d_to_D_relation')
plt.ylabel('alphak_ld_2')
plt.show()

# calculate alpha for round bar shoulder Torsion loading using the graph equation loading = 3 , case = 4
# Outer loop for r_to_t
for index_ld3, r_to_t in enumerate(r_to_t_relation):
    # Inner loop for d_to_D
  for j_ld3, d_to_D in enumerate(d_to_D_relation):
       smalld = bigd * d_to_D
       t = (bigd - smalld) / 2
       r = r_to_t * t
       #call function calk_aplhak_ld3
       alphak_ld_3 = calk_alphak_ld3(bigd, smalld, t, r)
       #alpha_array[index_ld3, j_ld3] = alphak_ld_3
       print("alphak_ld_3 = ", alphak_ld_3)
plt.plot(d_to_D_relation, np.transpose(alpha_array))
plt.xlabel('d_to_D_relation')
plt.ylabel('alphak_ld_3')
plt.show()

#calculate rated stress slope G' (1/mm)
#G' for a flat & round bar shoulder case

def clk_G_slope(case, loading):
      # loading: 1: tension/compression or bending
      #          2: torsion

      # case:    1: flat bar (notched)
      #          2: round bar (notched)
      #          3: flat bar (shoulder)
      #          4: round bar (shoulder)
      #          5: flat bar (not notched)

      # flat bar (notched)
      if case == 1:             # 1: flat bar (notched)
         if loading == 1:       # loading: 1: tension/compression or bending

             bigd = 20
             smalld = bigd * d_to_D
             t = (bigd - smalld) / 2
             r = r_to_t * t
             if t <= 0.5:
                 phi = np.sprt(8 * (bigd - smalld / r)) + 2
             else:
                 phi = 0
             r = r_to_t * t
             G_slope = 2 / r * (1 + phi)
         elif loading == 2:       #flat bar (notched) loading (torsion)

             G_slope = 0


      elif case == 2:          #round bar (notched)
         if loading == 1:      #loading (tension/compression, bending)
             bigd = 20
             smalld = bigd * d_to_D
             t = (bigd - smalld) / 2
             r = r_to_t * t
             if t <= 0.5:
                 phi = np.sprt(8 * (bigd - smalld / r)) + 2
             else:
                 phi = 0
             r = r_to_t * t
             G_slope = 2 / r * (1 + phi)
         elif loading == 2:   #loading torsion

             bigd = 20
             smalld = bigd * d_to_D
             t = (bigd - smalld) / 2
             r = r_to_t * t
             if t <= 0.5:
                 phi = np.sprt(8 * (bigd - smalld / r)) + 2
             else:
                 phi = 0
             r = r_to_t * t
             G_slope = 1 / r

      elif case == 3:  #flat bar shoulder
         if loading == 1: #loading tension/compression, bending
             bigd = 20
             smalld = bigd * d_to_D
             t = (bigd - smalld) / 2
             r = r_to_t * t
             if t <= 0.5:
                 phi = np.sprt(8 * (bigd - smalld / r)) + 2
             else:
                 phi = 0
             r = r_to_t * t
             G_slope = 2.3 / r * (1 + phi)

         elif loading == 2: #loading torsion

             G_slope = 0

      elif case == 4: #round bar shoulder
         if loading == 1: #loading tension/compression, bending
             bigd = 20
             smalld = bigd * d_to_D
             t = (bigd - smalld) / 2
             r = r_to_t * t
             if t <= 0.5:
                 phi = np.sprt(8 * (bigd - smalld / r)) + 2
             else:
                 phi = 0
             r = r_to_t * t
             G_slope = 2.3 / r * (1 + phi)

         elif loading == 2: #loading torsion
             bigd = 20
             smalld = bigd * d_to_D
             t = (bigd - smalld) / 2
             r = r_to_t * t
             if t <= 0.5:
                 phi = np.sprt(8 * (bigd - smalld / r)) + 2
             else:
                 phi = 0
             r = r_to_t * t
             G_slope = 1.15 / r
      elif case == 5:  #  bar (not notched)
         if loading == 1: #loading tension/compression, bending
             bigd = 20
             smalld = bigd * d_to_D
             t = (bigd - smalld) / 2
             r = r_to_t * t
             if t <= 0.5:
                 phi = np.sprt(8 * (bigd - smalld / r)) + 2
             else:
                 phi = 0
             r = r_to_t * t
             G_slope = 2 / smalld

         elif loading == 2: #loading torsion
             bigd = 20
             smalld = bigd * d_to_D
             t = (bigd - smalld) / 2
             r = r_to_t * t
             if t <= 0.5:
                 phi = np.sprt(8 * (bigd - smalld / r)) + 2
             else:
                 phi = 0
             r = r_to_t * t
             G_slope = 2 / smalld

      return G_slope


#calculate support figure n
# For bending, tension/compression loading
def calc_nb(clk_G_slope_case, Rp_0_2):
    if clk_G_slope_case ==1: #flat bar (notcehd)
        G_slope_bending= clk_G_slope(1,1)
        n_b = 1 + np.sqrt(G_slope_bending) * 10 ** (-0.33 - (Rp_0_2 / 712))
        print("Support Figure n-b for flat bar (notched): ", n_b)
    elif clk_G_slope_case ==2: #round bar (notcehd)
        G_slope_bending = clk_G_slope(2, 1)
        n_b = 1 + np.sqrt(G_slope_bending) * 10 ** (-0.33 - (Rp_0_2 / 712))
        print("Support Figure n-b for round bar (notched): ", n_b)
    elif clk_G_slope_case ==3: #flat bar (shoulder)
        G_slope_bending = clk_G_slope(3, 1)
        n_b = 1 + np.sqrt(G_slope_bending) * 10 ** (-0.33 - (Rp_0_2 / 712))
        print("Support Figure n-b for flat bar (shoulder): ", n_b)
    elif clk_G_slope_case ==4: #round bar (shoulder)
        G_slope_bending = clk_G_slope(4, 1)
        n_b = 1 + np.sqrt(G_slope_bending) * 10 ** (-0.33 - (Rp_0_2 / 712))
        print("Support Figure n-b for round bar (shoulder): ", n_b)
    elif clk_G_slope_case ==5: #flat bar (not notched)
        G_slope_bending = clk_G_slope(5, 1)
        n_b = 1 + np.sqrt(G_slope_bending) * 10 ** (-0.33 - (Rp_0_2 / 712))
        print("Support Figure n-b for flat bar (not notched): ", n_b)
    return n_b
#For torsion loading
def calc_nt(clk_G_slope_case, Rp_0_2):
        if clk_G_slope_case == 1:  # flat bar (notcehd)
            G_slope_torsion = clk_G_slope(1, 2)
            n_t = 1 + np.sqrt(G_slope_torsion) * 10 ** (-0.33 - (Rp_0_2 / 712))
            print("Support Figure n-t for flat bar (notched): ", n_t)
        elif clk_G_slope_case == 2:  # round bar (notcehd)
            G_slope_torsion = clk_G_slope(2, 2)
            n_t = 1 + np.sqrt(G_slope_torsion) * 10 ** (-0.33 - (Rp_0_2 / 712))
            print("Support Figure n-t for round bar (notched): ", n_t)
        elif clk_G_slope_case == 3:  # flat bar (shoulder)
            G_slope_torsion = clk_G_slope(3, 2)
            n_t = 1 + np.sqrt(G_slope_torsion) * 10 ** (-0.33 - (Rp_0_2 / 712))
            print("Support Figure n-t for flat bar (shoulder): ", n_t)
        elif clk_G_slope_case == 4:  # round bar (shoulder)
            G_slope_torsion = clk_G_slope(4, 2)
            n_t = 1 + np.sqrt(G_slope_torsion) * 10 ** (-0.33 - (Rp_0_2 / 712))
            print("Support Figure n-t for round bar (shoulder): ", n_t)
        elif clk_G_slope_case == 5:  # flat bar (not notched)
            G_slope_torsion = clk_G_slope(5, 2)
            n_t = 1 + np.sqrt(G_slope_torsion) * 10 ** (-0.33 - (Rp_0_2 / 712))
            print("Support Figure n-t for flat bar (not notched): ", n_t)
        return n_t


G_slope_relation = np.linspace(0, 6, 15)
Rp_relation = np.array(
  [200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 800, 900, 1000, 1200])
np.seterr(invalid='ignore')
n_array = np.zeros((15, 15))


for i, Rp in enumerate(Rp_relation):

  for j, G_slope in enumerate(G_slope_relation):
       #calculate support figure n
       n = 1 + np.sqrt(G_slope) * 10 ** (-0.33 - (Rp / 712))
       n_array[i, j] = n
       print("support_figure_factor = ", n)
plt.plot(G_slope_relation, np.transpose(n_array))
plt.xlabel('G_rated_slope')
plt.ylabel('support_figure_factor')
plt.show()

#Example case solving procedure
def calc_alphak_dash(loading, case):
    if loading == 1:   #bending
        if case == 1:  #flat bar (notched)
           alphak_bending= calc_alphak(2, 35, 30, 8, 1)  #flat bar (notched), loading bending
           print("alphak-bending = ", alphak_bending)
           support_figure_bending= calc_nb(1, 1050)
           alphak_dash_bending = alphak_bending / support_figure_bending
           alphak_dash = alphak_dash_bending
           print("Alphak_dash_bending =", alphak_dash_bending)
        elif case == 2: #round bar (notched)
            alphak_bending = calc_alphak(2, 35, 30, 8, 3)  #round bar (notched), loading bending
            print("alphak-bending = ", alphak_bending)
            support_figure_bending = calc_nb(2, 1050)
            alphak_dash_bending = alphak_bending / support_figure_bending
            alphak_dash = alphak_dash_bending
            print("Alphak_dash_bending =", alphak_dash_bending)
        elif case == 3:  #flat bar (shoudler)
            alphak_bending = calc_alphak(2, 35, 30, 8, 2)  # flat bar (shoudler), loading bending
            print("alphak-bending = ", alphak_bending)
            support_figure_bending = calc_nb(3, 1050)
            alphak_dash_bending = alphak_bending / support_figure_bending
            alphak_dash = alphak_dash_bending
            print("Alphak_dash_bending =", alphak_dash_bending)
        elif case == 4:  #round bar (shoulder) example case
            alphak_bending = calc_alphak(2, 35, 30, 8, 3)  # round bar (shoudler), loading bending
            print("alphak-bending = ", alphak_bending)
            support_figure_bending = calc_nb(4, 1050)
            alphak_dash_bending = alphak_bending / support_figure_bending
            alphak_dash = alphak_dash_bending
            print("Alphak_dash_bending =", alphak_dash_bending)

    elif loading == 2:  #loading torsion
        if case == 1:   #flat bar (notched)
           alphak_torsion = calc_alphak(2, 35, 30, 8, 1)  # flat bar (notched), loading torsion
           print("alphak-torsion = ", alphak_torsion)
           support_figure_torsion = calc_nt(1, 1050)
           alphak_dash_torsion = alphak_torsion / support_figure_torsion
           alphak_dash = alphak_dash_torsion
           print("Alphak_dash_torsion =", alphak_dash_torsion)
        elif case == 2:   #round bar (notched)
           alphak_torsion = calc_alphak(2, 35, 30, 8, 4)  # round bar (notched), loading torsion
           print("alphak-torsion = ", alphak_torsion)
           support_figure_torsion = calc_nt(2, 1050)
           alphak_dash_torsion = alphak_torsion / support_figure_torsion
           alphak_dash = alphak_dash_torsion
           print("Alphak_dash_torsion =", alphak_dash_torsion)
        elif case == 3:  #flat bar (shoulder)
            alphak_torsion = calc_alphak(2, 35, 30, 8, 2)  # flat bar (shoudler), loading torsion
            print("alphak-torsion = ", alphak_torsion)
            support_figure_torsion = calc_nt(3, 1050)
            alphak_dash_torsion = alphak_torsion / support_figure_torsion
            alphak_dash = alphak_dash_torsion
            print("Alphak_dash_torsion =", alphak_dash_torsion)
        elif case == 4:  #round bar (shoudler)
            alphak_torsion = calc_alphak(2, 35, 30, 8, 3)  # round bar (shoudler), loading torsion  example case
            print("alphak-torsion = ", alphak_torsion)
            support_figure_torsion = calc_nt(4, 1050)
            alphak_dash_torsion = alphak_torsion / support_figure_torsion
            alphak_dash = alphak_dash_torsion
            print("Alphak_dash_torsion =", alphak_dash_torsion)
    return alphak_dash






#calculation for notch senstivity factor
def func(r, Rm, Rp0_2):
    dom_nk = 1 + 8 / r * ( 1 - Rm / Rp0_2)**3
    nk = 1 / dom_nk
    print("notch sensitivity figure nk =", nk)
    return nk

#Stress concentration factor calculation
def func_beta_K(loading, case):
    # loading: 1: tension/compression, bending
    #          2: torsion

    # case     1: flat bar (notched)
    #          2: round bar (notched)
    #          3: flat bar (shoulder)
    #          4: round bar (shoulder)
    #          5: flat bar (not notched)

    if loading == 1:   #bending, tension/compression
        if case == 1:   #flat bar (notched)
            alphak = calc_alphak(1,35,30,8,1)
            n_b = calc_nb(1,1050)
            alphak_dash_tension_compression = alphak/n_b
            nk = func(8, 1450, 1050)
            beta_k = (1 + (alphak_dash_tension_compression - 1) * nk)
        elif case == 2:  #round bar (notched)
            alphak = calc_alphak(2, 35, 30, 8, 4)
            n_b = calc_nb(2, 1050)
            alphak_dash_bending = alphak / n_b
            nk = func(8, 1450, 1050)
            beta_k = (1 + (alphak_dash_bending - 1) * nk)
        elif case == 3:   #flat bar (shoulder)
            alphak = calc_alphak(2, 35, 30, 8, 2)
            n_b = calc_nb(3, 1050)
            alphak_dash_bending = alphak / n_b
            nk = func(8, 1450, 1050)
            beta_k = (1 + (alphak_dash_bending - 1) * nk)
        elif case == 4:  #round bar (shoulder)
            alphak = calc_alphak(2, 35, 30, 8, 3)
            n_b = calc_nb(4, 1050)
            alphak_dash_bending = alphak / n_b
            nk = func(8, 1450, 1050)
            beta_k = (1 + (alphak_dash_bending - 1) * nk)
        elif case == 5:   #flat bar (not notched)
            alphak = calc_alphak(2, 35, 30, 8, 1)      #for alphak there is no case for flat bar (not notched), so I considered it as case 1
            n_b = calc_nb(5, 1050)
            alphak_dash_bending = alphak / n_b
            nk = func(8, 1450, 1050)
            beta_k = (1 + (alphak_dash_bending - 1) * nk)
    elif loading == 2:  #torsion
        if case == 1:  #flat bar (notched)
            alphak = calc_alphak(3, 35, 30, 8, 1)
            n_t = calc_nt(1, 1050)
            alphak_dash_torsion = alphak / n_t
            nk = func(8, 1450, 1050)
            beta_k = (1 + (alphak_dash_torsion - 1) * nk)
        elif case == 2:   #round bar (notched)
            alphak = calc_alphak(3, 35, 30, 8, 3)
            n_t = calc_nt(2, 1050)
            alphak_dash_torsion = alphak / n_t
            nk = func(8, 1450, 1050)
            beta_k = (1 + (alphak_dash_torsion - 1) * nk)
        elif case == 3:   #flat bar (shoulder)
            alphak = calc_alphak(3, 35, 30, 8, 2)
            n_t = calc_nt(3, 1050)
            alphak_dash_torsion = alphak / n_t
            nk = func(8, 1450, 1050)
            beta_k = (1 + (alphak_dash_torsion - 1) * nk)
        elif case == 4:   #round bar (shoulder)
            alphak = calc_alphak(3, 35, 30, 8, 3)
            n_t = calc_nt(4, 1050)
            alphak_dash_torsion = alphak / n_t
            nk = func(8, 1450, 1050)
            beta_k = (1 + (alphak_dash_torsion - 1) * nk)
        elif case == 5:    #flat bar (not notched)
            alphak = calc_alphak(3, 35, 30, 8, 1)      #took alphak for flat bar notched
            n_t = calc_nt(5, 1050)
            alphak_dash_torsion = alphak / n_t
            nk = func(8, 1450, 1050)
            beta_k = (1 + (alphak_dash_torsion - 1) * nk)
    return beta_k


#notch size factor
alphak_bending= calc_alphak(2, 35, 30, 8, 4)  #round bar shoulder, loading bending
print("alphak-bending = ", alphak_bending)
alphak_torsion= calc_alphak(3, 35, 30, 8, 4)  #round bar shoulder, loading torsion
print("alphak-torsion = ", alphak_torsion)
#notch sensitivity factor
func(8, 1450, 1050)
#Rated stress slopes
G_slope_bending = clk_G_slope(4, 1)
print("Gb'=",G_slope_bending)
G_slope_torsion = clk_G_slope(4, 2)
print("Gt'=",G_slope_torsion)
#support figure
support_figure_bending= calc_nb(4, 1050)
print("Support figure for bending", support_figure_bending)
support_figure_torsion= calc_nt(4, 1050)
print("Support figure for torsion", support_figure_torsion)

alphak_dash_bending = calc_alphak_dash(1,4)  #bending loading, round bar shoudler
alphak_dash_torsion = calc_alphak_dash(2,4)  #torsion loading, round bar shoudler
print("Alphak_dash_bending =",alphak_dash_bending)
print("Alphak_dash_torsion =",alphak_dash_torsion)

#Solving the example case for beta_K_bending for round bar shoudler fillet (stepping shaft)
#for bending
beta_k_bending=func_beta_K(1,4)
print("Stress concentration factor for bending =", beta_k_bending)
#for torsion
beta_k_torsion=func_beta_K(2, 4)
print("Stress concentration factor for torsion =", beta_k_torsion)










