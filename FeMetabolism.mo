within ;
package FeMetabolism
  model FeMetabolismModel
    "Implementation of Enulescu et al. PLOS Comput. Biol. 2017"

    //################################################
    // Main varibles in ODES
    //################################################

    // Hepcidin expression
    //################################################
    Real hep(
      start=0.66540,
      fixed=true,
      unit="ug") "Hepcidin ammount (s16)";
    Real r_hep_syn(unit="ug.h-1") "Hepcidin synthesis rate";
    Real r_hep_deg(unit="ug.h-1") "Hepcidin degradation rate";

    Real Bmp6(
      start=15.834,
      fixed=true,
      unit="ug") "Bmp6 ammount (s18)";
    Real r_Bmp6_syn(unit="ug.h-1") "Bmp6 synthesis rate";
    Real r_Bmp6_deg(unit="ug.h-1") "Bmp6 degradation rate";

    Real LPS(
      start=1.000,
      fixed=true,
      unit="ug") "LPS (Lipopolysaccharide) ammount (s23)";

    Real Il6mRNA(
      start=0.000,
      fixed=true,
      unit="ug") "IL-6 mRNA ammount (s21)";
    Real r_Il6mRNA_syn(unit="ug.h-1") "IL-6 mRNA synthesis rate";
    //NESEDI JEDNOTKY !!!
    Real r_Il6mRNA_deg(unit="ug.h-1") "IL-6 mRNA degradation rate";

    Real Il6(
      start=0.000,
      fixed=true,
      unit="ug") "IL-L ammount (s22)";
    Real r_Il6_syn(unit="ug.h-1") "IL-6 synthesis rate";
    Real r_Il6_deg(unit="ug.h-1") "IL-6 degradation rate";

    // Ferroportin regulation
    //################################################
    Real Fpn_liv_mRNA(
      start=0.92241,
      fixed=true,
      unit="ug") "Fpn mRNA ammount in liver (s7)";
    Real Fpn_spl_mRNA(
      start=0.92241,
      fixed=true,
      unit="ug") "Fpn mRNA ammount in spleen (s9)";
    Real Fpn_duo_mRNA(
      start=0.92241,
      fixed=true,
      unit="ug") "Fpn mRNA ammount in duodenum (s8)";
    Real Fpn_res_mRNA(
      start=0.92241,
      fixed=true,
      unit="ug") "Fpn mRNA ammount in other organs (s11)";

    Real Fpn_liv(
      start=1.000,
      fixed=true,
      unit="ug") "Fpn ammount in liver (s12)";
    Real Fpn_spl(
      start=1.000,
      fixed=true,
      unit="ug") "Fpn ammount in spleen (s14)";
    Real Fpn_duo(
      start=1.000,
      fixed=true,
      unit="ug") "Fpn ammount in duodenum (s13)";
    Real Fpn_res(
      start=1.000,
      fixed=true,
      unit="ug") "Fpn ammount in other organs (s15)";

    // Dynamics of the iron pools
    //################################################
    Real Fe_liv(
      start=76.998,
      fixed=true,
      unit="ug") "Fe ammount in liver (s2)";
    Real Fe_ser_liv(unit="ug.h-1") "Fe flow from serum to liver";
    Real Fe_liv_ser(unit="ug.h-1") "Fe flow from liver to serum";

    Real Fe_spl(
      start=17.741,
      fixed=true,
      unit="ug") "Fe ammount in spleen (s5)";
    Real Fe_bm(
      start=63.171,
      fixed=true,
      unit="ug") "Fe ammount in bones (s3)";
    Real Fe_RBC(
      start=1016.900,
      fixed=true,
      unit="ug") "Fe ammount in red blood cells (s6)";
    Real Fe_duo(
      start=2.9726,
      fixed=true,
      unit="ug") "Fe ammount in duodenum (s4)";
    Real Fe_res(
      start=466.990,
      fixed=true,
      unit="ug") "Fe ammount in other organs (s10)";

    Real Fe_ser(
      start=1.51330,
      fixed=true,
      unit="ug") "Fe ammount in serum (s1)";
    Real Fe_input(unit="ug.h-1") "Fe ammount in serum: input";
    Real Fe_output(unit="ug.h-1") "Fe ammount in serum: output";

    //################################################
    // Model parameters
    // Format: rate constants (k): k^{X}_{Y} -->> k_{X}_{Y}; (K): K^{X}_{Y} -->> K_{X}_{Y}; MIN/MAX value -->> MIN/MAX in description
    //################################################

    // Hepcidin expression
    //################################################
    parameter Real k_hep_deg(unit="h-1") = 0.07
      "Hepcidin degradation rate (k20), 0.067/0.070";
    parameter Real v_Bmp6_max(unit="h-1") = 31.47
      "Bmp6 maximal synthesis rate (k21*k37 = k21*K_Bmp6) k21 = 1.6015, K_Bmp6 = 19.65, 14.2/126.5";
    parameter Real K_Bmp6(unit="ug") = 19.65
      "Michaelis-Menten constant Bmp6 synthesis (k37), 16.5/55.7";
    parameter Real Tf(unit="ug") = 1000.00
      "Paremeter determining the maximal amount of iron that can be bound to transferrin (k32)";
    parameter Real k_Bmp6_deg(unit="h-1") = 2.3962
      "Bmp6 degradation rate (k22), 1.0/9.5";
    parameter Real k_LPS_deg(unit="h-1") = 5.8560
      "LPS degradation rate, 5.9/5.9";
    parameter Real K_Il6mRNA(unit="ug") = 2.6e-6
      "Michaelis-Menten constant Il6mRNA synthesis, 2.6-e6/2.6e-6";
    parameter Real k_Il6mRNA_deg(unit="h-1") = 0.2814
      "Il6mRNA degradation rate, 0.28/0.28";
    parameter Real k_Il6_syn(unit="h-1") = 646.3946
      "Il6 synthesis rate, 4.1067*k23, k23 = 157.4, 136/872";
    // NESEDI JEDNOTKY!!!
    parameter Real k_Il6_deg(unit="h-1") = 4.4465
      "Il6 degradation rate, 4.45/4.45";

    // Ferroportin regulation
    //################################################
    parameter Real K_2(unit="a.u.") = 1.2836e-3
      "Constant FpnmRNA production (k41), 3.0e-4/2.0e-3";
    parameter Real k_FpnmRNA_deg(unit="h-1") = 1.0841
      "FpnmRNA degradation rate (k25), 1.02/1.20";
    parameter Real K_liv_1(unit="a.u.") = 30.6600
      "Constant Fpn_liv_mRNA production (k24), 28.0/32.5";
    parameter Real K_spl_1(unit="a.u.") = 33.3182
      "Constant Fpn_spl_mRNA production (k36*K_liv_1, k36 = 1.086700), 33.1/34.5";
    parameter Real K_duo_1(unit="a.u.") = 0.6132
      "Constant Fpn_duo_mRNA production (k35*K_liv_1, k35 = 0.020001), 0.56/0.79";
    parameter Real K_res_1(unit="a.u.") = 11.2305
      "Constant Fpn_res_mRNA production (k42*K_liv_1, k42 = 0.366290), 7.80/43.7";
    parameter Real k_Fpnliv_syn(unit="h-1") = 0.1297
      "Fpnliv synthesis rate (k28), 0.07/0.14";
    parameter Real k_Fpnspl_syn(unit="h-1") = 0.0227
      "Fpnspl synthesis rate (k30), 0.015/0.027";
    parameter Real k_Fpnduo_syn(unit="h-1") = 0.030299
      "Fpnduo synthesis rate (k29), 0.01/0.25";
    parameter Real k_Fpnres_syn(unit="h-1") = 0.0050598
      "Fpnres synthesis rate (k40), 0.004/0.108";
    parameter Real k_liv_1(unit="ug-1") = 0.0033177
      "Constant Fpn_liv production (k17), 0.002/0.006";
    parameter Real k_spl_1(unit="ug-1") = 0.014027
      "Constant Fpn_spl production (k27), 0.005/0.028";
    parameter Real k_duo_1(unit="ug-1") = 0.16007
      "Constant Fpn_duo production (k46), 0.06/0.49";
    parameter Real k_res_1(unit="ug-1") = 0.11376
      "Constant Fpn_res production (k47), 0.004/0.152";
    parameter Real k_Fpnliv_deg(unit="h-1") = 0.055363
      "Fpnliv degradation rate (k12), 0.01/0.06";
    parameter Real k_Fpnspl_deg(unit="h-1") = 0.003024
      "Fpnspl degradation rate (k12*k18, k18 = 0.054621), 0.0007/0.0038";
    parameter Real k_Fpnduo_deg(unit="h-1") = 0.021121
      "Fpnduo degradation rate (k12*k15, k15 = 0.3815), 0.0056/0.147";
    parameter Real k_Fpnres_deg(unit="h-1") = 0.028917
      "Fpnres degradation rate (k12*k38, k38 = 0.52232), 0.025/0.129";
    parameter Real k_liv_2(unit="a.u.") = 2.574300
      "Constant Fpnliv degradation (k13), 2.11/12.93";
    parameter Real k_spl_2(unit="a.u.") = 11.505576
      "Constant Fpnspl degradation (k13*k19, k19 = 4.4694), 9.2/73.7";
    parameter Real k_duo_2(unit="a.u.") = 1.432109
      "Constant Fpnduo degradation (k13*k16, k16 = 0.55631), 0.78/4.16";
    parameter Real k_res_2(unit="a.u.") = 11.626311
      "Constant Fpnres degradation (k13*k39, k39 = 4.5163), 3.7/38.1";

    // Dynamics of the iron pools
    //################################################
    parameter Real v_liv_1(unit="ug.h-1") = 3.9607
      "Low liver iron uptake (k4), 2.78/9.81";
    parameter Real v_liv_2(unit="ug.h-1") = 14.3810
      "High liver iron uptake (k45) - !!! Inconsistecy xml model vs. paper !!!, xml model value adopted";
    parameter Real th(unit="ug") = 2.6870
      "Threshold serum iron value (k26), 2.08/3.00";
    parameter Real u_liv(unit="h-1") = 0.07784
      "Liver iron export rate (k1), 0.05/0.19";
    parameter Real Fe_liv_max(unit="uq") = 119.55
      "Threshold value liver iron export (k43), 100.0/159.0";
    parameter Real v_spl_1(unit="h-1") = 0.0036035
      "Spleen iron uptake rate from RBC (k10), 0.002 = value for trace exp., 0.004/0.005";
    parameter Real v_spl_2(unit="h-1") = 0.0096817
      "Spleen iron uptake rate from bones (k7), 0.008/0.019";
    parameter Real u_spl(unit="h-1") = 0.24102
      "Spleen export rate (k3), 0.21/0.36";
    parameter Real Fe_spl_max(unit="ug") = 88.216
      "Threshold value spleen iron export (k44), 57.0/95.0";
    parameter Real v_bm(unit="h-1") = 2.8256
      "Bone marrow uptake rate (k5), 2.66/4.16";
    parameter Real v_RBC(unit="h-1") = 0.058008
      "RBC uptake rate (k8), 0.055/0.075";
    // !!! Problem u rovnice (14) - doplnit prvni term na prave strane
    parameter Real k14(unit="1") = 0.19419 "k14 constant";
    parameter Real k31(unit="1") = 33.333 "k31 constant";
    parameter Real K_duo(unit="ug.h-1") = 177.34
      "Saturation parameter duodenal uptake, corresponds to k11, units??";
    // !!! Problem u rovnice (14) - doplnit prvni term na prave strane
    parameter Real v_duo(unit="h-1") = 0.70485
      "Duodenal uptake rate from blood (k6), 0.6/1.3";
    parameter Real u_duo(unit="h-1") = 0.88356
      "Duodenal export rate (k2), 0.66/1.38";
    parameter Real v_res(unit="h-1") = 6.3206
      "Other organs uptake rate (k34), 5.4/10.0";
    parameter Real u_res(unit="h-1") = 0.017143
      "Other organs export rate (k33), 0.014/0.030";
    parameter Real v_duo_lost(unit="h-1") = 0.091919
      "Iron lost rate duodenum (k48), 0.001/0.320";
    parameter Real u_res_lost(unit="h-1") = 0.0033401
      "Iron lost rate rest (k9), 0.002/0.004";
    parameter Real Fe_res_max(unit="ug") = 510.68
      "Limit value, iron lost rest (k49), 510/947";

    //################################################
    // EQUATIONS
    //################################################

  equation

    // Hepcidin expression
    //################################################
    der(hep) = r_hep_syn - r_hep_deg;
    // eq 1.
    r_hep_syn = Promoter(Bmp6, Il6);
    // calling set of in-house built functions
    r_hep_deg = k_hep_deg*hep;

    der(Bmp6) = r_Bmp6_syn - r_Bmp6_deg;
    // eq 2.
    r_Bmp6_syn = v_Bmp6_max*(Fe_liv/(K_Bmp6 + Fe_liv))*min(Fe_ser, Tf);
    r_Bmp6_deg = k_Bmp6_deg*Bmp6;

    der(LPS) = -k_LPS_deg*LPS;
    // eq 3.

    der(Il6mRNA) = r_Il6mRNA_syn - r_Il6mRNA_deg;
    // eq 4.
    r_Il6mRNA_syn = LPS/(LPS + K_Il6mRNA);
    r_Il6mRNA_deg = k_Il6mRNA_deg*Il6mRNA;

    der(Il6) = r_Il6_syn - r_Il6_deg;
    // eq 5.
    r_Il6_syn = k_Il6_syn*(Il6mRNA^4);
    r_Il6_deg = k_Il6_deg*Il6;

    // Ferroportin regulation
    //################################################
    der(Fpn_liv_mRNA) = 1/(1 + (K_liv_1*Il6)/(K_2 + Il6)) - k_FpnmRNA_deg*
      Fpn_liv_mRNA;
    // eq 6.1
    der(Fpn_spl_mRNA) = 1/(1 + (K_spl_1*Il6)/(K_2 + Il6)) - k_FpnmRNA_deg*
      Fpn_spl_mRNA;
    // eq 6.2
    der(Fpn_duo_mRNA) = 1/(1 + (K_duo_1*Il6)/(K_2 + Il6)) - k_FpnmRNA_deg*
      Fpn_duo_mRNA;
    // eq 6.3
    der(Fpn_res_mRNA) = 1/(1 + (K_res_1*Il6)/(K_2 + Il6)) - k_FpnmRNA_deg*
      Fpn_res_mRNA;
    // eq 6.4

    der(Fpn_liv) = k_Fpnliv_syn*(1 + k_liv_1*Fe_liv)*Fpn_liv_mRNA -
      k_Fpnliv_deg*(1 + k_liv_2*hep)*Fpn_liv;
    // eq 7.1
    der(Fpn_spl) = k_Fpnspl_syn*(1 + k_spl_1*Fe_spl)*Fpn_spl_mRNA -
      k_Fpnspl_deg*(1 + k_spl_2*hep)*Fpn_spl;
    // eq 7.2
    der(Fpn_duo) = k_Fpnduo_syn*(1 + k_duo_1*Fe_duo)*Fpn_duo_mRNA -
      k_Fpnduo_deg*(1 + k_duo_2*hep)*Fpn_duo;
    // eq 7.3
    der(Fpn_res) = k_Fpnres_syn*(1 + k_res_1*Fe_res)*Fpn_res_mRNA -
      k_Fpnres_deg*(1 + k_res_2*hep)*Fpn_res;
    // eq 7.4

    // Dynamics of the iron pools
    //################################################
    der(Fe_liv) = Fe_ser_liv - Fe_liv_ser;
    // eq 8.
    Fe_ser_liv = v_liv_1*min(Fe_ser, th) + v_liv_2*max(Fe_ser - th, 0);
    // eq 9.
    Fe_liv_ser = u_liv*min(Fe_liv, Fe_liv_max)*Fpn_liv;
    // eq 10.

    der(Fe_spl) = v_spl_1*Fe_RBC + v_spl_2*Fe_bm - u_spl*min(Fe_spl, Fe_spl_max)
      *Fpn_spl;
    // eq 11.

    der(Fe_bm) = v_bm*Fe_ser - (v_RBC + v_spl_2)*Fe_bm;
    //eq 12.

    der(Fe_RBC) = v_RBC*Fe_bm - v_spl_1*Fe_RBC;
    //eq 13.

    der(Fe_duo) = k14*(k31/(1 + k31/K_duo))*1/Fe_duo + v_duo*Fe_ser - u_duo*
      Fe_duo*Fpn_duo - v_duo_lost*Fe_duo;
    // eq 14. !!! First term issue

    der(Fe_res) = v_res*Fe_ser - u_res*Fe_res*Fpn_res - u_res_lost*min(Fe_res,
      Fe_res_max);
    // eq 15.

    der(Fe_ser) = Fe_input - Fe_output;
    // eq 16.
    Fe_input = u_liv*min(Fe_liv, Fe_liv_max)*Fpn_liv + u_spl*min(Fe_spl,
      Fe_spl_max)*Fpn_spl + u_duo*Fe_duo*Fpn_duo + u_res*Fe_res*Fpn_res;
    Fe_output = v_liv_1*min(Fe_ser, th) + v_liv_2*max(Fe_ser - th, 0) + (v_bm
       + v_duo + v_res)*Fe_ser;

    annotation (
      Icon(coordinateSystem(preserveAspectRatio=false)),
      Diagram(coordinateSystem(preserveAspectRatio=false)),
      experiment(
        StopTime=100,
        __Dymola_NumberOfIntervals=1500,
        Tolerance=1e-009));
  end FeMetabolismModel;

  function facB1 "facB1 function (id f9)"
    input Real facB1_i1;
    input Real facB1_i2;
    output Real facB1_o1 "result";
  algorithm
    facB1_o1 := (pSMAD(pSMAD_i1=facB1_i1, pSMAD_i2=facB1_i2)^1.7807)/0.4391;
  end facB1;

  function facB2 "facB2 function (id f10)"
    input Real facB2_i1;
    input Real facB2_i2;
    output Real facB2_o1 "result";
  algorithm
    facB2_o1 := (pSMAD(pSMAD_i1=facB2_i1, pSMAD_i2=facB2_i2)^1.7807)/16.8738;
  end facB2;

  function facSig "facSig function (id f5)"
    input Real facSig_i1;
    input Real facSig_i2;
    output Real facSig_o1 "result";
  algorithm
    facSig_o1 := (0.1285*2.852*hill(
        hill_i1=facSig_i1,
        hill_i2=1.0242,
        hill_i3=7.7388) - (1 + 0.4135*(0.0583 + 1.949*hill(
        hill_i1=facSig_i2,
        hill_i2=1.4481,
        hill_i3=140.244))))/(1 + 0.4135*0.0583);
  end facSig;

  function facSig1 "facSig1 function (id f6)"
    input Real facSig1_i1;
    output Real facSig1_o1 "result";
  algorithm
    facSig1_o1 := 0.1285*2.852*hill(
        hill_i1=facSig1_i1,
        hill_i2=1.0242,
        hill_i3=7.7388)/(1 + 0.4135*0.0583);
  end facSig1;

  function facST "facST function (id f11)"
    input Real facST_i1;
    input Real facST_i2;
    output Real facST_o1 "result";
  algorithm
    facST_o1 := pSTAT(pSTAT_i1=facST_i1, pSTAT_i2=facST_i2)/206.3988;
  end facST;

  function Freg "Freg function (id f14)"
    input Real Freg_i1;
    input Real Freg_i2;
    output Real Freg_o1 "result";
  algorithm
    Freg_o1 := Freg1(Freg1_i1=Freg_i1, Freg1_i2=Freg_i2)/Freg2(Freg2_i1=Freg_i1,
      Freg2_i2=Freg_i2);
  end Freg;

  function Freg1 "Freg1 function (id f12)"
    input Real Freg1_i1;
    input Real Freg1_i2;
    output Real Freg1_o1 "result";
  algorithm
    Freg1_o1 := 1 + FeMetabolism.facB1(facB1_i1=Freg1_i1, facB1_i2=Freg1_i2)*
      537.649 + FeMetabolism.facB2(facB2_i1=Freg1_i1, facB2_i2=Freg1_i2)*4972.6
       + FeMetabolism.facST(facST_i1=Freg1_i1, facST_i2=Freg1_i2)*584.75 +
      FeMetabolism.facB1(facB1_i1=Freg1_i1, facB1_i2=Freg1_i2)*
      FeMetabolism.facB2(facB2_i1=Freg1_i1, facB2_i2=Freg1_i2)*537.649*4972.6
       + FeMetabolism.facB1(facB1_i1=Freg1_i1, facB1_i2=Freg1_i2)*
      FeMetabolism.facST(facST_i1=Freg1_i1, facST_i2=Freg1_i2)*537.649*584.75*
      5.3869 + FeMetabolism.facB2(facB2_i1=Freg1_i1, facB2_i2=Freg1_i2)*
      FeMetabolism.facST(facST_i1=Freg1_i1, facST_i2=Freg1_i2)*4972.6*584.75 +
      FeMetabolism.facB1(facB1_i1=Freg1_i1, facB1_i2=Freg1_i2)*
      FeMetabolism.facB2(facB2_i1=Freg1_i1, facB2_i2=Freg1_i2)*
      FeMetabolism.facST(facST_i1=Freg1_i1, facST_i2=Freg1_i2)*537.649*4972.6*
      584.75*5.3869;
  end Freg1;

  function Freg2 "Freg2 function (id f13)"
    input Real Freg2_i1;
    input Real Freg2_i2;
    output Real Freg2_o1 "result";
  algorithm
    Freg2_o1 := 1 + FeMetabolism.facB1(facB1_i1=Freg2_i1, facB1_i2=Freg2_i2) +
      FeMetabolism.facB2(facB2_i1=Freg2_i1, facB2_i2=Freg2_i2) +
      FeMetabolism.facST(facST_i1=Freg2_i1, facST_i2=Freg2_i2) +
      FeMetabolism.facB1(facB1_i1=Freg2_i1, facB1_i2=Freg2_i2)*
      FeMetabolism.facB2(facB2_i1=Freg2_i1, facB2_i2=Freg2_i2) +
      FeMetabolism.facB1(facB1_i1=Freg2_i1, facB1_i2=Freg2_i2)*
      FeMetabolism.facST(facST_i1=Freg2_i1, facST_i2=Freg2_i2)*5.3869 +
      FeMetabolism.facB2(facB2_i1=Freg2_i1, facB2_i2=Freg2_i2)*
      FeMetabolism.facST(facST_i1=Freg2_i1, facST_i2=Freg2_i2) +
      FeMetabolism.facB1(facB1_i1=Freg2_i1, facB1_i2=Freg2_i2)*
      FeMetabolism.facB2(facB2_i1=Freg2_i1, facB2_i2=Freg2_i2)*
      FeMetabolism.facST(facST_i1=Freg2_i1, facST_i2=Freg2_i2)*5.3869;
  end Freg2;

  function hill "Hill function (id f4)"
    input Real hill_i1;
    input Real hill_i2;
    input Real hill_i3;
    output Real hill_o1 "result";
    constant Real eps=1e-9;
    Real i1_lim;
    Real i3_lim;
  algorithm
    if hill_i1 < eps then
      i1_lim := eps;
    else
      i1_lim := hill_i1;
    end if;

    if hill_i3 < eps then
      i3_lim := eps;
    else
      i3_lim := hill_i3;
    end if;

    hill_o1 := (i1_lim^hill_i2)/(i1_lim^hill_i2 + i3_lim^hill_i2);
  end hill;

  function Promoter "Promoter occupancy function (id f15)"
    input Real Promoter_i1;
    input Real Promoter_i2;
    output Real Promoter_o1 "result";
  algorithm
    Promoter_o1 := FeMetabolism.Freg(Freg_i1=Promoter_i1, Freg_i2=Promoter_i2)/
      (FeMetabolism.Freg(Freg_i1=Promoter_i1, Freg_i2=Promoter_i2) + 6804.7);
  end Promoter;

  function pSMAD "pSMAD function (id f8)"
    input Real pSMAD_i1;
    input Real pSMAD_i2;
    output Real pSMAD_o1 "result";
  algorithm
    pSMAD_o1 := 0.0583 + 1.949*FeMetabolism.hill(
        hill_i1=pSMAD_i2,
        hill_i2=1.4481,
        hill_i3=140.244)/(1 + 0.1285*pSTAT(pSTAT_i1=pSMAD_i1, pSTAT_i2=pSMAD_i2));
  end pSMAD;

  function pSTAT "pSTAT function (id f7)"
    input Real pSTAT_i1;
    input Real pSTAT_i2;
    output Real pSTAT_o1 "result";
  algorithm
    pSTAT_o1 := 0.5/0.1285*(FeMetabolism.facSig(facSig_i1=pSTAT_i1, facSig_i2=
      pSTAT_i2) + sqrt(FeMetabolism.facSig(facSig_i1=pSTAT_i1, facSig_i2=
      pSTAT_i2)^2 + 4*FeMetabolism.facSig1(facSig1_i1=pSTAT_i1)));
  end pSTAT;
end FeMetabolism;
