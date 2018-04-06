within ;
model FeRegulation "Basic model of Fe regulation"
  // Main variables in differential equations

  Real NTBI(
    start=0e-6,
    fixed=true,
    unit="mol.l-1") "Non-transferin-bound iron";
  Real FeTf(
    start=0e-6,
    fixed=true,
    unit="mol.l-1") "Diferric transferrin";
  Real Hepatocyte(
    start=1.52e-6,
    fixed=true,
    unit="mol") "Iron in hepatocytes";
  Real Macrophage(
    start=2.63e-6,
    fixed=true,
    unit="mol") "Iron in macrophages";
  Real RBC(
    start=23.0e-6,
    fixed=true,
    unit="mol") "Iron in red blood cell (RBC) pool";

  // Hepcidine regulated rate constants, format k_{from}_{to}

  Real k_gut_NTBI(unit="mol.h-1") = b*Fpn_enter
    "Rate of iron transfer from gut to NTBI";
  Real k_hepat_NTBI(unit="atoms.cell-1.h-1") = c*Fpn_hepat
    "Rate of iron transfer from hepatocytes to NTBI";
  Real k_macro_NTBI(unit="mol.h-1") = d*Fpn_macro
    "Rate of iron transfer from macrophages to NTBI";

  Real Fpn_enter(unit="1") "Fpn level in enterocytes";
  Real Fpn_hepat(unit="1") "Fpn level in hepatocytes";
  Real Fpn_macro(unit="1") "Fpn level in macrophages";

  Real Fpn_enter_conc2(unit="1") "Fpn level in enterocytes for [TfR2] = conc2";
  Real Fpn_hepat_conc2(unit="1") "Fpn level in hepatocytes for [TfR2] = conc2";
  Real Fpn_macro_conc2(unit="1") "Fpn level in macrophages for [TfR2] = conc2";

  // Other time-dependent variables

  Real k_FeTf_hepat(unit="atoms.cell-1.h-1") = a*TfR2
    "Rate of iron transfer from FeTf to hepatocytes";
  Real k_NTBI_hepat(unit="mol.cell-1.h-1") = Vmax_NTBI_hepat*NTBI/(
    Km_NTBI_hepat + NTBI)
    "Rate of iron transfer from NTBI to hepatocytes, depends on [NTBI], Vmax_NTBI_hepat, Km_NTBI_hepat";
  Real TfR2(unit="1") = TfR2_max*(1 - exp(-time/tau))
    "Transferrin receptor type 2 level, value from interval (0,1), unitless TODO solution of DAE instead of DAE?";
  Real TfR2_max(unit="1") = FeTf/(Km_FeTf_hepat + FeTf)
    "Transferrin receptor type 2 maximum level, depends on [FeTf] and Km_FeTf_hepat";

  // time-independent rate constants, format k_{from}_{to}
  parameter Real k_NTBI_gut(unit="mol.h-1") = 0.266e-9
    "Rate of iron transfer from NTBI to gut, taken from [1]";
  parameter Real k_NTBI_FeTf(unit="l.mol-1.h-1") = 1e10
    "Rate of iron transfer from NTBI to FeTf, taken from [1]";
  parameter Real k_FeTf_RBC(unit="mol.h-1") = 36.0e-9
    "Rate of iron transfer from FeTf to RBC, value taken from [4]";
  parameter Real k_RBC_macro(unit="h-1") = 1.6e-3
    "Rate of iron transfer from RBC to macrophages, value taken from [4]";

  // Other scecies-dependent parameters

  parameter Real Vmax_NTBI_hepat(unit="mol.cell-1.h-1") = 2.49e-15
    "Michaelis-Mentan maximum rate Vmax for NTBI->hepat iron transfer, value taken from [1]";
  parameter Real Km_NTBI_hepat(unit="mol.l-1") = 7.2e-6
    "Michaelis-Menten constant Km for NTBI->hepat iron transfer, value taken from [1]";
  parameter Real Km_FeTf_hepat(unit="mol.l-1") = 2.5e-6
    "Michaelis-Menten constant Km for FeTf->hepat iron transfer, value taken from [1]";
  parameter Real a(unit="atoms.cell-1.h-1") = 2.5e7
    "Proportionality constant relating TfR2 to k_FeTf_hepat, value taken from [1,2]";
  parameter Real tau(unit="h") = 1.4
    "Time constant for approach of TfR2 to TfR2_max, value taken from [1,2]";
  parameter Real conc1(unit="mol.l-1") = 2.5e-6
    "FeTf concentration for TfR2_max(conc1) calculation";
  parameter Real conc2(unit="mol.l-1") = 5.0e-6
    "TfR2 concentration for Fpn(conc2) calculation";
  parameter Real TfR2_max_conc1(unit="1") = conc1/(Km_FeTf_hepat + conc1)
    "Transferrin receptor type 2 maximum level at [FeTf] = conc1";
  parameter Real TfR2_max_conc2(unit="1") = conc2/(Km_FeTf_hepat + conc2)
    "Transferrin receptor type 2 maximum level at [FeTf] = conc2";
  parameter Real Tf_max(unit="mol.l-1") = 50.0e-6
    "Total amount of Tf in circulation, taken from [3]";
  parameter Real n_hepat(unit="cell") = 1e8
    "Hepatocyte cell number, 1e8 for mouse, taken from [5]";
  parameter Real decay_rate(unit="l.mol-1") = 37.465
    "Fitted exponential decay rate constant, hepcidin being downregulated by NTBI [1]";

  // Sensitivities pf Fpn to changes in...

  parameter Real alpha(unit="1") = 1
    "Sensitivity of Fpn to changes in TfR2 for enterocytes";
  parameter Real beta(unit="1") = 1
    "Sensitivity of Fpn to changes in TfR2 for hepatocytes";
  parameter Real gamma(unit="1") = 1
    "Sensitivity of Fpn to changes in TfR2 for macrophages";
  parameter Real delta(unit="1") = 1
    "Sensitivity of Fpn to changes in intracellular iron stores in hepatocytes";
  parameter Real epsilon(unit="1") = 1
    "Sensitivity of Fpn to changes in intracellular iron stores in macrophages";

  // Initial values of intracellular iron levels in...

  parameter Real Hepatocyte_0(unit="mol") = 1.52e-6
    "Initial value of hepatocyte intracellular iron level, value for mouse, taken from [6]";
  parameter Real Macrophage_0(unit="mol") = 2.63e-6
    "Initial value of macrophage intracellular iron level, value for mouse, taken from [7,8]";

  // Proportionality constants for iron transport in...

  Real b(unit="mol.h-1") = k_NTBI_gut/Fpn_enter_conc2
    "The proportionality constant for enterocyte iron transport";
  Real c(unit="mol.h-1") = a*TfR2_max_conc2/Fpn_hepat_conc2
    "The proportionality constant for hepatocyte iron transport";
  Real d(unit="mol.h-1") = k_FeTf_RBC/Fpn_macro_conc2
    "The proportionality constant for macrophage iron transport";

  // Constants
  constant Real V_plasma(unit="dm3") = 1.3e-3
    "Plasma volume, 1.3e-3 dm3 for mouse, taken from [4]";
  constant Real Na(unit="atoms.mol-1") = 6.022e23 "Avogadro constant";


equation

  // eq for NTBI
  der(NTBI) = (k_gut_NTBI - k_NTBI_gut)/V_plasma - 2*k_NTBI_FeTf*NTBI*(Tf_max
     - FeTf) + k_hepat_NTBI*n_hepat/(V_plasma*Na) - k_NTBI_hepat*(n_hepat/
    V_plasma) + k_macro_NTBI/V_plasma;

  // eq for FeTf
  der(FeTf) = k_NTBI_FeTf*NTBI*(Tf_max - FeTf) - k_FeTf_hepat*n_hepat/(2*
    V_plasma*Na) - k_FeTf_RBC/(2*V_plasma);

  // eq for Hepatocyte
  der(Hepatocyte) = (k_FeTf_hepat - k_hepat_NTBI)*(n_hepat/Na) + k_NTBI_hepat*
    n_hepat;

  // eq for RBC
  der(RBC) = k_FeTf_RBC - k_RBC_macro*RBC;

  // eq for Macrophage
  der(Macrophage) = k_RBC_macro*RBC - k_macro_NTBI;


  // !!!!! ONLY ONE LEVEL OF THEORY PERMITTED !!!!!!
  // eqs for Fpn_enter (3 levels of complexity) - !!! select one level and comment the other
  Fpn_enter = -alpha*(TfR2 - TfR2_max_conc1) + TfR2_max_conc1;
  // level 1: Fpn regulation by hepcidine
  Fpn_enter_conc2 = -alpha*(conc2 - TfR2_max_conc1) + TfR2_max_conc1;
  // level 1: Fpn regulation by hepcidine, [TfR2] = conc2

  //Fpn_enter = -alpha*(exp(-decay_rate*NTBI)*TfR2 - TfR2_max_conc1) + TfR2_max_conc1;                                               // level 3: Fpn regulation by NTBI
  //Fpn_enter_conc2 = -alpha*(exp(-decay_rate*NTBI)*conc2 - TfR2_max_conc1) + TfR2_max_conc1;                                        // level 3: Fpn regulation by NTBI, [TfR2] = conc2


  // eqs for Fpn_hepat (3 levels of complexity[ - !!! select one level and comment the other
  Fpn_hepat = -beta*(TfR2 - TfR2_max_conc1) + TfR2_max_conc1;
  // level 1: regulation by hepcidine
  Fpn_hepat_conc2 = -beta*(conc2 - TfR2_max_conc1) + TfR2_max_conc1;
  // level 1: regulation by hepcidine, [TfR2] = conc2

  //Fpn_hepat = ((Hepatocyte/Hepatocyte_0)^delta)*(-beta*(TfR2 - TfR2_max_conc1)  + TfR2_max_conc1);                             // level 2: regulation by intracellular iron stores
  //Fpn_hepat_conc2 = ((Hepatocyte/Hepatocyte_0)^delta)*(-beta*(conc2 - TfR2_max_conc1)  + TfR2_max_conc1);                      // level 2: regulation by intracellular iron stores, [TfR2] = conc2

  //Fpn_hepat = ((Hepatocyte/Hepatocyte_0)^delta)*(-beta*(exp(-decay_rate*NTBI)*TfR2 - TfR2_max_conc1)  + TfR2_max_conc1);           // level 3: Fpn regulation by NTBI
  //Fpn_hepat_conc2 = ((Hepatocyte/Hepatocyte_0)^delta)*(-beta*(exp(-decay_rate*NTBI)*conc2 - TfR2_max_conc1)  + TfR2_max_conc1);    // level 3: Fpn regulation by NTBI, [TfR2] = conc2


  // eqs for Fpn_macro (3 levels of complexity) - !!! select one level and comment the other
  Fpn_macro = -gamma*(TfR2 - TfR2_max_conc1) + TfR2_max_conc1;
  // level 1: regulation by hepcidine
  Fpn_macro_conc2 = -gamma*(conc2 - TfR2_max_conc1) + TfR2_max_conc1;
  // level 1: regulation by hepcidine, [TfR2] = conc2

  //Fpn_macro = ((Macrophage/Macrophage_0)^epsilon)*(-gamma*(TfR2 - TfR2_max_conc1) + TfR2_max_conc1);                           // level 2: regulation by intracellular iron stores
  //Fpn_macro_conc2 = ((Macrophage/Macrophage_0)^epsilon)*(-gamma*(conc2 - TfR2_max_conc1) + TfR2_max_conc1);                    // level 2: regulation by intracellular iron stores, [TfR2] = conc2

  //Fpn_macro = ((Macrophage/Macrophage_0)^epsilon)*(-gamma*(exp(-decay_rate*NTBI)*TfR2 - TfR2_max_conc1) + TfR2_max_conc1);         // level 3: Fpn regulation by NTBI
  //Fpn_macro_conc2 = ((Macrophage/Macrophage_0)^epsilon)*(-gamma*(exp(-decay_rate*NTBI)*conc2 - TfR2_max_conc1) + TfR2_max_conc1);  // level 3: Fpn regulation by NTBI, [TfR2] = conc2

  annotation (Documentation(info="<html>
<pre>REFERENCES
[1]&nbsp;Lao&nbsp;et&nbsp;al.,&nbsp;J.&nbsp;Theor.&nbsp;Biol.&nbsp;243,&nbsp;2006
[2]&nbsp;Cole&nbsp;et&nbsp;al.,&nbsp;Biochim.&nbsp;Biophys.&nbsp;Acta&nbsp;762,&nbsp;1983
[3] Johnson et al., Blood. 104, 2004
[4] Rivera et al., Blood. 106, 2005
[5] Molodykh et al., Bull. Exp. Biol. Med. 130, 2000
[6] Simpson et al., Blood. 101, 2003
[7] Makui et al., Blood. 106, 2005
[8] West et al., Am. J. Physiol. 275, 1998</pre>
</html>"));

end FeRegulation;
