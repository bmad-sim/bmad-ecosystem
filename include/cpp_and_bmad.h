#ifndef CPP_AND_BMAD

#include <string>
#include <valarray>
#include "bmad_parameters.h"

using namespace std;

class C_coord;
class C_taylor_term;
class C_sr1_wake;
class C_sr2_wake;
class C_lr_wake;
class C_control;
class C_ele;
class C_wig_term;
class C_taylor;

typedef const double    Re;
typedef const int       Int;
typedef const char*     Char;
typedef const bool      Bool;
typedef const double*   ReArr;
typedef const int*      IntArr;

typedef valarray<double>                Real_Array;
typedef valarray<bool>                  Bool_Array;
typedef valarray<int>                   Int_Array;
typedef valarray<C_taylor_term>         C_taylor_term_array;
typedef valarray<C_wig_term>            C_wig_term_array;
typedef valarray<C_sr1_wake>            C_sr1_wake_array;
typedef valarray<C_sr2_wake>            C_sr2_wake_array;
typedef valarray<C_lr_wake>             C_lr_wake_array;
typedef valarray<C_control>             C_control_array;
typedef valarray<C_ele>                 C_ele_array;
typedef valarray<C_taylor>              C_taylor_array;
typedef valarray<Real_Array>            Real_Matrix;
typedef valarray<Bool_Array>            Bool_Matrix;

const Real_Array V2_array(double(0), 2);
const Real_Array V3_array(double(0), 3);
const Real_Array V6_array(double(0), 6);
const Real_Matrix M2_mat(V2_array, 2);
const Real_Matrix M3_mat(V3_array, 3);
const Real_Matrix M6_mat(V6_array, 6);
const double     V0(0);

//--------------------------------------------------------------------

bool is_all_true (Bool_Array v);
bool is_any_true (Bool_Array v);

template <class T> void operator<< (valarray<T>& arr, const T* ptr);
template <class T> void operator<< (valarray< valarray<T> >& mat, const T* ptr);
template <class T> void operator<< (valarray<T>& arr1, const valarray<T>& arr2);
template <class T> void operator<< (valarray< valarray<T> >& mat1, 
                              const valarray< valarray<T> >& mat2);

void matrix_to_array (const Real_Matrix& mat, double* arr);

//--------------------------------------------------------------------
// Coord 

class coord_struct {};
extern "C" void init_coord_struct_(coord_struct*&);

class C_coord {
public:
  Real_Array vec;   // size = 6

  C_coord(Re v[6]) : vec(v, 6) {}

  C_coord(double v0, double v1, double v2, double v3, double v4, double v5)
     {double v[] = {v0, v1, v2, v3, v4, v5}; vec = Real_Array(v, 6);}

  C_coord(Re v = 0) : vec(v, 6) {}

  C_coord(Int i) : vec(double(i), 6) {}

  C_coord(Real_Array v) : vec(v) {}

};    // End Class

extern "C" void coord_to_c_(coord_struct*, C_coord&);
extern "C" void coord_to_f_(C_coord&, coord_struct*);

bool operator== (const C_coord&, const C_coord&);

void operator>> (C_coord&, coord_struct*);
void operator>> (coord_struct*, C_coord&);

//--------------------------------------------------------------------
// orbit

class orbit_struct {};

class C_orbit {
public:
  valarray<C_coord> at;  // size = 6

  C_orbit(Int ix = 1) : at(C_coord(0), ix) {}

};    // End Class

extern "C" void orbit_to_c_(orbit_struct*, C_orbit&);
extern "C" void orbit_to_f_(C_orbit&, orbit_struct*);

bool operator== (const C_orbit&, const C_orbit&);

void operator>> (C_orbit&, orbit_struct*);
void operator>> (orbit_struct*, C_orbit&);

//--------------------------------------------------------------------
// Twiss 

class twiss_struct {}; 

class C_twiss {
public:
  double beta, alpha, gamma, phi, eta, etap;
  double eta_lab, etap_lab;
  double sigma;

  C_twiss(double b, double a, double g, double p, double e, double ep, 
                                         double el, double epl, double s) : 
    beta(b), alpha(a), gamma(g), phi(p), eta(e), etap(ep), eta_lab(el),
    etap_lab(epl), sigma(s) {}

  C_twiss(double z = 0) : 
    beta(z), alpha(z), gamma(z), phi(z), eta(z), 
    etap(z), eta_lab(z), etap_lab(z), sigma(z) {}

};    // End Class

extern "C" void twiss_to_c_(twiss_struct*, C_twiss&);
extern "C" void twiss_to_f_(const C_twiss&, twiss_struct*);

bool operator== (const C_twiss&, const C_twiss&);

void operator>> (C_twiss&, twiss_struct*);
void operator>> (twiss_struct*, C_twiss&);


//--------------------------------------------------------------------
// Floor_position 

class floor_position_struct {};

class C_floor_position {
public:
  double x, y, z;            // offset from origin
  double theta, phi, psi;    // angular orientation

  C_floor_position (double xx, double yy, double zz, double th, double ph, double ps) :
    x(xx), y(yy), z(zz), theta(th), phi(ph), psi(ps) {}

  C_floor_position (double z = 0) : 
    x(z), y(z), z(z), theta(z), phi(z), psi(z) {}

};    // End Class

extern "C" void floor_position_to_c_(floor_position_struct*, C_floor_position&);
extern "C" void floor_position_to_f_(C_floor_position&, floor_position_struct*);

bool operator== (const C_floor_position&, const C_floor_position&);

void operator>> (C_floor_position&, floor_position_struct*);
void operator>> (floor_position_struct*, C_floor_position&);

//--------------------------------------------------------------------
// Wig_term 

class wig_term_struct {};

class C_wig_term {
public:
  double coef;
  double kx, ky, kz;
  double phi_z;
  int type;      // HYPER_Y, HYPER_XY, OR HYPER_X

  C_wig_term (double c, double x, double y, double z, double p, int t) :
    coef(c), kx(x), ky(y), kz(z), phi_z(p), type(t) {}

  C_wig_term (double z = 0, int t = 1) : 
    coef(z), kx(z), ky(z), kz(z), phi_z(z), type(t) {}

};    // End Class

extern "C" void wig_term_to_c_(wig_term_struct*, C_wig_term&);
extern "C" void wig_term_to_f_(C_wig_term&, wig_term_struct*);

bool operator== (const C_wig_term&, const C_wig_term&);

void operator>> (C_wig_term&, wig_term_struct*);
void operator>> (wig_term_struct*, C_wig_term&);

//--------------------------------------------------------------------
// Taylor_term 

class taylor_term_struct {};

class C_taylor_term {
public:
  double coef;
  Int_Array exp;  // size = 6

  C_taylor_term (double c, int e[6]) : coef(c), exp(e, 6) {}

  C_taylor_term (double c, Int_Array e) : coef(c), exp(e) {}

  C_taylor_term (double c, int e0, int e1, int e2, int e3, int e4, int e5) :
    coef(c) {int e[] = {e0, e1, e2, e3, e4, e5}; exp = Int_Array(e, 6);}

  C_taylor_term (double c = 0) : coef(c), exp(0, 6) {}

};    // End Class

extern "C" void taylor_term_to_c_(taylor_term_struct*, C_taylor_term&);
extern "C" void taylor_term_to_f_(C_taylor_term&, taylor_term_struct*);

bool operator== (const C_taylor_term&, const C_taylor_term&);

void operator>> (C_taylor_term&, taylor_term_struct*);
void operator>> (taylor_term_struct*, C_taylor_term&);

//--------------------------------------------------------------------
// Taylor 

class taylor_struct {};

class C_taylor {
public:
  double ref;
  C_taylor_term_array term;  // size = variable

  C_taylor(int n_term = 0, double r = 0) : 
      ref(r), term (C_taylor_term(), n_term) {}

  C_taylor& operator= (const C_taylor&);

};    // End Class

extern "C" void taylor_to_c_(taylor_struct*, C_taylor&);
extern "C" void taylor_to_f_(C_taylor&, taylor_struct*);

bool operator== (const C_taylor&, const C_taylor&);

void operator>> (C_taylor&, taylor_struct*);
void operator>> (taylor_struct*, C_taylor&);

//--------------------------------------------------------------------
// SR1_wake 

class sr1_wake_struct {};

class C_sr1_wake {
public:
  double z;                 // Longitudinal distance
  double longitudinal;      // Longitudinal wake in V/C/m
  double transverse;        // Transverse wake in V/C/m^2

  C_sr1_wake (double zz, double lw, double tw) :
      z(zz), longitudinal(lw), transverse(tw) {}

  C_sr1_wake (double zz = 0) :
      z(zz), longitudinal(0), transverse(0) {}
};    // End Class

extern "C" void sr1_wake_to_c_(sr1_wake_struct*, C_sr1_wake&);
extern "C" void sr1_wake_to_f_(C_sr1_wake&, sr1_wake_struct*);

bool operator== (const C_sr1_wake&, const C_sr1_wake&);

void operator>> (C_sr1_wake&, sr1_wake_struct*);
void operator>> (sr1_wake_struct*, C_sr1_wake&);

//--------------------------------------------------------------------
// SR2_wake 

class sr2_wake_struct {};

class C_sr2_wake {
public:
  double amp;         // Amplitude
  double damp;        // damping factor
  double freq;        // Freq in Hz
  double phi;         // Phase in radians/2pi
  double norm_sin;    // non-skew sin-like component of the wake
  double norm_cos;    // non-skew cos-like component of the wake
  double skew_sin;    // skew sin-like component of the wake
  double skew_cos;    // skew cos-like component of the wake


  C_sr2_wake (double a, double d, double f, double p, double n_sin = 0, 
                  double n_cos = 0, double s_sin = 0, double s_cos = 0) :
      amp(a), damp(d), freq(f), phi(p), norm_sin(n_sin), 
      norm_cos(n_cos), skew_sin(s_sin), skew_cos(s_cos) {}

  C_sr2_wake (double a = 0) :
      amp(0), damp(0), freq(0), phi(0), norm_sin(0), norm_cos(0), 
      skew_sin(0), skew_cos(0) {}

};    // End Class

extern "C" void sr2_wake_to_c_(sr2_wake_struct*, C_sr2_wake&);
extern "C" void sr2_wake_to_f_(C_sr2_wake&, sr2_wake_struct*);

bool operator== (const C_sr2_wake&, const C_sr2_wake&);

void operator>> (C_sr2_wake&, sr2_wake_struct*);
void operator>> (sr2_wake_struct*, C_sr2_wake&);

//--------------------------------------------------------------------
// LR_wake  

class lr_wake_struct {};

class C_lr_wake {
public:
  double freq;       // frequency in Hz
  double freq_in;    // freq in input file. strength in V/C/m^2
  double R_over_Q;   // wake strength.
  double Q;          // Quality factor
  double angle;      // Polarization angle
  double norm_sin;
  double norm_cos;
  double skew_sin;
  double skew_cos;
  int m;             // Order number (1 = dipole, etc.)
  bool polarized;

  C_lr_wake (double f, double f_in, double rq, double q, double ang,
          double n_sin, double n_cos, double s_sin, double s_cos, 
          int mm, bool pol) :
      freq(f), freq_in(f_in), R_over_Q(rq), Q(q), angle(ang), norm_sin(n_sin),
      norm_cos(n_cos), skew_sin(s_sin), skew_cos(s_cos), m(mm), polarized(pol){}

  C_lr_wake (double f = 0) :
      freq(f), freq_in(0), R_over_Q(0), Q(0), angle(0), norm_sin(0), norm_cos(0),
      skew_sin(0), skew_cos(0), m(0), polarized(0){}
};    // End Class

extern "C" void lr_wake_to_c_(lr_wake_struct*, C_lr_wake&);
extern "C" void lr_wake_to_f_(C_lr_wake&, lr_wake_struct*);

bool operator== (const C_lr_wake&, const C_lr_wake&);

void operator>> (C_lr_wake&, lr_wake_struct*);
void operator>> (lr_wake_struct*, C_lr_wake&);

//--------------------------------------------------------------------
// Wake 

class wake_struct {};

class C_wake {
public:
  string sr_file;
  string lr_file;
  C_sr1_wake_array sr1;        // size = variable
  C_sr2_wake_array sr2_long;   // size = variable
  C_sr2_wake_array sr2_trans;  // size = variable
  C_lr_wake_array lr;          // size = variable
  double z_cut_sr;             // Cutoff between sr1 and sr2

  C_wake (const char* srf, const char* lrf, int n_sr1, int n_sr2_long, int n_sr2_trans, int n_lr) : 
      sr1(C_sr1_wake(), n_sr1), sr2_long(C_sr2_wake(), n_sr2_long), 
      sr2_trans(C_sr2_wake(), n_sr2_trans), lr(C_lr_wake(), n_lr),
      sr_file(string(srf, strlen(srf))),
      lr_file(string(lrf, strlen(lrf))) {}

  C_wake (string srf, string lrf, int n_sr1, int n_sr2_long, int n_sr2_trans, int n_lr) : 
      sr_file(srf), lr_file(lrf), sr1(C_sr1_wake(), n_sr1), sr2_long(C_sr2_wake(), n_sr2_long), 
      sr2_trans(C_sr2_wake(), n_sr2_trans), lr(C_lr_wake(), n_lr) {}

  C_wake () : sr_file(""), lr_file("") {}

  C_wake& operator= (const C_wake&);

};    // End Class

extern "C" void wake_to_c_(wake_struct*, C_wake&);
extern "C" void wake_to_f_(C_wake&, wake_struct*);

bool operator== (const C_wake&, const C_wake&);

void operator>> (C_wake&, wake_struct*);
void operator>> (wake_struct*, C_wake&);

//--------------------------------------------------------------------
// Control 

class control_struct {};

class C_control {
public:
  double coef;                // control coefficient
  int ix_lord;                // index to lord element
  int ix_slave;               // index to slave element
  int ix_attrib;              // index of attribute controlled

  C_control (double c, int il, int is, int ia) :
      coef(c), ix_lord(il), ix_slave(is), ix_attrib(ia) {}

  C_control () :
      coef(0), ix_lord(0), ix_slave(0), ix_attrib(0) {}
};    // End Class

extern "C" void control_to_c_(control_struct*, C_control&);
extern "C" void control_to_f_(C_control&, control_struct*);

bool operator== (const C_control&, const C_control&);

void operator>> (C_control&, control_struct*);
void operator>> (control_struct*, C_control&);

//--------------------------------------------------------------------
// Param 

class param_struct {};

class C_param {
public:
  double n_part;            // Particles/bunch (for BeamBeam elements).
  double total_length;      // total_length of ring
  double growth_rate;       // growth rate/turn if not stable
  Real_Matrix t1_with_RF;   // Full 1-turn matrix with RF on.
  Real_Matrix t1_no_RF;     // Full 1-turn matrix with RF off.
  int particle;             // positron$, electron$, etc.
  int ix_lost;              // Index of element particle was lost at.
  int end_lost_at;          // entrance_end$ or exit_end$
  int lattice_type;         // linear_lattice$, etc...
  int ixx;                  // Int for general use
  int ran_seed;             // Random number generator seed
  bool stable;              // is closed ring stable?
  bool aperture_limit_on;   // use apertures in tracking?
  bool lost;                // for use in tracking

  C_param () : n_part(0), total_length(0), growth_rate(0),
      t1_with_RF(V6_array, 6), t1_no_RF(V6_array, 6), 
      particle(0), ix_lost(0), end_lost_at(0), lattice_type(0), ixx(0),
      ran_seed(0), stable(1), aperture_limit_on(1), lost(0) {}

  C_param (double np, double tl, double gr, Real_Matrix t1w,
    Real_Matrix t1n, int pa, int il, int ela, int lt, int ix, 
    int r_seed, int st, int alo, int lo) :
        n_part(np), total_length(tl), growth_rate(gr),
        t1_with_RF(t1w), t1_no_RF(t1n), particle(pa), 
        ix_lost(il), end_lost_at(ela), lattice_type(lt), ixx(ix),
        ran_seed(r_seed), stable(st), aperture_limit_on(alo), lost(lo) {}
};    // End Class

extern "C" void param_to_c_(param_struct*, C_param&);
extern "C" void param_to_f_(C_param&, param_struct*);

bool operator== (const C_param&, const C_param&);

void operator>> (C_param&, param_struct*);
void operator>> (param_struct*, C_param&);

//--------------------------------------------------------------------
// Amode 

class amode_struct {};

class C_amode {
public:
  double emittance;        // Beam emittance
  double synch_int4;       // I4 Synchrotron integral
  double synch_int5;       // I5 Synchrotron integral
  double j_damp;           // damping partition number
  double alpha_damp;       // damping per turn
  double chrom;            // Chromaticity
  double tune;             // "Fractional" tune in radians

  C_amode () :
      emittance(0), synch_int4(0), synch_int5(0), j_damp(0), 
      alpha_damp(0), chrom(0), tune(0) {}

  C_amode (double em, double si4, double si5, double jd, 
                                  double ad, double ch, double tu) :
      emittance(em), synch_int4(si4), synch_int5(si5), j_damp(jd), 
      alpha_damp(ad), chrom(ch), tune(tu) {}

};    // End Class

extern "C" void amode_to_c_(amode_struct*, C_amode&);
extern "C" void amode_to_f_(C_amode&, amode_struct*);

bool operator== (const C_amode&, const C_amode&);

void operator>> (C_amode&, amode_struct*);
void operator>> (amode_struct*, C_amode&);

//--------------------------------------------------------------------
// linac_mode 

class linac_mode_struct {};

class C_linac_mode {
public:
  double i2_E4;        // Integral: g^2*  gamma^4
  double i3_E7;        // Integral: g^3*  gamma^7
  double i5a_E6;       // Integral: (g^3*  H_a)*  gamma^6
  double i5b_E6;       // Integral: (g^3*  H_b)*  gamma^6
  double sig_E1;       // Energy spread after 1 pass (eV)
  double emittance_a;  // a mode emittance at end of linac
  double emittance_b;  // b mode emittance at end of linac

  C_linac_mode () :
      i2_E4(0), i3_E7(0), i5a_E6(0), i5b_E6(0), sig_E1(0),
      emittance_a(0), emittance_b(0) {}

  C_linac_mode (double i2, double i3, double i5a, double i5b, double sig,
                                                       double a, double b) :
      i2_E4(i2), i3_E7(i3), i5a_E6(i5a), i5b_E6(i5b), sig_E1(sig),
      emittance_a(a), emittance_b(b) {}

};    // End Class

extern "C" void linac_mode_to_c_(linac_mode_struct*, C_linac_mode&);
extern "C" void linac_mode_to_f_(C_linac_mode&, linac_mode_struct*);

bool operator== (const C_linac_mode&, const C_linac_mode&);

void operator>> (C_linac_mode&, linac_mode_struct*);
void operator>> (linac_mode_struct*, C_linac_mode&);

//--------------------------------------------------------------------
// Modes 

class modes_struct {};

class C_modes {
public:
  double synch_int1;    // Synchrotron integrals I1
  double synch_int2;    // Synchrotron integrals I2
  double synch_int3;    // Synchrotron integrals I3
  double sigE_E;        // SigmaE/E
  double sig_z;         // Sigma_Z
  double e_loss;        // Energy loss / turn (eV)
  C_amode  a, b, z;
  C_linac_mode lin;

  C_modes () :
      synch_int1(0), synch_int2(0), synch_int3(0), sigE_E(0), 
      sig_z(0), e_loss(0), a(), b(), z(), lin() {}

  C_modes (double i1, double i2, double i3, double se, double sz, 
          double el, C_amode aa, C_amode bb, C_amode zz, C_linac_mode l) :
      synch_int1(i1), synch_int2(i2), synch_int3(i3), sigE_E(se), 
      sig_z(sz), e_loss(el), a(aa), b(bb), z(zz), lin(l) {}

};    // End Class

extern "C" void modes_to_c_(modes_struct*, C_modes&);
extern "C" void modes_to_f_(C_modes&, modes_struct*);

bool operator== (const C_modes&, const C_modes&);

void operator>> (C_modes&, modes_struct*);
void operator>> (modes_struct*, C_modes&);

//--------------------------------------------------------------------
// Bmad_com 

class bmad_com_struct {};
class C_bmad_com;

extern "C" void bmad_com_to_c_(C_bmad_com&);
extern "C" void bmad_com_to_f_(C_bmad_com&);

class C_bmad_com {
public:
  Real_Array d_orb;              // for the make_mat6_tracking routine
  double max_aperture_limit;   
  double grad_loss_sr_wake;                 // Internal var for LCavities.
  double rel_tollerance; 
  double abs_tollerance; 
  int taylor_order;              // 3rd order is default
  int default_integ_order;       // PTC integration order
  int default_num_steps;         // Number integration steps
  bool canonical_coords;         // Use (x, px) [not (x, x')]
  bool use_liar_lcavity;         // Liar like tracking?
  bool sr_wakes_on;              // Short range wakefields?
  bool lr_wakes_on;              // Long range wakefields
  bool mat6_track_symmetric;     // symmetric offsets
  bool auto_bookkeeper;          // Automatic bookkeeping when elements change?

  C_bmad_com () : d_orb(double(0), 6) {bmad_com_to_c_(*this);}

  C_bmad_com (ReArr orb, double max_ap, double kl, double rel_t,
                        double abs_t, int to, int io, int steps, int cc,
                        int liar, int sr, int lr, int sym, int a_book) :
      d_orb(orb, 6), max_aperture_limit(max_ap), grad_loss_sr_wake(kl), 
      rel_tollerance(rel_t), abs_tollerance(abs_t), taylor_order(to), 
      default_integ_order(io), default_num_steps(steps), canonical_coords(cc), 
      use_liar_lcavity(liar), sr_wakes_on(sr), lr_wakes_on(lr), 
      mat6_track_symmetric(sym), auto_bookkeeper(a_book) {}

  C_bmad_com (Real_Array orb, double max_ap, double kl, double rel_t,
                        double abs_t, int to, int io, int steps, int cc,
                        int liar, int sr, int lr, int sym, int a_book) :
      d_orb(orb), max_aperture_limit(max_ap), grad_loss_sr_wake(kl), 
      rel_tollerance(rel_t), abs_tollerance(abs_t), taylor_order(to), 
      default_integ_order(io), default_num_steps(steps), canonical_coords(cc), 
      use_liar_lcavity(liar), sr_wakes_on(sr), lr_wakes_on(lr), 
      mat6_track_symmetric(sym), auto_bookkeeper(a_book) {}

};    // End Class

bool operator== (const C_bmad_com&, const C_bmad_com&);

//--------------------------------------------------------------------
// EM_field 

class em_field_struct {};

class C_em_field {
public:
  Real_Array E;       // electric field, size = 3
  Real_Array B;       // magnetic field, size = 3
  Real_Array kick;    // kick, size = 3
  Real_Matrix dE;     // electric field gradient, size = 3x3
  Real_Matrix dB;     // magnetic field gradient, size = 3x3
  Real_Matrix dkick;  // kick gradiant, size = 3x3
  int type;           // kick_field$ or em_field$

  C_em_field () :
    E(V0, 3), B(V0, 3), kick(V0, 3), dE(V3_array, 3), dB(V3_array, 3),
    dkick(V3_array, 3), type(0) {}

  C_em_field (Real_Array e, Real_Array b, Real_Array k, 
                              Real_Matrix ee, Real_Matrix bb, Real_Matrix dk, int tp) :
      E(e), B(b), kick(k), dE(ee), dB(bb), dkick(dk), type(tp) {}

};    // End Class

extern "C" void em_field_to_c_(em_field_struct*, C_em_field&);
extern "C" void em_field_to_f_(C_em_field&, em_field_struct*);

bool operator== (const C_em_field&, const C_em_field&);

void operator>> (C_em_field&, em_field_struct*);
void operator>> (em_field_struct*, C_em_field&);

//--------------------------------------------------------------------
// Ele 

class ele_struct {};

class C_ele {

public:
  string name;                  // name of element
  string type;                  // type name
  string alias;                 // Another name
  string attribute_name;        // Used by overlays
  C_twiss  x, y, z;             // Twiss parameters at end of element
  C_floor_position floor;       // Global floor position at end of ele.
  Real_Array value;             // attribute values. size = N_ATTRIB_MAXX
  Real_Array gen0;              // constant part of the genfield map. size = 6
  Real_Array vec0;              // 0th order transport vector. size = 6
  Real_Matrix mat6;             // 1st order transport matrix. size = 6x6
  Real_Matrix c_mat;            // 2x2 C coupling matrix. size = 2x2
  double gamma_c;               // gamma associated with C matrix
  double s;                     // longitudinal position at the end
  Real_Matrix r;                // For general use. Not used by Bmad.
  Real_Array a;                 // multipole
  Real_Array b;                 // multipoles
  Real_Array const_arr;         // Working constants.
  string descrip;               // For general use
  void* gen_field;              // Pointer to a PTC genfield
  C_taylor_array taylor;        // Taylor terms
  C_wig_term_array wig_term;    // Wiggler Coefs
  C_wake wake;                  // Wakefields
  int key;                      // key value
  int sub_key;                  // For wigglers: map_type$, periodic_type$
  int control_type;             // SUPER_SLAVE$, OVERLAY_LORD$, etc.
  int ix_value;                 // Pointer for attribute to control
  int n_slave;                  // Number of slaves
  int ix1_slave;                // Start index for slave elements
  int ix2_slave;                // Stop  index for slave elements
  int n_lord;                   // Number of lords
  int ic1_lord;                 // Start index for lord elements
  int ic2_lord;                 // Stop  index for lord elements
  int ix_pointer;               // For general use. Not used by Bmad.
  int ixx;                      // Index for Bmad internal use
  int ix_ele;                   // Index in ring%ele_(:) array
  int mat6_calc_method;         // bmad_standard$, taylor$, etc.
  int tracking_method;          // bmad_standard$, taylor$, etc.
  int field_calc;               // Used with Boris, Runge-Kutta integrators.
  int num_steps;                // number of slices for DA_maps
  int integrator_order;         // For Etiennes' PTC: 2, 4, or 6.
  int ptc_kind;                 // For setting the ptc kind type.
  int taylor_order;             // Order of the taylor series.
  int aperture_at;              // Where aperture is applied. exit_end$, ...
  bool symplectify;             // Symplectify mat6 matrices.
  bool mode_flip;               // Have the normal modes traded places?
  bool multipoles_on;           // For turning multipoles on/off
  bool exact_rad_int_calc;      // Exact radiation integral calculation?
  bool field_master;            // Calculate strength from the field value?
  bool is_on;                   // For turning element on/off.
  bool internal_logic;          // For Bmad internal use only.
  bool logic;                   // For general use. Not used by Bmad.
  bool on_an_i_beam;            // Have an I_Beam overlay_lord?

  C_ele () : value(double(0), Bmad::N_ATTRIB_MAXX+1), taylor(C_taylor(0), 6),
    mat6(M6_mat), c_mat(M2_mat), vec0(V6_array), gen0(V6_array) {}

  C_ele& operator= (const C_ele&);

};    // End Class

extern "C" void ele_to_c_(ele_struct*, C_ele&);
extern "C" void ele_to_f_(C_ele&, ele_struct*);

bool operator== (const C_ele&, const C_ele&);

void operator>> (C_ele&, ele_struct*);
void operator>> (ele_struct*, C_ele&);

//--------------------------------------------------------------------
// Mode_info 

class mode_info_struct {};

class C_mode_info {
public:
  double tune;      // "fractional" tune in radians: 0 < tune < 2pi
  double emit;      // Emittance
  double chrom;     // Chromaticity

  C_mode_info () : tune(0), emit(0), chrom(0) {}
  C_mode_info (Re t, Re e, Re c) : tune(t), emit(e), chrom(c) {}

};    // End Class

extern "C" void mode_info_to_c_(mode_info_struct*, C_mode_info&);
extern "C" void mode_info_to_f_(C_mode_info&, mode_info_struct*);

bool operator== (const C_mode_info&, const C_mode_info&);

void operator>> (C_mode_info&, mode_info_struct*);
void operator>> (mode_info_struct*, C_mode_info&);

//--------------------------------------------------------------------
// Ring 

class ring_struct {};

class C_ring {
public:
  string name;               // Name of ring given by USE statement
  string lattice;            // Lattice
  string input_file_name;    // Name of the lattice input file
  string title;              // General title
  C_mode_info x, y, z;       // Tunes, etc.
  C_param param;             // Parameters
  int version;               // Version number
  int n_ele_use;             // Number of regular ring elements
  int n_ele_max;             // Index of last element used
  int n_control_max;         // Last index used in CONTROL_ array
  int n_ic_max;              // Last index used in IC_ array
  int input_taylor_order;    // As set in the input file
  C_ele ele_init;            // For use by any program
  C_ele_array ele;           // Array of elements
  C_control_array control;   // control list
  Int_Array ic;              // index to %control_(:)

  C_ring () {}

  C_ring& operator= (const C_ring&);

};    // End Class

extern "C" void ring_to_c_(ring_struct*, C_ring&);
extern "C" void ring_to_f_(C_ring&, ring_struct*);

bool operator== (const C_ring&, const C_ring&);

void operator>> (C_ring&, ring_struct*);
void operator>> (ring_struct*, C_ring&);

#define CPP_AND_BMAD
#endif
