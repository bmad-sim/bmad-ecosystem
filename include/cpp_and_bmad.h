#ifndef CPP_AND_BMAD

#include <string>
#include <string.h>
#include <valarray>
#include <complex>
#include "bmad_parameters.h"

using namespace std;

class C_coord;            // Needed for typedef...
class C_taylor_term;
class C_sr_table_wake;
class C_sr_mode_wake;
class C_lr_wake;
class C_control;
class C_ele;
class C_wig_term;
class C_taylor;
class C_branch;
class C_rf_field_mode;
class C_rf_field_mode_term;
class C_wall3d_vertex;
class C_wall3d_section;

typedef const double    Re;
typedef const int       Int;
typedef const char*     Char;
typedef const bool      Bool;
typedef const double*   ReArr;
typedef const int*      IntArr;

typedef complex<double>                 Complx;
typedef const complex<double>           CComplx;

typedef valarray<double>                Real_Array;
typedef valarray<Complx>                Complx_Array;
typedef valarray<bool>                  Bool_Array;
typedef valarray<int>                   Int_Array;

typedef valarray<Real_Array>            Real_Matrix;
typedef valarray<Bool_Array>            Bool_Matrix;

typedef valarray<C_taylor_term>         C_taylor_term_array;
typedef valarray<C_wig_term>            C_wig_term_array;
typedef valarray<C_sr_table_wake>       C_sr_table_wake_array;
typedef valarray<C_sr_mode_wake>        C_sr_mode_wake_array;
typedef valarray<C_lr_wake>             C_lr_wake_array;
typedef valarray<C_control>             C_control_array;
typedef valarray<C_branch>              C_branch_array;
typedef valarray<C_ele>                 C_ele_array;
typedef valarray<C_taylor>              C_taylor_array;
typedef valarray<C_rf_field_mode>       C_rf_field_mode_array;
typedef valarray<C_rf_field_mode_term>  C_rf_field_mode_term_array;
typedef valarray<C_wall3d_vertex>       C_wall3d_vertex_array;
typedef valarray<C_wall3d_section>      C_wall3d_section_array;

const Real_Array V2_array(double(0), 2);
const Real_Array V3_array(double(0), 3);
const Real_Array V6_array(double(0), 6);
const Real_Matrix M2_mat(V2_array, 2);
const Real_Matrix M3_mat(V3_array, 3);
const Real_Matrix M6_mat(V6_array, 6);
const double       V0(0);
const Complx       C0(0);

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
  Real_Array vec;      // size = 6
  double s;
  double t;
  Complx_Array spin;   // size = 2
  double e_field_x;
  double e_field_y;
  double phase_x;
  double phase_y;

  C_coord(Re v[6], double ss = 0, double tt = 0, CComplx spn[2] = 0, double field_x = 0, double field_y = 0, 
          double p_x = 0, double p_y = 0) : vec(v, 6), s(ss), t(tt), spin(spn, 2), 
          e_field_x(field_x), e_field_y(field_y), phase_x(p_x), phase_y(p_y) {}

  C_coord(double v0, double v1, double v2, double v3, double v4, double v5) :
     s(0), t(0), spin(2), e_field_x(0), e_field_y(0), phase_x(0), phase_y(0)
     {double v[] = {v0, v1, v2, v3, v4, v5}; vec = Real_Array(v, 6);}

  C_coord(Re v = 0) : vec(v, 6), s(0), t(0), spin(2), e_field_x(0), e_field_y(0), phase_x(0), phase_y(0) {}

  C_coord(Int i) : vec(double(i), 6), s(0), t(0), spin(2), e_field_x(0), e_field_y(0), phase_x(0), phase_y(0) {}

  C_coord(Real_Array v) : vec(v), s(0), t(0), spin(2), e_field_x(0), e_field_y(0), phase_x(0), phase_y(0) {}

};    // End Class

extern "C" void coord_to_c_(coord_struct*, C_coord&);
extern "C" void coord_to_f_(C_coord&, coord_struct*);

bool operator== (const C_coord&, const C_coord&);

void operator>> (C_coord&, coord_struct*);
void operator>> (coord_struct*, C_coord&);

//--------------------------------------------------------------------
// Twiss 

class twiss_struct {}; 

class C_twiss {
public:
  double beta, alpha, gamma, phi, eta, etap;
  double sigma, sigma_p, emit, norm_emit;

  C_twiss(double b, double a, double g, double p, double e, double ep, 
                                  double s, double sp, double em, double n_em) : 
    beta(b), alpha(a), gamma(g), phi(p), eta(e), etap(ep), 
    sigma(s), sigma_p(sp), emit(em), norm_emit(n_em) {}

  C_twiss(double z = 0) : 
    beta(z), alpha(z), gamma(z), phi(z), eta(z), 
    etap(z), sigma(z), sigma_p(z), emit(z), norm_emit(z) {}

};    // End Class

extern "C" void twiss_to_c_(twiss_struct*, C_twiss&);
extern "C" void twiss_to_f_(C_twiss&, twiss_struct*);

bool operator== (const C_twiss&, const C_twiss&);

void operator>> (C_twiss&, twiss_struct*);
void operator>> (twiss_struct*, C_twiss&);

//--------------------------------------------------------------------
// Xy_Disp 

class xy_disp_struct {}; 

class C_xy_disp {
public:
  double eta, etap;

  C_xy_disp(double e, double ep) : 
    eta(e), etap(ep) {}

  C_xy_disp(double z = 0) : 
    eta(z), etap(z) {}

};    // End Class

extern "C" void xy_disp_to_c_(xy_disp_struct*, C_xy_disp&);
extern "C" void xy_disp_to_f_(C_xy_disp&, xy_disp_struct*);

bool operator== (const C_xy_disp&, const C_xy_disp&);

void operator>> (C_xy_disp&, xy_disp_struct*);
void operator>> (xy_disp_struct*, C_xy_disp&);

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
  Int_Array expn;  // size = 6

  C_taylor_term (double c, int e[6]) : coef(c), expn(e, 6) {}

  C_taylor_term (double c, Int_Array e) : coef(c), expn(e) {}

  C_taylor_term (double c, int e0, int e1, int e2, int e3, int e4, int e5) :
    coef(c) {int e[] = {e0, e1, e2, e3, e4, e5}; expn = Int_Array(e, 6);}

  C_taylor_term (double c = 0) : coef(c), expn(0, 6) {}

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
// sr_table_wake 

class sr_table_wake_struct {};

class C_sr_table_wake {
public:
  double z;                 // Longitudinal distance
  double longitudinal;      // Longitudinal wake in V/C/m
  double transverse;        // Transverse wake in V/C/m^2

  C_sr_table_wake (double zz, double lw, double tw) :
      z(zz), longitudinal(lw), transverse(tw) {}

  C_sr_table_wake (double zz = 0) :
      z(zz), longitudinal(0), transverse(0) {}
};    // End Class

extern "C" void sr_table_wake_to_c_(sr_table_wake_struct*, C_sr_table_wake&);
extern "C" void sr_table_wake_to_f_(C_sr_table_wake&, sr_table_wake_struct*);

bool operator== (const C_sr_table_wake&, const C_sr_table_wake&);

void operator>> (C_sr_table_wake&, sr_table_wake_struct*);
void operator>> (sr_table_wake_struct*, C_sr_table_wake&);

//--------------------------------------------------------------------
// sr_mode_wake 

class sr_mode_wake_struct {};

class C_sr_mode_wake {
public:
  double amp;         // Amplitude
  double damp;        // damping factor
  double k;           // k factor
  double phi;         // Phase in radians
  double b_sin;       // non-skew sin-like component of the wake
  double b_cos;       // non-skew cos-like component of the wake
  double a_sin;       // skew sin-like component of the wake
  double a_cos;       // skew cos-like component of the wake


  C_sr_mode_wake (double a, double d, double kk, double p, double n_sin = 0, 
                  double n_cos = 0, double s_sin = 0, double s_cos = 0) :
      amp(a), damp(d), k(kk), phi(p), b_sin(n_sin), 
      b_cos(n_cos), a_sin(s_sin), a_cos(s_cos) {}

  C_sr_mode_wake (double a = 0) :
      amp(0), damp(0), k(0), phi(0), b_sin(0), b_cos(0), 
      a_sin(0), a_cos(0) {}

};    // End Class

extern "C" void sr_mode_wake_to_c_(sr_mode_wake_struct*, C_sr_mode_wake&);
extern "C" void sr_mode_wake_to_f_(C_sr_mode_wake&, sr_mode_wake_struct*);

bool operator== (const C_sr_mode_wake&, const C_sr_mode_wake&);

void operator>> (C_sr_mode_wake&, sr_mode_wake_struct*);
void operator>> (sr_mode_wake_struct*, C_sr_mode_wake&);

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
  double b_sin;
  double b_cos;
  double a_sin;
  double a_cos;
  double t_ref;
  int m;             // Order number (1 = dipole, etc.)
  bool polarized;

  C_lr_wake (double f, double f_in, double rq, double q, double ang,
          double n_sin, double n_cos, double s_sin, double s_cos, double t_ref,
          int mm, bool pol) :
      freq(f), freq_in(f_in), R_over_Q(rq), Q(q), angle(ang), b_sin(n_sin),
      b_cos(n_cos), a_sin(s_sin), a_cos(s_cos), t_ref(t_ref), m(mm), polarized(pol){}

  C_lr_wake (double f = 0) :
      freq(f), freq_in(0), R_over_Q(0), Q(0), angle(0), b_sin(0), b_cos(0),
      a_sin(0), a_cos(0), t_ref(0), m(0), polarized(0){}
};    // End Class

extern "C" void lr_wake_to_c_(lr_wake_struct*, C_lr_wake&);
extern "C" void lr_wake_to_f_(C_lr_wake&, lr_wake_struct*);

bool operator== (const C_lr_wake&, const C_lr_wake&);

void operator>> (C_lr_wake&, lr_wake_struct*);
void operator>> (lr_wake_struct*, C_lr_wake&);

//--------------------------------------------------------------------
// rf_wake 

class rf_wake_struct {};

class C_rf_wake {
public:
  string sr_file;
  string lr_file;
  C_sr_table_wake_array sr_table;        // size = variable
  C_sr_mode_wake_array sr_mode_long;   // size = variable
  C_sr_mode_wake_array sr_mode_trans;  // size = variable
  C_lr_wake_array lr;          // size = variable
  double z_sr_mode_max;             // Cutoff between sr_table and sr_mode

  C_rf_wake (const char* srf, const char* lrf, int n_sr_table, int n_sr_mode_long, int n_sr_mode_trans, int n_lr) : 
      sr_table(C_sr_table_wake(), n_sr_table), sr_mode_long(C_sr_mode_wake(), n_sr_mode_long), 
      sr_mode_trans(C_sr_mode_wake(), n_sr_mode_trans), lr(C_lr_wake(), n_lr),
      sr_file(string(srf, strlen(srf))),
      lr_file(string(lrf, strlen(lrf))) {}

  C_rf_wake (string srf, string lrf, int n_sr_table, int n_sr_mode_long, int n_sr_mode_trans, int n_lr) : 
      sr_file(srf), lr_file(lrf), sr_table(C_sr_table_wake(), n_sr_table), sr_mode_long(C_sr_mode_wake(), n_sr_mode_long), 
      sr_mode_trans(C_sr_mode_wake(), n_sr_mode_trans), lr(C_lr_wake(), n_lr) {}

  C_rf_wake () : sr_file(""), lr_file("") {}

  C_rf_wake& operator= (const C_rf_wake&);

};    // End Class

extern "C" void rf_wake_to_c_(rf_wake_struct*, C_rf_wake&);
extern "C" void rf_wake_to_f_(C_rf_wake&, rf_wake_struct*);

bool operator== (const C_rf_wake&, const C_rf_wake&);

void operator>> (C_rf_wake&, rf_wake_struct*);
void operator>> (rf_wake_struct*, C_rf_wake&);

//--------------------------------------------------------------------
// RF field mode

class rf_field_mode_term_struct {};

class C_rf_field_mode_term {
public:
  Complx e, b;

  C_rf_field_mode_term (int ee, int bb) : e(ee), b(bb) {}

  C_rf_field_mode_term () : e(0), b(0) {}

};    // End Class

extern "C" void rf_field_mode_term_to_c_(rf_field_mode_term_struct*, C_rf_field_mode_term&);
extern "C" void rf_field_mode_term_to_f_(C_rf_field_mode_term&, rf_field_mode_term_struct*);

bool operator== (const C_rf_field_mode_term&, const C_rf_field_mode_term&);

void operator>> (C_rf_field_mode_term&, rf_field_mode_term_struct*);
void operator>> (rf_field_mode_term_struct*, C_rf_field_mode_term&);

//--------------------------------------------------------------------
// RF field mode

class rf_field_mode_struct {};

class C_rf_field_mode {
public:
  int m;
  double freq;
  double f_damp;
  double theta_t0;
  double stored_energy;
  double phi_0;
  double dz;
  double field_scale;
  C_rf_field_mode_term_array term;

};    // End Class

extern "C" void rf_field_mode_to_c_(rf_field_mode_struct*, C_rf_field_mode&);
extern "C" void rf_field_mode_to_f_(C_rf_field_mode&, rf_field_mode_struct*);

bool operator== (const C_rf_field_mode&, const C_rf_field_mode&);

void operator>> (C_rf_field_mode&, rf_field_mode_struct*);
void operator>> (rf_field_mode_struct*, C_rf_field_mode&);

//--------------------------------------------------------------------
// RF_field

class rf_field_struct {};

class C_rf_field {
public:
  C_rf_field_mode_array mode;

};    // End Class

extern "C" void rf_field_to_c_(rf_field_struct*, C_rf_field&);
extern "C" void rf_field_to_f_(C_rf_field&, rf_field_struct*);

bool operator== (const C_rf_field&, const C_rf_field&);

void operator>> (C_rf_field&, rf_field_struct*);
void operator>> (rf_field_struct*, C_rf_field&);

//--------------------------------------------------------------------
// RF 

class rf_struct {};

class C_rf {
public:
  C_rf_wake wake;                 // RF Wakefields
  C_rf_field field;               // RF fields

};    // End Class

extern "C" void rf_to_c_(rf_struct*, C_rf&);
extern "C" void rf_to_f_(C_rf&, rf_struct*);

bool operator== (const C_rf&, const C_rf&);

void operator>> (C_rf&, rf_struct*);
void operator>> (rf_struct*, C_rf&);

//--------------------------------------------------------------------
// wall3d_vertex

class wall3d_vertex_struct {};

class C_wall3d_vertex {
public:
  double x, y;
  double radius_x;
  double radius_y;
  double tilt;
  double angle;
  double x0, y0;

};    // End Class

extern "C" void wall3d_vertex_to_c_(wall3d_vertex_struct*, C_wall3d_vertex&);
extern "C" void wall3d_vertex_to_f_(C_wall3d_vertex&, wall3d_vertex_struct*);

bool operator== (const C_wall3d_vertex&, const C_wall3d_vertex&);

void operator>> (C_wall3d_vertex&, wall3d_vertex_struct*);
void operator>> (wall3d_vertex_struct*, C_wall3d_vertex&);

//--------------------------------------------------------------------
// wall3d_section

class wall3d_section_struct {};

class C_wall3d_section {
public:
  int type;
  double s;
  Real_Array s_spline;
  int n_slice_spline;
  C_wall3d_vertex_array v;
  int n_vertex_input;

};    // End Class

extern "C" void wall3d_section_to_c_(wall3d_section_struct*, C_wall3d_section&);
extern "C" void wall3d_section_to_f_(C_wall3d_section&, wall3d_section_struct*);

bool operator== (const C_wall3d_section&, const C_wall3d_section&);

void operator>> (C_wall3d_section&, wall3d_section_struct*);
void operator>> (wall3d_section_struct*, C_wall3d_section&);

//--------------------------------------------------------------------
// wall3d

class wall3d_struct {};

class C_wall3d {
public:
  C_wall3d_section_array section;

};    // End Class

extern "C" void wall3d_to_c_(wall3d_struct*, C_wall3d&);
extern "C" void wall3d_to_f_(C_wall3d&, wall3d_struct*);

bool operator== (const C_wall3d&, const C_wall3d&);

void operator>> (C_wall3d&, wall3d_struct*);
void operator>> (wall3d_struct*, C_wall3d&);

//--------------------------------------------------------------------
// Control 

class control_struct {};

class C_control {
public:
  double coef;                // control coefficient
  int ix_lord;                // index to lord element
  int ix_slave;               // index to slave element
  int ix_branch;               // index to a photon line
  int ix_attrib;              // index of attribute controlled

  C_control (double c, int il, int is, int ib, int ia) :
      coef(c), ix_lord(il), ix_slave(is), ix_branch(ib), ix_attrib(ia) {}

  C_control () :
      coef(0), ix_lord(0), ix_slave(0), ix_branch(0), ix_attrib(0) {}
};    // End Class

extern "C" void control_to_c_(control_struct*, C_control&);
extern "C" void control_to_f_(C_control&, control_struct*);

bool operator== (const C_control&, const C_control&);

void operator>> (C_control&, control_struct*);
void operator>> (control_struct*, C_control&);

//--------------------------------------------------------------------
// Param 

class lat_param_struct {};

class C_lat_param {
public:
  double n_part;            // Particles/bunch (for BeamBeam elements).
  double total_length;      // total_length of ring
  double unstable_factor;   // growth rate/turn if not stable
  Real_Matrix t1_with_RF;   // Full 1-turn matrix with RF on.
  Real_Matrix t1_no_RF;     // Full 1-turn matrix with RF off.
  int particle;             // positron$, electron$, etc.
  int ix_lost;              // Index of element particle was lost at.
  int end_lost_at;          // entrance_end$ or exit_end$
  int plane_lost_at;        // x_plane$ or y_plane$
  int lattice_type;         // linear_lattice$, etc...
  int ixx;                  // Int for general use
  bool stable;              // is closed ring stable?
  bool aperture_limit_on;   // use apertures in tracking?
  bool lost;                // for use in tracking

  C_lat_param () : n_part(0), total_length(0), unstable_factor(0),
      t1_with_RF(V6_array, 6), t1_no_RF(V6_array, 6), 
      particle(0), ix_lost(0), end_lost_at(0), plane_lost_at(0), lattice_type(0), ixx(0),
      stable(1), aperture_limit_on(1), lost(0) {}

  C_lat_param (double np, double tl, double gr, Real_Matrix t1w,
    Real_Matrix t1n, int pa, int il, int ela, int pla, int lt, int ix, 
    int st, int alo, int lo) :
        n_part(np), total_length(tl), unstable_factor(gr),
        t1_with_RF(t1w), t1_no_RF(t1n), particle(pa), 
        ix_lost(il), end_lost_at(ela), plane_lost_at(pla), lattice_type(lt), ixx(ix),
        stable(st), aperture_limit_on(alo), lost(lo) {}
};    // End Class

extern "C" void lat_param_to_c_(lat_param_struct*, C_lat_param&);
extern "C" void lat_param_to_f_(C_lat_param&, lat_param_struct*);

bool operator== (const C_lat_param&, const C_lat_param&);

void operator>> (C_lat_param&, lat_param_struct*);
void operator>> (lat_param_struct*, C_lat_param&);

//--------------------------------------------------------------------
// Amode 

class anormal_mode_struct {};

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

extern "C" void amode_to_c_(anormal_mode_struct*, C_amode&);
extern "C" void amode_to_f_(C_amode&, anormal_mode_struct*);

bool operator== (const C_amode&, const C_amode&);

void operator>> (C_amode&, anormal_mode_struct*);
void operator>> (anormal_mode_struct*, C_amode&);

//--------------------------------------------------------------------
// linac_mode 

class linac_normal_mode_struct {};

class C_linac_mode {
public:
  double i2_E4;        // Integral: g^2*  gamma^4
  double i3_E7;        // Integral: g^3*  gamma^7
  double i5a_E6;       // Integral: (g^3*  H_a)*  gamma^6
  double i5b_E6;       // Integral: (g^3*  H_b)*  gamma^6
  double sig_E1;       // Energy spread after 1 pass (eV)
  double a_emittance_end;  // a mode emittance at end of linac
  double b_emittance_end;  // b mode emittance at end of linac

  C_linac_mode () :
      i2_E4(0), i3_E7(0), i5a_E6(0), i5b_E6(0), sig_E1(0),
      a_emittance_end(0), b_emittance_end(0) {}

  C_linac_mode (double i2, double i3, double i5a, double i5b, double sig,
                                                       double a, double b) :
      i2_E4(i2), i3_E7(i3), i5a_E6(i5a), i5b_E6(i5b), sig_E1(sig),
      a_emittance_end(a), b_emittance_end(b) {}

};    // End Class

extern "C" void linac_mode_to_c_(linac_normal_mode_struct*, C_linac_mode&);
extern "C" void linac_mode_to_f_(C_linac_mode&, linac_normal_mode_struct*);

bool operator== (const C_linac_mode&, const C_linac_mode&);

void operator>> (C_linac_mode&, linac_normal_mode_struct*);
void operator>> (linac_normal_mode_struct*, C_linac_mode&);

//--------------------------------------------------------------------
// Modes 

class normal_modes_struct {};

class C_modes {
public:
  double synch_int1;    // Synchrotron integrals I1
  double synch_int2;    // Synchrotron integrals I2
  double synch_int3;    // Synchrotron integrals I3
  double sigE_E;        // SigmaE/E
  double sig_z;         // Sigma_Z
  double e_loss;        // Energy loss / turn (eV)
  double rf_voltage;    // Total rfcavity voltage (eV)
  double pz_aperture;   // Momentem aperture
  C_amode  a, b, z;
  C_linac_mode lin;

  C_modes () :
      synch_int1(0), synch_int2(0), synch_int3(0), sigE_E(0), 
      sig_z(0), e_loss(0), rf_voltage(0), pz_aperture(0), a(), b(), z(), lin() {}

  C_modes (double i1, double i2, double i3, double se, double sz, 
          double el, double rf_volt, double pz,
          C_amode aa, C_amode bb, C_amode zz, C_linac_mode l) :
      synch_int1(i1), synch_int2(i2), synch_int3(i3), sigE_E(se), 
      sig_z(sz), e_loss(el), rf_voltage(rf_volt), pz_aperture(pz),
      a(aa), b(bb), z(zz), lin(l) {}

};    // End Class

extern "C" void modes_to_c_(normal_modes_struct*, C_modes&);
extern "C" void modes_to_f_(C_modes&, normal_modes_struct*);

bool operator== (const C_modes&, const C_modes&);

void operator>> (C_modes&, normal_modes_struct*);
void operator>> (normal_modes_struct*, C_modes&);

//--------------------------------------------------------------------
// Bmad_com 

class bmad_common_struct {};
class C_bmad_com;

extern "C" void bmad_com_to_c_(C_bmad_com&);
extern "C" void bmad_com_to_f_(C_bmad_com&);

class C_bmad_com {
public:
  double max_aperture_limit;   
  Real_Array d_orb;              // for the make_mat6_tracking routine
  double grad_loss_sr_wake;      // Internal var for LCavities.
  double default_ds_step;        // Number integration steps
  double significant_longitudinal_length;
  double rel_tolerance; 
  double abs_tolerance; 
  double rel_tol_adaptive_tracking;
  double abs_tol_adaptive_tracking;
  int taylor_order;              // 3rd order is default
  int default_integ_order;       // PTC integration order
  bool canonical_coords;         // Use (x, px) [not (x, x')]
  bool sr_wakes_on;              // Short range wakefields?
  bool lr_wakes_on;              // Long range wakefields
  bool mat6_track_symmetric;     // symmetric offsets
  bool auto_bookkeeper;          // Automatic bookkeeping when elements change?
  bool trans_space_charge_on;
  bool coherent_synch_rad_on;
  bool spin_tracking_on;
  bool radiation_damping_on;
  bool radiation_fluctuations_on;
  bool compute_ref_energy;
  bool conserve_taylor_maps;

  C_bmad_com () : d_orb(double(0), 6) {bmad_com_to_c_(*this);}

  C_bmad_com (double max_ap, ReArr orb, double grad_loss, double default_ds, 
                double significant, double rel_t, double abs_t, double rel_t_adapt, 
                double abs_t_adapt, int to, int io, int cc,
                int sr, int lr, int sym, int a_book,
                int tsc_on, int csr_on, int st_on, int rad_d, int rad_f, int cre,
                int con_t) :
      max_aperture_limit(max_ap), d_orb(orb, 6), grad_loss_sr_wake(grad_loss), 
      significant_longitudinal_length(significant),
      rel_tolerance(rel_t), abs_tolerance(abs_t), 
      rel_tol_adaptive_tracking(rel_t_adapt), abs_tol_adaptive_tracking(abs_t_adapt), 
      taylor_order(to), 
      default_integ_order(io), default_ds_step(default_ds), canonical_coords(cc), 
      sr_wakes_on(sr), lr_wakes_on(lr), 
      mat6_track_symmetric(sym), auto_bookkeeper(a_book), 
      trans_space_charge_on(tsc_on), coherent_synch_rad_on(csr_on),
      spin_tracking_on(st_on), radiation_damping_on(rad_d), 
      radiation_fluctuations_on(rad_f), compute_ref_energy(cre), 
      conserve_taylor_maps(con_t) {}

  C_bmad_com (double max_ap, Real_Array orb, double grad_loss, double default_ds, 
                double significant, double rel_t, double abs_t, double rel_t_adapt, 
                double abs_t_adapt, int to, int io, int cc,
                int sr, int lr, int sym, int a_book, int tsc_on, 
                int csr_on, int st_on, int rad_d, int rad_f, int cre, int con_t) :
      max_aperture_limit(max_ap), d_orb(orb), grad_loss_sr_wake(grad_loss), 
      significant_longitudinal_length(significant),
      rel_tolerance(rel_t), abs_tolerance(abs_t), 
      rel_tol_adaptive_tracking(rel_t_adapt), abs_tol_adaptive_tracking(abs_t_adapt), 
      taylor_order(to), 
      default_integ_order(io), default_ds_step(default_ds), canonical_coords(cc), 
      sr_wakes_on(sr), lr_wakes_on(lr), 
      mat6_track_symmetric(sym), auto_bookkeeper(a_book), 
      trans_space_charge_on(tsc_on), coherent_synch_rad_on(csr_on),
      spin_tracking_on(st_on), radiation_damping_on(rad_d), 
      radiation_fluctuations_on(rad_f), compute_ref_energy(cre), 
      conserve_taylor_maps(con_t) {}

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
// mode3 

class mode3_struct {};

class C_mode3 {
public:
  Real_Matrix v;
  C_twiss a, b, c;
  C_twiss x, y;

};    // End Class

extern "C" void mode3_to_c_(mode3_struct*, C_mode3&);
extern "C" void mode3_to_f_(C_mode3&, mode3_struct*);

bool operator== (const C_mode3&, const C_mode3&);

void operator>> (C_mode3&, mode3_struct*);
void operator>> (mode3_struct*, C_mode3&);

//--------------------------------------------------------------------
// space_charge 

class space_charge_struct {};

class C_space_charge {
public:
  C_coord closed_orb;      // beam orbit
  double kick_const;
  double sig_x;
  double sig_y;
  double phi;
  double sin_phi;
  double cos_phi;
  double sig_z;

};    // End Class

extern "C" void space_charge_to_c_(space_charge_struct*, C_space_charge&);
extern "C" void space_charge_to_f_(C_space_charge&, space_charge_struct*);

bool operator== (const C_space_charge&, const C_space_charge&);

void operator>> (C_space_charge&, space_charge_struct*);
void operator>> (space_charge_struct*, C_space_charge&);

//--------------------------------------------------------------------
// Ele 

class ele_struct {};

class C_ele {

public:
  string name;                  // name of element
  string type;                  // type name
  string alias;                 // Another name
  string component_name;        // Used by overlays
  string descrip;               // For general use
  C_twiss  a, b, z;             // Twiss parameters at end of element
  C_xy_disp x, y;               // Projected dispersion
  C_floor_position floor;       // Global floor position at end of ele.
  C_mode3 mode3;
  C_coord map_ref_orb_in;
  C_coord map_ref_orb_out;
  void* gen_field;              // Pointer to a PTC genfield
  C_taylor_array taylor;        // Taylor terms
  C_rf rf;                      // Fields and wakes
  C_wig_term_array wig_term;    // Wiggler Coefs
  C_space_charge space_charge;
  C_wall3d wall3d;
  Real_Array value;             // attribute values. size = N_ATTRIB_MAXX
  Real_Array gen0;              // constant part of the genfield map. size = 6
  Real_Array vec0;              // 0th order transport vector. size = 6
  Real_Matrix mat6;             // 1st order transport matrix. size = 6x6
  Real_Matrix c_mat;            // 2x2 C coupling matrix. size = 2x2
  double gamma_c;               // gamma associated with C matrix
  double s;                     // longitudinal position at the end
  double ref_time;              // Time ref particle passes exit end.
  Real_Matrix r;                // For general use. Not used by Bmad.
  Real_Array a_pole;            // multipole
  Real_Array b_pole;            // multipoles
  Real_Array const_arr;         // Working constants.
  int key;                      // key value
  int sub_key;                  // For wigglers: map_type$, periodic_type$
  int ix_ele;                   // Index in ring%ele(:) array
  int ix_branch;                // Photon line index
  int ix_value;                 // Pointer for attribute to control
  int slave_status;             // super_slave$, etc.
  int n_slave;                  // Number of slaves
  int ix1_slave;                // Start index for slave elements
  int ix2_slave;                // Stop  index for slave elements
  int lord_status;              // overlay_lord$, etc.
  int n_lord;                   // Number of lords
  int ic1_lord;                 // Start index for lord elements
  int ic2_lord;                 // Stop  index for lord elements
  int ix_pointer;               // For general use. Not used by Bmad.
  int ixx;                      // Index for Bmad internal use
  int mat6_calc_method;         // bmad_standard$, taylor$, etc.
  int tracking_method;          // bmad_standard$, taylor$, etc.
  int field_calc;               // Used with Boris, Runge-Kutta integrators.
  int ref_orbit;                // For setting the ptc kind type.
  int taylor_order;             // Order of the taylor series.
  int aperture_at;              // Where aperture is applied. exit_end$, ...
  int aperture_type;            // Where aperture is applied. exit_end$, ...
  int attribute_status;         // Element attributes have been modified?
  bool symplectify;             // Symplectify mat6 matrices.
  bool mode_flip;               // Have the normal modes traded places?
  bool multipoles_on;           // For turning multipoles on/off
  bool scale_multipoles;
  bool map_with_offsets;        // Exact radiation integral calculation?
  bool field_master;            // Calculate strength from the field value?
  bool is_on;                   // For turning element on/off.
  bool old_is_on;               // For Bmad internal use only.
  bool logic;                   // For general use. Not used by Bmad.
  bool on_a_girder;             // Have an I_Beam overlay_lord?
  bool csr_calc_on;             // Coherent synchrotron radiation calculation.
  bool offset_moves_aperture;   // element offsets affects aperture?

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
  double sigma;    
  double sigmap;

  C_mode_info () : tune(0), emit(0), chrom(0) {}
  C_mode_info (Re t, Re e, Re c, Re s, Re sp) : 
          tune(t), emit(e), chrom(c), sigma(s), sigmap(sp) {}

};    // End Class

extern "C" void mode_info_to_c_(mode_info_struct*, C_mode_info&);
extern "C" void mode_info_to_f_(C_mode_info&, mode_info_struct*);

bool operator== (const C_mode_info&, const C_mode_info&);

void operator>> (C_mode_info&, mode_info_struct*);
void operator>> (mode_info_struct*, C_mode_info&);

//--------------------------------------------------------------------
// branch

class branch_struct {};

class C_branch {
public:
  string name;
  int key;
  int ix_branch;
  int ix_from_branch;
  int ix_from_ele;
  int n_ele_track;
  int n_ele_max;
  C_ele_array ele;
  C_wall3d wall3d;
  C_lat_param param;

  C_branch () : ix_branch(0), n_ele_track(0), n_ele_max(0) {}

};

extern "C" void branch_to_c_(branch_struct*, C_branch&);
extern "C" void branch_to_f_(C_branch&, branch_struct*);

bool operator== (const C_branch&, const C_branch&);

void operator>> (C_branch&, branch_struct*);
void operator>> (branch_struct*, C_branch&);

//--------------------------------------------------------------------
// Lattice

class lat_struct {};

class C_lat {
public:
  string name;               // Name of ring given by USE statement
  string lattice;            // Lattice
  string input_file_name;    // Name of the lattice input file
  string title;              // General title
  C_mode_info a, b, z;       // Tunes, etc.
  C_lat_param param;         // Parameters
  int version;               // Version number
  int n_ele_track;           // Number of regular ring elements
  int n_ele_max;             // Index of last element used
  int n_control_max;         // Last index used in CONTROL_ array
  int n_ic_max;              // Last index used in IC_ array
  int input_taylor_order;    // As set in the input file
  C_ele ele_init;            // For use by any program
  C_ele_array ele;           // Array of elements
  C_branch_array branch;
  C_control_array control;   // control list
  C_wall3d wall3d;
  Int_Array ic;              // index to %control(:)

  C_lat () {}

  C_lat& operator= (const C_lat&);

};    // End Class

extern "C" void lat_to_c_(lat_struct*, C_lat&);
extern "C" void ring_to_f_(C_lat&, lat_struct*);

bool operator== (const C_lat&, const C_lat&);

void operator>> (C_lat&, lat_struct*);
void operator>> (lat_struct*, C_lat&);

#define CPP_AND_BMAD
#endif
