module namelist_moga

! provides:
! moga_output_file, set_chrom_x, set_chrom_y, initial_pop, seed, breeder_params, max_gen, co_limit, 
! linear_vec_cutoff, x_fp_min, x_fp_max, y_fp_min, y_fp_max, 
! nux_min, nux_max, nuy_min, nuy_max, fp_dE_neg, fp_dE_pos, n_fp_steps, mags_in

use bmad
use pisa_mod, only: breeder_params_struct
use dynap_mod, only: mag_struct

implicit none

integer, parameter :: max_mags = 200

character*100 moga_output_file
real(rp) set_chrom_x
real(rp) set_chrom_y
character*100 initial_pop
integer seed
integer generate_feasible_seeds_only
type(breeder_params_struct) breeder_params
integer max_gen
real(rp) co_limit
real(rp) linear_vec_cutoff
real(rp) x_fp_min, x_fp_max
real(rp) y_fp_min, y_fp_max
real(rp) tr_a_min, tr_a_max
real(rp) tr_b_min, tr_b_max
real(rp) fp_dE_neg, fp_dE_pos
real(rp) la_overshoot
real(rp) work_pt_x_min, work_pt_x_max
real(rp) work_pt_y_min, work_pt_y_max
integer n_fp_steps
type(mag_struct) mags_in(max_mags)

namelist / nl_moga /       moga_output_file, &
                        generate_feasible_seeds_only, &
                        set_chrom_x, &
                        set_chrom_y, &
                        initial_pop, &
                        seed, &
                        breeder_params, &
                        max_gen, &
                        co_limit, &
                        linear_vec_cutoff, &
                        la_overshoot, &
                        work_pt_x_min, &
                        work_pt_x_max, &
                        work_pt_y_min, &
                        work_pt_y_max, &
                        x_fp_min, &
                        x_fp_max, &
                        y_fp_min, &
                        y_fp_max, &
                        fp_dE_neg, &
                        fp_dE_pos, &
                        n_fp_steps, &
                        tr_a_min, &
                        tr_a_max, &
                        tr_b_min, &
                        tr_b_max, &
                        mags_in

end module
