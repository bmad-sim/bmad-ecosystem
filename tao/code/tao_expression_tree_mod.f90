module tao_expression_tree_mod

use tao_interface

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine tao_type_expression_tree (tree, indent)
!
! Routine to print an expression tree in tree form.
! Good for debugging.
!
! Input:
!   tree        -- tao_eval_node_struct: Tree to print.
!   indent      -- integer, optional: Initial indent. Default is zero.
!-

recursive subroutine tao_type_expression_tree(tree, indent)

type (tao_eval_node_struct) tree
integer, optional :: indent
integer n, ind
character(40) fmt

!

ind = integer_option(0, indent)

write(fmt, '(a, i0, a, i0, a)') '(', 4*ind+2, 'x, a, t', 4*ind+30, ', a, i0, z12)'
print fmt, trim(tree%name), ':', tree%type   !, loc(tree)

if (.not. associated(tree%node)) return

do n = 1, size(tree%node)
  call tao_type_expression_tree(tree%node(n), ind+1)
enddo

end subroutine tao_type_expression_tree

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function tao_expression_tree_to_string (tree, n_node, parent) result(str_out)
!
! Routine to convert an expression tree to a expression string.
!
! Input:
!   tree        -- tao_eval_node_struct: Tree to print.
!   n_node      -- integer, optional: Internal use only. Used with recursive calls.
!   parent      -- tao_eval_node_struct, optional: Internal use only. Used with recusive calls.
!
! Output:
!   str_out         -- character(*): Expression string.
!-

recursive function tao_expression_tree_to_string (tree, n_node, parent) result (str_out)

type (tao_eval_node_struct) tree
type (tao_eval_node_struct), optional :: parent
integer, optional :: n_node
integer n, iss, ns, n_sub
character(:), allocatable :: str_out
character(2000) str, ss(10)

!

str = ''
if (.not. associated(tree%node)) then
  n_sub = 0
else
  n_sub = size(tree%node)
endif

select case (tree%type)
case (root$, compound$)
  ! No printing

case (square_brackets$, parens$, func_parens$, curly_brackets$)
  str = tree%name(1:1)
  do n = 1, n_sub
    str = trim(str) //  tao_expression_tree_to_string(tree%node(n), n, tree)
  enddo
  str = trim(str) // tree%name(2:2)
  allocate(character(len_trim(str)) :: str_out)
  str_out = trim(str)
  if (tree%type == func_parens$) str_out = trim(parent%node(n_node+1)%name) // str_out
  return

case (function$)
  ! Handled by func_parens$ case.

case (comma$, equal$)
  if (integer_option(2, n_node) > 1) str = trim(str) // tree%name

case default
  str = trim(str) // tree%name 
end select

!

iss = 0
ns = size(ss)

do n = 1, n_sub
  if (iss == ns) then
    str = trim(str) // ss(1)
    ss(1:ns-1) = ss(2:ns)
  else
    iss = iss + 1
  endif

  ss(iss) = tao_expression_tree_to_string(tree%node(n), n, tree)

  select case (tree%node(n)%type)
  case (plus$, minus$, times$, divide$, power$)
    ss(iss-2) = trim(ss(iss-2)) // trim(ss(iss)) // ss(iss-1)
    iss = iss - 2
  case (unary_plus$, unary_minus$)
    ss(iss-1) = trim(ss(iss)) // ss(iss-1)
    iss = iss - 1
  end select
enddo

do n = 1, iss
  str = trim(str) // ss(n)
enddo

allocate(character(len_trim(str)) :: str_out)
str_out = trim(str)

end function tao_expression_tree_to_string

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine tao_re_associate_node_array(tree, n, exact)
!
! Routine to resize the tree%node(:) array.
!
! Note: The data of the array is preserved but data at the end of the
! array will be lost if n is less than the original size of the array
!
! Input:
!   tree       -- tao_eval_node_struct:
!   n          -- integer: Size wanted.
!   exact      -- logical, optional:  Default is False. If False, the size of 
!                   the output array is permitted to be larger than n. 
!
! Output:
!   tree       -- tao_eval_node_struct:
!-

subroutine tao_re_associate_node_array(tree, n, exact)

type (tao_eval_node_struct), target :: tree, temp_tree
integer n, n_old, n_save, in
logical, optional :: exact

!

if (associated(tree%node)) then
  n_old = size(tree%node)
  if (n == n_old) return
  if (.not. logic_option(.false., exact) .and. n < n_old) return
  n_save = min(n, n_old)
  temp_tree%node => tree%node
  allocate (tree%node(n))
  tree%node(1:n_save) = temp_tree%node(1:n_save)
  do in = n_save+1, n_old
    call tao_deallocate_tree(temp_tree%node(in))
  enddo
  deallocate (temp_tree%node)  
else
  allocate (tree%node(n))
endif

end subroutine tao_re_associate_node_array

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine tao_deallocate_tree (tree)
!
! Routine to deallocate tree%node(:) and everything below it
!
! Input:
!   tree      -- tao_eval_node_struct: Root of tree to deallocate.
!
! Output:
!   tree      -- tao_eval_node_struct: Deallocated tree.
!-

recursive subroutine tao_deallocate_tree(tree)

type (tao_eval_node_struct) tree
integer in

!

if (.not. associated(tree%node)) return

do in = 1, size(tree%node)
  call tao_deallocate_tree(tree%node(in))
enddo

deallocate(tree%node)

end subroutine 

end module
