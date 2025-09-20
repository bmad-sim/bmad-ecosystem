module expression_mod

use bmad_struct

implicit none

type (expression_tree_struct), pointer :: tree_root

private pushit

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine expression_string_to_tree (string, root_tree, err_flag, err_str)
!
! Routine to create an expression tree array which can be used 
! to evaluate an arithmethic expression.
!
! Also see:
!   expression_value
!   expression_tree_value
!   deallocate_expression_tree
!
! Important! trees use pointers as opposed to allocatable arrays due to the ifort compiler not being able to 
! handle node%node(:) being an allocatable array. Thus deallocate_expression_tree must be called before
! any tree instance goes out of scope.
!
! Note types used:
!   plus$, minus$, times$, divide$, power$, unary_minus$, unary_plus$
!   constant$, numeric$, variable$, function$
!   root$, parens$, func_parens$, square_brackets$, curly_brackets$
!   arrow$, equal$, colon$, double_colon$, vertical_bar$, compound$
! 
! An expression string will be split on:
!   Two character operators: "->", "::" 
!   operators: +-*/^=:
!   brackets: [](){}
!   comma: ,
!
! Root node name is "root" and is of type root$
! Brackets in the expression string must be matched.
! The corresponding tree node will have a name / type of:
!   "[]" / square_brackets$,    "()" / parens$ or func_parens$,   "{}" / curley_brackets$
!
! The root node, equal nodes, and all bracket nodes, will have an array of child nodes all of which will be comma nodes.
! Example:
!   "[A, B]" 
! will translate to a "[]" node with two comma children and the first comma child will have a 
! single child "A" and the second comma child will have a single child "B".
! Example:
!   "(A)"
! will translate to a "()" node with one comma child and this comma child will have a single child "A".
!
! Exception: If the string is an equation. For example, "A, B = C, D, Z". In this case the root node
! will have two equal node children (and not comma children), The first equal node represents the left 
! hand side of the equation and this node will have two comma children. The second equal node child will
! have three comma children.
!
! Compound variables are something like "data::orxit.x" (this is a Tao construct) which get 
! translated to a compound$ node with children:
!   "data",  "::", "orbit.x"
! Also functions line "atan()" are considered compound vars with children
!   "atan",  "()"
!
! The funciton argument of a species related function like "He++" in the construct "mass_of(He++)",
!  will not get split and will get marked as a species_const$. 
!
! Input:
!   root_tree -- expression_tree_struct: Only used when recursively called.
!   string    -- character(*): Expression to be converted.
!
! Output:
!   tree      -- expression_tree_struct: Expression evaluation tree.
!   err_flag  -- logical: Set True if there is an error (EG divide by 0).
!   err_str   -- character(*): String describing the error. Make length large to hold the expression.
!-

subroutine expression_string_to_tree (string, root_tree, err_flag, err_str)

type (expression_tree_struct), target :: root_tree

logical err_flag

character(*) string, err_str
character(len(string)) parse_line

!

tree_root => root_tree  ! For error messages.

err_flag = .false.
err_str = ''
parse_line = string

! parse_pass: Parse expression string into node chunks and where there are brackets, 
!   form a subtree with a bracket node as the root.
! comma_pass: Where there are comma deliminated chunks, reform into a set of subtrees, one for each chunk
!   and with a comma node as the root for each chunk.
! node_markup_pass: Assign the appropriate node%type for all nodes.
! reverse_polish_pass: For any tree%node(:) array that represents an expression, Reformulate in reverse Polish.
!   Also form compound_var$ subtrees as needed.

call deallocate_expression_tree(root_tree)
call parse_pass(parse_line, root_tree, err_flag, err_str, 'root'); if (err_flag) return
call node_markup_pass(root_tree, err_flag, err_str); if (err_flag) return
call comma_pass(root_tree, err_flag, err_str); if (err_flag) return
call reverse_polish_pass(root_tree, err_flag, err_str); if (err_flag) return

!------------------------------------------------------------------------
contains

recursive subroutine parse_pass(parse_line, tree, err_flag, err_str, tree_name)

type (expression_tree_struct), target :: tree, t2
type (expression_tree_struct), pointer :: node 

integer i, id, ixe, ix_word, n_node

logical delim_found, do_combine
logical err_flag

character(*) parse_line, err_str, tree_name
character(2) delim, cc
character(80) word, word2

! 

select case (tree_name)
case ('=');     tree%type = equal$;            tree%name = '='
case ('[');     tree%type = square_brackets$;  tree%name = '[]'
case ('(');     tree%type = parens$;           tree%name = '()'
case ('{');     tree%type = curly_brackets$;   tree%name = '{}'
case ('root');  tree%type = root$;             tree%name = 'root'
end select

tree%value = 0
n_node = 0

!

parsing_loop: do
  call get_next_chunk (parse_line, word, ix_word, ' []+-*/()^,{}=:|', delim, delim_found)
  if (ix_word == 0 .and. .not. delim_found) then
    if (tree%type /= root$ .and. tree%type /= equal$) then
      call set_this_err('Mismatched brackets. Cannot find closing bracket for opening: ' // quote(tree%name(1:1)) // &
                        '  In: ' // quote(string), err_str, err_flag)
    endif
    call re_associate_node_array(tree, n_node, exact = .true.)
    return
  endif

  select case (delim)
  case ('[', '(', '{')
    call push_node(tree, n_node, word)
    call increment_n_node(tree, n_node)
    call parse_pass(parse_line, tree%node(n_node), err_flag, err_str, delim); if (err_flag) return

  case ('+', '-')
    do_combine = (ix_word > 1)

    if (do_combine) then
      cc = upcase(word(ix_word:ix_word))
      if (cc /= 'E' .and. cc /= 'D') do_combine = .false.
      if (.not. is_integer(parse_line, delims = '+-*/()^,[]{}=:| ', ix_word = ixe)) do_combine = .false.

      do i = 1, ix_word-1
        if (index('.0123456789', word(i:i)) == 0) do_combine = .false.
      enddo
    endif

    if (do_combine) then
      word = trim(word) // trim(delim) // parse_line(1:ixe)
      parse_line = parse_line(ixe+1:)
      call push_node(tree, n_node, word)

    else
      call push_node(tree, n_node, word)
      call push_node(tree, n_node, delim)
    endif

  case (']', ')', '}')
    call push_node(tree, n_node, word)

    if (tree%name(2:2) /= delim) then
      if (tree%name == 'root') then
        call set_this_err('End bracket: ' // quote(delim) // ' has no beginning bracket.' // &
                          '  In: ' // quote(string), err_str, err_flag)
      else
        call set_this_err('Brackets mismatch. Opening: ' // quote(tree%name) // ' has ending: ' // quote(delim) // &
                          '  In: ' // quote(string), err_str, err_flag)
      endif
    endif
    exit

  case ('=')
    call push_node(tree, n_node, word)
    call re_associate_node_array(tree, n_node, exact = .true.)
    t2%node => tree%node
    allocate(tree%node(2))
    n_node = 2

    node => tree%node(1)
    node%type = equal$
    node%name = '='
    node%node => t2%node

    call parse_pass(parse_line, tree%node(2), err_flag, err_str, delim); if (err_flag) return

  case default
    call push_node(tree, n_node, word)
    call push_node(tree, n_node, delim, keep_blank = delim_found)
  end select
enddo parsing_loop

call re_associate_node_array(tree, n_node, exact = .true.)

end subroutine parse_pass

!------------------------------------------------------------------------
! contains

recursive subroutine node_markup_pass(tree, err_flag, err_str)

type (expression_tree_struct), target :: tree
type (expression_tree_struct), pointer :: node

integer in, n_node
logical err_flag
character(*) err_str

!

if (.not. associated(tree%node)) return
n_node = size(tree%node)

do in = 1, n_node
  node => tree%node(in)
  call node_markup_pass(node, err_flag, err_str); if (err_flag) return

  ! Set type

  select case (upcase(node%name))
  case ('MIN', 'MAX', 'COT', 'CSC', 'SEC', 'SIN', 'ACOSH', 'ATANH', &
           'SINC', 'COS', 'TAN', 'ASIN', 'ACOS', 'ATAN', 'ATAN2', 'MODULO', &
           'ABS', 'SQRT', 'LOG', 'EXP', 'FACTORIAL', 'RAN', 'RAN_GAUSS', 'INT', &
           'SIGN', 'NINT', 'FLOOR', 'CEILING', 'CHARGE_OF', 'MASS_OF', 'SPECIES', 'ANTIPARTICLE', &
           'ANOMALOUS_MOMENT_OF', 'COTH', 'SINH', 'COSH', 'TANH', 'ACOTH', 'ASINH')
    node%type = function$
  case ('->');                   node%type = arrow$
  case ('*');                    node%type = times$
  case ('/');                    node%type = divide$
  case ('^');                    node%type = power$
  case (':');                    node%type = colon$
  case ('::');                   node%type = double_colon$
  case ('|');                    node%type = vertical_bar$
  case (',');                    node%type = comma$
  case (' ');                    node%type = blank$
  case ('{}', '[]', '=')        ! Already marked

  case ('()')
    if (is_alphabetic(tree%node(max(1,in-1))%name(1:1))) node%type = func_parens$

  case ('+')
    node%type = plus$
    if (in == 1) then
      node%type = unary_plus$
    elseif (index('+-*/^', trim(tree%node(in-1)%name)) > 0) then
      node%type = unary_plus$
    endif

  case ('-')
    node%type = minus$
    if (in == 1) then
      node%type = unary_minus$
    elseif (index('+-*/^', trim(tree%node(in-1)%name)) > 0) then
      node%type = unary_minus$
    endif

  case default
    if (is_real(node%name)) then
      node%type = constant$
    else
      node%type = variable$
    endif
  end select

  ! Some error checks

  select case (node%type)
  case (function$)
    if (in == n_node) then
      call set_this_err('Function: ' // quote(node%name) // &
                           ' is not followed by parenteses character "(" in: ' // string, err_str, err_flag)
    elseif (tree%node(in+1)%type /= parens$) then
      call set_this_err('Function: ' // quote(node%name) // &
                           ' is not followed by parenteses character "(" in: ' // string, err_str, err_flag)
    endif
  end select
enddo

! Combine unary minus or plus node with following constant node.

do in = n_node-1, 1, -1
  if (tree%node(in+1)%type /= constant$) cycle
  select case (tree%node(in)%type)
  case (unary_minus$, unary_plus$)
    tree%node(in+1)%name = tree%node(in)%name(1:1) // tree%node(in+1)%name
    n_node = n_node - 1
    tree%node(in:n_node) = tree%node(in+1:n_node+1)
    nullify(tree%node(n_node+1)%node)
  end select
enddo

call re_associate_node_array(tree, n_node, exact = .true.)


end subroutine node_markup_pass

!------------------------------------------------------------------------
! contains

recursive subroutine comma_pass(tree, err_flag, err_str)

type (expression_tree_struct), target :: tree, t2
type (expression_tree_struct), pointer :: t2a

integer n_comma, na, n0, ia, nn
logical err_flag
character(*) err_str

! Exception: root node with equal subnodes does not get a comma layer.

if (.not. associated(tree%node)) return
nn = size(tree%node)
if (nn == 0) return

if (tree%type == root$ .and. tree%node(1)%type == equal$) then
  do na = 1, nn
    call comma_pass(tree%node(na), err_flag, err_str); if (err_flag) return
  enddo
  return
endif

! Add comma layer

n_comma = 0
n0 = 1

do na = 1, nn
  select case (tree%node(na)%type)
  case (square_brackets$, parens$, func_parens$, curly_brackets$, equal$)
    call comma_pass(tree%node(na), err_flag, err_str); if (err_flag) return

  case (comma$)
    call increment_n_node(t2, n_comma)
    t2a => t2%node(n_comma)
    t2a%name = ','
    t2a%type = comma$

    allocate(t2a%node(na-n0))
    do ia = 1, na - n0
      t2a%node(ia) = tree%node(ia + n0 - 1)
    enddo
    n0 = na + 1
  end select
end do

call increment_n_node(t2, n_comma)
t2a => t2%node(n_comma)
t2a%name = ','
t2a%type = comma$
allocate(t2a%node(nn + 1 - n0))
do ia = 1, nn + 1 - n0
  t2a%node(ia) = tree%node(ia + n0 - 1)
enddo

call re_associate_node_array(t2, n_comma, exact = .true.)
deallocate(tree%node)
tree%node => t2%node

end subroutine comma_pass

!------------------------------------------------------------------------
! contains

recursive subroutine reverse_polish_pass(tree, err_flag, err_str)

type (expression_tree_struct), target :: tree, t2, op(20)
type (expression_tree_struct), pointer :: node2, snode

integer i, it2, it, n_node, i_op, n_nonop, level
logical err_flag, has_op, callit
character(*) err_str

! If tree%node(:) array does not represent an expression, skip reverse Polish step.

if (.not. associated(tree%node)) return
n_node = size(tree%node)

has_op = .false.
do it2 = 1, n_node
  node2 => tree%node(it2)

  ! Species names "He++" are not to be put are to be consoladated
  callit = .true.
  if (node2%type == func_parens$) then
    select case (tree%node(it2-1)%name)
    case ('mass_of', 'charge_of', 'anomalous_moment_of', 'species')
      callit = .false.
      snode => node2%node(1)   ! Comma node
      t2%node => snode%node
      allocate (snode%node(1))
      snode%node(1)%type = species$
      do i = 1, size(t2%node)
        snode%node(1)%name = trim(snode%node(1)%name) // t2%node(i)%name
      enddo
      call deallocate_expression_tree(t2)
    end select
  endif

  if (callit) call reverse_polish_pass(node2, err_flag, err_str); if (err_flag) return

  select case (node2%type)
  case (plus$, minus$, times$, divide$, power$, unary_plus$, unary_minus$, func_parens$)
    has_op = .true.
  end select
enddo

if (.not. has_op) return

!

t2%node => tree%node
allocate(tree%node(n_node))

i_op = 0
it = 0
n_nonop = 0

do it2 = 1, n_node
  node2 => t2%node(it2)

  select case (node2%type)
  case (plus$, minus$, times$, divide$, power$)
    if (n_nonop > 1) call make_compound_node(tree, it, n_nonop)
    n_nonop = 0

    ! See if there are operations on the OP stack that need to be transferred
    ! to the tree%node array.

    do i = i_op, 1, -1
      if (expression_eval_level(op(i)%type) < expression_eval_level(node2%type)) exit
      it = it + 1
      tree%node(it) = op(i)
    enddo
    i_op = i

    i_op = i_op + 1
    op(i_op) = node2

  case (unary_plus$, unary_minus$)
    i_op = i_op + 1
    op(i_op) = node2

  case (func_parens$)
    it = it + 1
    tree%node(it) = tree%node(it-1)
    tree%node(it-1) = node2
    nullify(node2%node)
    n_nonop = n_nonop + 1

  case default
    it = it + 1
    tree%node(it) = node2
    nullify(node2%node)
    n_nonop = n_nonop + 1
  end select
enddo

if (n_nonop > 1) call make_compound_node(tree, it, n_nonop)

do i = i_op, 1, -1
  it = it + 1
  tree%node(it) = op(i)
enddo

! Error check

level = 0

do i = 1, it
  select case (tree%node(i)%type)
  case (plus$, minus$, times$, divide$, power$)
    level = level - 1
    if (level < 1) then
      call set_this_err('Bad expression1:' // quote(string), err_str, err_flag)
      return
    endif
  
  case (unary_plus$, unary_minus$)
    ! No op

  case (func_parens$)
    ! No op

  case default
    level = level + 1
  end select
enddo

if (level /= 1) then
  call set_this_err('Bad expression2:' // quote(string), err_str, err_flag)
  return
endif

!

deallocate(t2%node)
call re_associate_node_array(tree, it, exact = .true.)

end subroutine reverse_polish_pass

!------------------------------------------------------------------------
! contains

subroutine make_compound_node(tree, it, n_nonop)

type (expression_tree_struct), target :: tree, t2, compound
integer it, n_nonop
integer j

!

allocate(compound%node(n_nonop))
compound%type = compound$
compound%name = 'compound'
do j = 1, n_nonop
  compound%node(j) = tree%node(it-n_nonop+j)
  nullify(tree%node(it-n_nonop+j)%node)
enddo

it = it - n_nonop + 1
tree%node(it) = compound

end subroutine make_compound_node

!------------------------------------------------------------------------
! contains

subroutine push_node(tree, n_node, str, keep_blank)

type (expression_tree_struct), target :: tree
integer n_node
logical, optional :: keep_blank
character(*) str

if (len_trim(str) == 0 .and. .not. logic_option(.false., keep_blank)) return
call increment_n_node(tree, n_node)
tree%node(n_node)%name = str

end subroutine push_node

!-------------------------------------------------------------------------
! contains

subroutine get_next_chunk (parse_line, word, ix_word, delim_list, delim, delim_found)

character(*) parse_line, word, delim_list, delim

integer ix_word

logical delim_found

!

call word_read (parse_line, delim_list, word, ix_word, delim, delim_found, parse_line)

if (delim == '-' .and. parse_line(1:1) == '>' .or. &
    delim == ':' .and. parse_line(1:1) == ':') then
  delim = delim(1:1) // parse_line(1:1)
  call string_trim(parse_line(2:), parse_line, ix_word)
endif



end subroutine get_next_chunk

!------------------------------------------------------------------------
! contains

subroutine increment_n_node(tree, n_node)

type (expression_tree_struct), target :: tree
integer n_node

n_node = n_node + 1
if (associated(tree%node)) then
  if (n_node > size(tree%node)) call re_associate_node_array(tree, n_node+10)
else
  call re_associate_node_array(tree, n_node+10)
endif

end subroutine increment_n_node

!------------------------------------------------------------------------
! contains

subroutine set_this_err(str, err_str, err_flag)

character(*) str, err_str
logical err_flag
integer n

!

err_flag = .true.
n = min(len(str), len(err_str))
err_str = str(1:n)

end subroutine set_this_err

end subroutine expression_string_to_tree

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine deallocate_expression_tree(tree)
!
! Routine to deallocate an expression tree.
!
! Input:
!   tree      -- expression_tree_struct: Tree to deallocate.
!
! Output:
!   tree      -- expression_tree_struct: Deallocated tree.
!-

recursive subroutine deallocate_expression_tree(tree)

type (expression_tree_struct) tree
integer in

!

if (associated(tree%node)) then
  do in = 1, size(tree%node)
    call deallocate_expression_tree(tree%node(in))
  enddo

  deallocate(tree%node)
endif

end subroutine deallocate_expression_tree

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine type_expression_tree (tree, indent)
!
! Routine to print an expression tree in tree form.
! Good for debugging.
!
! Input:
!   tree        -- expression_tree_struct: Tree to print.
!   indent      -- integer, optional: Initial indent. Default is zero.
!-

recursive subroutine type_expression_tree(tree, indent)

type (expression_tree_struct) tree
integer, optional :: indent
integer n, ind
character(40) fmt

!

ind = integer_option(0, indent)

write(fmt, '(a, i0, a, i0, a)') '(', 4*ind+2, 'x, a, t', 4*ind+30, ', a, i0, z12)'
print fmt, trim(tree%name), ':', tree%type   !, loc(tree)

if (.not. associated(tree%node)) return

do n = 1, size(tree%node)
  call type_expression_tree(tree%node(n), ind+1)
enddo

end subroutine type_expression_tree

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function expression_tree_to_string (tree, include_root, n_node, parent) result(str_out)
!
! Routine to convert an expression tree to a expression string.
!
! Input:
!   tree          -- expression_tree_struct: Root of tree to print.
!   include_root  -- logical, optional: Default is True. If True, do not inculde in the output 
!                     string the root node. Note: If the root node is of type root$, this node is
!                     always ignored. 
!   n_node        -- integer, optional: Node index. parent%node(n_node) === tree.
!                       Internal use only. Used with recursive calls.
!   parent        -- expression_tree_struct, optional: Internal use only. Used with recusive calls.
!
! Output:
!   str_out       -- character(*): Expression string.
!-

recursive function expression_tree_to_string (tree, include_root, n_node, parent) result (str_out)

type (expression_tree_struct) tree
type (expression_tree_struct), optional :: parent
integer, optional :: n_node
integer n, iss, ns, n_sub
character(:), allocatable :: str_out
character(2000) str, ss(10)
logical, optional :: include_root
logical rt_inc

!

str = ''
rt_inc = logic_option(.true., include_root)

if (.not. associated(tree%node)) then
  n_sub = 0
else
  n_sub = size(tree%node)
endif

!

select case (tree%type)
case (root$, compound$)
  ! No printing

case (square_brackets$, parens$, func_parens$, curly_brackets$)
  if (rt_inc) str = tree%name(1:1)
  do n = 1, n_sub
    str = trim(str) //  expression_tree_to_string(tree%node(n), .true., n, tree)
  enddo
  if (rt_inc) str = trim(str) // tree%name(2:2)
  allocate(character(len_trim(str)) :: str_out)
  str_out = trim(str)
  if (tree%type == func_parens$ .and. rt_inc) str_out = trim(parent%node(n_node+1)%name) // str_out
  return

case (function$)
  ! Handled by func_parens$ case.

case (comma$, equal$)
  if (integer_option(2, n_node) > 1 .and. rt_inc) str = tree%name

case default
  if (rt_inc) str = tree%name 
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

  ss(iss) = expression_tree_to_string(tree%node(n), .true., n, tree)

  select case (tree%node(n)%type)
  case (plus$, minus$, times$, divide$, power$)
    if (iss > 2) then
      ss(iss-2) = trim(ss(iss-2)) // trim(ss(iss)) // ss(iss-1)
      iss = iss - 2
    endif
  case (unary_plus$, unary_minus$)
    ss(iss-1) = trim(ss(iss)) // ss(iss-1)
    iss = iss - 1
  end select
enddo

do n = 1, iss
  if (ss(max(1,n-1)) == ' ') then
    str = trim(str) // ' ' //  ss(n)
  else
    str = trim(str) // ss(n)
  endif
enddo

allocate(character(len_trim(str)) :: str_out)
str_out = trim(str)

end function expression_tree_to_string

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine re_associate_node_array(tree, n, exact)
!
! Routine to resize the tree%node(:) array.
!
! Note: The data of the array is preserved but data at the end of the
! array will be lost if n is less than the original size of the array
!
! Input:
!   tree       -- expression_tree_struct:
!   n          -- integer: Size wanted.
!   exact      -- logical, optional:  Default is False. If False, the size of 
!                   the output array is permitted to be larger than n. 
!
! Output:
!   tree       -- expression_tree_struct:
!-

subroutine re_associate_node_array(tree, n, exact)

type (expression_tree_struct), target :: tree, temp_tree
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
    call deallocate_tree(temp_tree%node(in))
  enddo
  deallocate (temp_tree%node)  
else
  allocate (tree%node(n))
endif

end subroutine re_associate_node_array

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine deallocate_tree (tree)
!
! Routine to deallocate tree%node(:) and everything below it
!
! Input:
!   tree      -- expression_tree_struct: Root of tree to deallocate.
!
! Output:
!   tree      -- expression_tree_struct: Deallocated tree.
!-

recursive subroutine deallocate_tree(tree)

type (expression_tree_struct) tree
integer in

!

if (.not. associated(tree%node)) return

do in = 1, size(tree%node)
  call deallocate_tree(tree%node(in))
enddo

deallocate(tree%node)

end subroutine 

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine expression_string_to_stack (string, stack, n_stack, err_flag, err_str)
!
! This routine creates an expression stack array which can be used 
! to evaluate an arithmethic expression.
!
! Stack end elements not used are marked stack(i)%type = end_stack$
!
! Stack elements with stack(i)%type = variable$ are elements that need
! to be evaluated before calling expression_stack_value.
!
! Also see:
!   expression_value
!   expression_stack_value
!
! Input:
!   string    -- character(*): Expression to be converted.
!
! Output:
!   stack(:)  -- expression_atom_struct, allocatable: Expression evaluation stack.
!   n_stack   -- integer: number of "atoms" used by the expression
!   err_flag  -- logical: Set True if there is an error (EG divide by 0).
!   err_str   -- character(*): String describing the error.
!-

subroutine expression_string_to_stack (string, stack, n_stack, err_flag, err_str)

type expression_func_struct
  character(12) :: name = ''      ! Name of function
  integer :: n_arg_target = 0     ! Number of arguments the function should have. -1 => 0 or 1 arg
  integer :: n_arg_count = 0      ! Number of arguments found.
end type

type (expression_atom_struct), allocatable :: stack(:)
type (expression_func_struct) func(0:20)

integer i_op, i, var_type, ixe
integer op(100), ix_word, ix, i_delim, i2, ix_word2, n_stack, op0, n_func

real(rp) value

character(*) string, err_str
character(1) delim, old_delim, cc
character(80) word, word2
character(len(string)) parse_line

logical delim_found, do_combine
logical err_flag, err

! The general idea is to rewrite the expression on a stack in reverse polish.
! Reverse polish means that the operand goes last so that 2 * 3 is written 
! on the stack as: [2, 3, *]

! Since operations move towards the end of the stack we need a separate
! stack called op which keeps track of what operations have not yet
! been put on stack.

! init

delim = ''
err_flag = .true.
n_stack = 0
n_func = 0
i_op = 0
if (.not. allocated(stack)) allocate(stack(10))
stack(:)%type = end_stack$
parse_line = string

! parsing loop to build up the stack.

parsing_loop: do

  ! Get a word. If last thing was mass_of, etc., then next word is a species name.
  ! In this case, word is everything up to next parens.

  old_delim = delim
  op0 = 0
  var_type = variable$
  if (i_op > 1) op0 = op(i_op-1) ! If parsed "mass_of(" then op(i_op) corresponds to "("
  select case (op0)
  case (mass_of$, charge_of$, anomalous_moment_of$, species$, antiparticle$) 
    call get_next_chunk (parse_line, word, ix_word, '()', delim, delim_found)
    var_type = species_const$
  case default
    call get_next_chunk (parse_line, word, ix_word, '[]+-*/()^,:} ', delim, delim_found)
  end select

  if (delim == '*' .and. word(1:1) == '*') then
    err_str = 'EXPONENTIATION SYMBOL IS "^" AS OPPOSED TO "**"!'
    return
  endif

  ! Args are counted counted at the beginning of the function and at each comma.

  if (n_func > 0 .and. ix_word /= 0) then
    if (func(n_func)%n_arg_count == 0) func(n_func)%n_arg_count = func(n_func)%n_arg_count + 1
  endif

  !--------------------------
  ! Preliminary: If we have split up something that should have not been split
  ! then put it back together again...

  ! just make sure we are not chopping a number in two, e.g. "3.5d-7" should not
  ! get split at the "-" even though "-" is a delimiter

  do_combine = ((delim == '-' .or. delim == '+') .and. ix_word > 1)

  if (do_combine) then
    cc = upcase(word(ix_word:ix_word))
    if (cc /= 'E' .and. cc /= 'D') do_combine = .false.

    if (.not. is_integer(parse_line, delims = '+-*/()^,:[]}', ix_word = ixe)) do_combine = .false.

    do i = 1, ix_word-1
      if (index('.0123456789', word(i:i)) == 0) do_combine = .false.
    enddo
  endif

  if (do_combine) then
    word = trim(word) // trim(delim) // parse_line(1:ixe)
    parse_line = parse_line(ixe+1:)
    ix_word = len_trim(word)
  
    call get_next_chunk (parse_line, word2, ix_word2, '+-*/()^,:[}', delim, delim_found)
    if (ix_word2 /= 0) then
      err_str = 'Malformed number: ' // trim(word) // word2
      return
    endif
  endif

  ! Something like "lcav[lr(2).freq]" will get split on the "["

  if (delim == '[') then
    call get_next_chunk (parse_line, word2, ix_word2, ']', delim, delim_found)
    if (delim /= ']') then
      err_str = 'No "]" found to match "["'
      return
    endif
    word = trim(word) // '[' // trim(word2) // ']'
    call get_next_chunk (parse_line, word2, ix_word2, '+-*/()^,:} ', delim, delim_found)
    word = trim(word) // word2
    ix_word = len_trim(word)
  endif

  !---------------------------
  ! Now see what we got...

  ! For a "(" delim we must have a function

  if (delim == '(') then
    if (ix_word /= 0) then
      n_func = n_func + 1
      func(n_func) = expression_func_struct(upcase(word), 1, 0)
      select case (upcase(word))
      case ('COT');           call pushit (op, i_op, cot$)
      case ('CSC');           call pushit (op, i_op, csc$)
      case ('SEC');           call pushit (op, i_op, sec$)
      case ('SIN');           call pushit (op, i_op, sin$)
      case ('SINC');          call pushit (op, i_op, sinc$)
      case ('COS');           call pushit (op, i_op, cos$)
      case ('TAN');           call pushit (op, i_op, tan$)
      case ('ASIN');          call pushit (op, i_op, asin$)
      case ('ACOS');          call pushit (op, i_op, acos$)
      case ('ATAN');          call pushit (op, i_op, atan$)
      case ('ATAN2')
        call pushit (op, i_op, atan2$)
        func(n_func)%n_arg_target = 2
      case ('MODULO')
        call pushit (op, i_op, modulo$)
        func(n_func)%n_arg_target = 2
      case ('ABS');           call pushit (op, i_op, abs$)
      case ('SQRT');          call pushit (op, i_op, sqrt$)
      case ('LOG');           call pushit (op, i_op, log$)
      case ('EXP');           call pushit (op, i_op, exp$)
      case ('FACTORIAL');     call pushit (op, i_op, factorial$)
      case ('RAN')
        call pushit (op, i_op, ran$)
        func(n_func)%n_arg_target = 0
      case ('RAN_GAUSS') 
        call pushit (op, i_op, ran_gauss$)
        func(n_func)%n_arg_target = -1        ! 0 or 1 args.
      case ('INT');           call pushit (op, i_op, int$)
      case ('SIGN');          call pushit (op, i_op, sign$)
      case ('NINT');          call pushit (op, i_op, nint$)
      case ('FLOOR');         call pushit (op, i_op, floor$)
      case ('CEILING');       call pushit (op, i_op, ceiling$)
      case ('CHARGE_OF');     call pushit (op, i_op, charge_of$)
      case ('MASS_OF');       call pushit (op, i_op, mass_of$)
      case ('SPECIES');       call pushit (op, i_op, species$)
      case ('ANTIPARTICLE');  call pushit (op, i_op, antiparticle$)
      case ('ANOMALOUS_MOMENT_OF'); call pushit (op, i_op, anomalous_moment_of$)
      case ('COTH');           call pushit (op, i_op, coth$)
      case ('SINH');           call pushit (op, i_op, sinh$)
      case ('COSH');           call pushit (op, i_op, cosh$)
      case ('TANH');           call pushit (op, i_op, tanh$)
      case ('ACOTH');          call pushit (op, i_op, acoth$)
      case ('ASINH');          call pushit (op, i_op, asinh$)
      case ('ACOSH');          call pushit (op, i_op, acosh$)
      case ('ATANH');          call pushit (op, i_op, atanh$)
      case default
        err_str = 'UNEXPECTED CHARACTERS IN EXPRESSION BEFORE "(": ' // word
        return
      end select

      call pushit (op, i_op, l_func_parens$)

    else
      call pushit (op, i_op, l_parens$)
    endif

    cycle parsing_loop

  ! for a unary "-"

  elseif (delim == '-' .and. ix_word == 0) then
    call pushit (op, i_op, unary_minus$)
    cycle parsing_loop

  ! for a unary "+"

  elseif (delim == '+' .and. ix_word == 0) then
    call pushit (op, i_op, unary_plus$)
    cycle parsing_loop

  ! for a ")" delim

  elseif (delim == ')') then
    if (ix_word == 0) then
      if (n_func == 0 .or. (func(n_func)%n_arg_target /= 0 .and. func(n_func)%n_arg_target /= -1)) then
        err_str = 'CONSTANT OR VARIABLE MISSING BEFORE ")"'
        return
      endif
    else
      call push_numeric_or_var(word, err, op, i_op, stack, n_stack, var_type); if (err) return
    endif

    do
      do i = i_op, 1, -1     ! release pending ops
        if (op(i) == l_parens$ .or. op(i) == l_func_parens$) exit          ! break do loop
        call pushit_stack (stack, n_stack, op(i))
      enddo

      if (i == 0) then
        err_str = 'UNMATCHED ")" IN EXPRESSION'
        return
      endif

      i_op = i - 1

      if (op(i) == l_func_parens$) then
        if (func(n_func)%n_arg_target == -1) then
          if (func(n_func)%n_arg_count /= 0 .and. func(n_func)%n_arg_count /= 1) then
            err_str = 'FUNCTION: ' // trim(func(n_func)%name) // ' DOES NOT HAVE 0 OR 1 ARGUMENTS.'
            return
          endif
          call pushit_stack (stack, n_stack, arg_count$)
          stack(n_stack)%value = func(n_func)%n_arg_count

        else
          if (func(n_func)%n_arg_count /= func(n_func)%n_arg_target) then
            err_str = 'FUNCTION: ' // trim(func(n_func)%name) // ' DOES NOT HAVE THE CORRECT NUMBER OF ARGUMENTS.'
            return
          endif
        endif

        n_func = n_func - 1
      endif

      call get_next_chunk (parse_line, word, ix_word, '+-*/()^,:}', delim, delim_found)
      if (ix_word /= 0) then
        err_str = 'UNEXPECTED CHARACTERS IN EXPRESSION AFTER ")"'
        return
      endif

      if (delim /= ')') exit  ! if no more ')' then no need to release more
    enddo


    if (delim == '(') then
      err_str = '")(" CONSTRUCT DOES NOT MAKE SENSE'
      return
    endif

  ! For binary "+-/*^" delims

  else
    if (ix_word == 0) then
      if (old_delim == "*" .and. delim == "*") then
        err_str = 'EXPONENTIATION "**" NEEDS TO BE REPLACED BY "^"'
      else
        err_str = 'CONSTANT OR VARIABLE MISSING'
      endif
      return
    endif
    call push_numeric_or_var (word, err, op, i_op, stack, n_stack, var_type); if (err) return
  endif

  ! If we are here then we have an operation that is waiting to be identified

  if (.not. delim_found) delim = ':'

  select case (delim)
  case ('+')
    i_delim = plus$
  case ('-')
    i_delim = minus$
  case ('*')
    i_delim = times$
  case ('/')
    i_delim = divide$
  case ('^')
    i_delim = power$
  case (',')
    i_delim = comma$
    func(n_func)%n_arg_count = func(n_func)%n_arg_count + 1
  case ('}', ':', ')')   ! End of expression delims
    i_delim = no_delim$
  case default
    err_str = 'MALFORMED EXPRESSION'
    return
  end select

  ! Now see if there are operations on the OP stack that need to be transferred
  ! to the STACK stack

  do i = i_op, 1, -1
    if (expression_eval_level(op(i)) < expression_eval_level(i_delim)) exit

    if (op(i) == l_parens$) then
      err_str = 'UNMATCHED "("'
      return
    endif

    if (op(i) == l_func_parens$) then
      if (i_delim /= comma$) then
        err_str = 'UNMATCHED "("'
        return
      endif
      i_op = i
      cycle parsing_loop
    endif

    call pushit_stack (stack, n_stack, op(i))
  enddo

  i_op = i

  ! put the pending operation on the OP stack

  if (i_delim == no_delim$ .or. i_delim == comma$) then
    exit parsing_loop
  else
    call pushit (op, i_op, i_delim)
  endif

enddo parsing_loop

!------------------------------------------------------------------
! Go through the stack and perform the operations

if (i_op /= 0) then
  err_str = 'UNMATCHED "("'
  return
endif

if (n_stack == 0) then
  err_str = 'NO VALUE FOUND'
  return
endif

err_flag = .false.

!-------------------------------------------------------------------------
contains

subroutine get_next_chunk (parse_line, word, ix_word, delim_list, delim, delim_found)

character(*) parse_line, word, delim_list, delim

integer ix_word

logical delim_found

!

old_delim = delim
call word_read (parse_line, delim_list, word, ix_word, delim, delim_found, parse_line)

end subroutine get_next_chunk

!-------------------------------------------------------------------------
! contains

subroutine push_numeric_or_var (word, err, op, i_op, stack, n_stack, var_type)

type (expression_atom_struct), allocatable :: stack(:)
logical err
integer op(:), i_op, n_stack, var_type, ios, ix
character(*) word

!

err = .true.

if (is_real(word)) then
  call pushit_stack (stack, n_stack, numeric$)
  read (word, *, iostat = ios) stack(n_stack)%value
  if (ios /= 0) then
    err_str = 'BAD NUMERIC: ' // trim(word)
    return
  endif
  stack(n_stack)%name = word

  if (i_op > 0) then
    if (op(i_op) == unary_minus$) then
      stack(n_stack)%name = '-' // stack(n_stack)%name
      stack(n_stack)%value = -stack(n_stack)%value
      i_op = i_op - 1
    endif
  endif

else
  call pushit_stack (stack, n_stack, var_type)
  ! "my_species" in "mass_of(my_species)" is considered a variable and not a species_const
  if (var_type == species_const$) then
    stack(n_stack)%name = word
    stack(n_stack)%value = species_id(word, print_err = .false.)
    if (stack(n_stack)%value == invalid$) stack(n_stack)%type = variable$
  elseif (var_type == variable$) then
    stack(n_stack)%name = upcase(word)
  else
    stack(n_stack)%name = word
  endif

  select case (downcase(word))
  case ('twopi');                 stack(n_stack) = expression_atom_struct(word, constant$, twopi)
  case ('fourpi');                stack(n_stack) = expression_atom_struct(word, constant$, fourpi)
  case ('pi');                    stack(n_stack) = expression_atom_struct(word, constant$, pi)
  case ('e', 'e_log');            stack(n_stack) = expression_atom_struct(word, constant$, exp(1.0_rp))
  case ('sqrt_2');                stack(n_stack) = expression_atom_struct(word, constant$, sqrt_2)
  case ('degrad');                stack(n_stack) = expression_atom_struct(word, constant$, 180 / pi)
  case ('degrees', 'raddeg');     stack(n_stack) = expression_atom_struct(word, constant$, pi / 180)
  case ('m_electron');            stack(n_stack) = expression_atom_struct(word, constant$, m_electron)
  case ('m_muon');                stack(n_stack) = expression_atom_struct(word, constant$, m_muon)
  case ('m_pion_0');              stack(n_stack) = expression_atom_struct(word, constant$, m_pion_0)
  case ('m_pion_charged');        stack(n_stack) = expression_atom_struct(word, constant$, m_pion_charged)
  case ('m_proton');              stack(n_stack) = expression_atom_struct(word, constant$, m_proton)
  case ('m_deuteron');            stack(n_stack) = expression_atom_struct(word, constant$, m_deuteron)
  case ('m_neutron');             stack(n_stack) = expression_atom_struct(word, constant$, m_neutron)
  case ('c_light');               stack(n_stack) = expression_atom_struct(word, constant$, c_light)
  case ('r_e');                   stack(n_stack) = expression_atom_struct(word, constant$, r_e)
  case ('r_p');                   stack(n_stack) = expression_atom_struct(word, constant$, r_p)
  case ('e_charge');              stack(n_stack) = expression_atom_struct(word, constant$, e_charge)
  case ('h_planck');              stack(n_stack) = expression_atom_struct(word, constant$, h_planck)
  case ('h_bar_planck');          stack(n_stack) = expression_atom_struct(word, constant$, h_bar_planck)
  case ('fine_struct_const');     stack(n_stack) = expression_atom_struct(word, constant$, fine_structure_constant)
  case ('anom_moment_electron');  stack(n_stack) = expression_atom_struct(word, constant$, anomalous_mag_moment_electron)
  case ('anom_moment_proton');    stack(n_stack) = expression_atom_struct(word, constant$, anomalous_mag_moment_proton)
  case ('anom_moment_neutron');   stack(n_stack) = expression_atom_struct(word, constant$, anomalous_mag_moment_neutron)
  case ('anom_moment_muon');      stack(n_stack) = expression_atom_struct(word, constant$, anomalous_mag_moment_muon)
  case ('anom_moment_deuteron');  stack(n_stack) = expression_atom_struct(word, constant$, anomalous_mag_moment_deuteron)
  case ('anom_moment_he3');       stack(n_stack) = expression_atom_struct(word, constant$, anomalous_mag_moment_he3)
  end select
endif

err = .false.

end subroutine push_numeric_or_var

!-------------------------------------------------------------------------
! contains

subroutine pushit_stack (stack, ix_stack, this_type)

type (expression_atom_struct), allocatable :: stack(:), temp_stack(:)
integer ix_stack, this_type

!

if (ix_stack == size(stack)) then
  call move_alloc(stack, temp_stack)
  allocate (stack(ix_stack+10))
  stack(1:ix_stack) = temp_stack
endif

call pushit (stack%type, ix_stack, this_type)

end subroutine pushit_stack

end subroutine expression_string_to_stack

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! subroutine pushit (array, ix_arr, this_type)
!
! Private routine used by expression_string_to_stack
!-

subroutine pushit (array, ix_arr, this_type)

integer array(:), ix_arr, this_type
character(*), parameter :: r_name = 'pushit'

!

ix_arr = ix_arr + 1

if (ix_arr > size(array)) then
  call out_io (s_fatal$, r_name, 'ARRAY OVERFLOW, EXPERT HELP IS NEEDED!')
  if (global_com%exit_on_error) call err_exit
  return
endif

array(ix_arr) = this_type

end subroutine pushit

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function expression_value (expression, err_flag, err_str, var, use_old) result (value)
!
! Routine to evaluate a mathematical expression encoded in a string.
!
! Also see:
!   expression_string_to_stack
!   expression_stack_value
!
! Input:
!   expression  -- character(*): Expression string.
!   var(:)      -- control_var1_struct, optional: Array of control variables. 
!                   Used with Bmad controller elements.
!   use_old     -- logical, optional: Use var%old_value? Must be present if var(:) is present.
!
! Output:
!   value       -- real(rp): Value of the expression.
!   err_flag    -- logical: True if there is an evaluation problem. False otherwise.
!   err_str     -- character(*), optional: Error string explaining error if there is one.
!-

function expression_value (expression, err_flag, err_str, var, use_old) result (value)

type (expression_atom_struct), allocatable :: stack(:)
type (control_var1_struct), optional :: var(:)

real(rp) value
integer i, i2, ix, n_stack

logical, optional :: use_old
logical err_flag

character(*) expression
character(*), optional :: err_str
character(100) err_string
character(*), parameter :: r_name = 'expression_value'

!

call expression_string_to_stack (expression, stack, n_stack, err_flag, err_string)
if (err_flag) then
  if (present(err_str)) then
    err_str = err_string
  else
    call out_io (s_error$, r_name, err_string, 'FOR EXPRESSION: ' // expression)
  endif
  value = real_garbage$
  return
endif

value = expression_stack_value (stack, err_flag, err_str, var, use_old) 
if (err_flag) then
  if (present(err_str)) then
    err_str = err_string
  else
    call out_io (s_error$, r_name, err_string, 'FOR EXPRESSION: ' // expression)
  endif
endif

end function expression_value

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function expression_stack_value (stack, err_flag, err_str, var, use_old) result (value)
!
! Routine to evaluate a mathematical expression represented by an "expression stack".
! Expression stacks are created by expression_string_to_stack.
!
! Note: Stack elements with stack(i)%type == variable$ need to be evalauated before
! calling this routine and the value placed in stack(i)%value.
!
! Also see:
!   expression_value
!   expression_string_to_stack
!
! Input:
!   stack(:)    -- expression_atom_struct: Expression to evaluate.
!   var(:)      -- control_var1_struct, optional: Array of control variables. 
!                   Used with Bmad controller elements.
!   use_old     -- logical, optional: Use var%old_value? Must be present if var(:) is present.
!
! Output:
!   value       -- real(rp): Value of the expression.
!   err_flag    -- logical: True if there is an evaluation problem. False otherwise.
!   err_str     -- character(*): Error string explaining error if there is one.
!-

function expression_stack_value (stack, err_flag, err_str, var, use_old) result (value)

use random_mod

type (expression_atom_struct) stack(:)
type (expression_atom_struct) stack2(size(stack))
type (control_var1_struct), optional :: var(:)

real(rp) value
integer i, i2, ix

logical, optional :: use_old
logical err_flag

character(*) err_str
character(*), parameter :: r_name = 'expression_stack_value'

!

err_flag = .true.

i2 = 0
do i = 1, size(stack)

  select case (stack(i)%type)

  case (end_stack$)
    exit

  case (arg_count$)
    cycle

  case (numeric$, variable$, constant$)
    i2 = i2 + 1
    stack2(i2) = stack(i)

  case (unary_minus$)
    stack2(i2)%value = -stack2(i2)%value

  case (unary_plus$)
    stack2(i2)%value = stack2(i2)%value

  case (plus$)
    stack2(i2-1)%value = stack2(i2-1)%value + stack2(i2)%value
    i2 = i2 - 1

  case (minus$)
    stack2(i2-1)%value = stack2(i2-1)%value - stack2(i2)%value
    i2 = i2 - 1

  case (times$)
    stack2(i2-1)%value = stack2(i2-1)%value * stack2(i2)%value
    i2 = i2 - 1

  case (divide$)
    if (stack2(i2)%value == 0) then
      err_str = 'DIVIDE BY 0 IN EXPRESSION'
      return
    endif
    stack2(i2-1)%value= stack2(i2-1)%value / stack2(i2)%value
    i2 = i2 - 1

  case (power$)
    stack2(i2-1)%value = stack2(i2-1)%value**stack2(i2)%value
    i2 = i2 - 1

  case (cot$)
    if (stack2(i2)%value == 0) then
      err_str = 'DIVIDE BY 0 IN EXPRESSION'
      return
    endif
    stack2(i2)%value = 1.0_rp / tan(stack2(i2)%value)

  case (csc$)
    if (stack2(i2)%value == 0) then
      err_str = 'DIVIDE BY 0 IN EXPRESSION'
      return
    endif
    stack2(i2)%value = 1.0_rp / sin(stack2(i2)%value)

  case (sec$)
    stack2(i2)%value = 1.0_rp / cos(stack2(i2)%value)

  case (sin$)
    stack2(i2)%value = sin(stack2(i2)%value)

  case (sinc$)
    stack2(i2)%value = sinc(stack2(i2)%value)

  case (cos$)
    stack2(i2)%value = cos(stack2(i2)%value)

  case (tan$)
    stack2(i2)%value = tan(stack2(i2)%value)

  case (asin$)
    if (stack2(i2)%value < -1 .or. stack2(i2)%value > 1) then
      err_str = 'ASIN ARGUMENT HAS MAGNITUDE GREATER THAN 1'
      return
    endif
    stack2(i2)%value = asin(stack2(i2)%value)

  case (acos$)
    if (stack2(i2)%value < -1 .or. stack2(i2)%value > 1) then
      err_str = 'ACOS ARGUMENT HAS MAGNITUDE GREATER THAN 1'
      return
    endif
    stack2(i2)%value = acos(stack2(i2)%value)

  case (factorial$)
    stack2(i2)%value = factorial(nint(stack2(i2)%value))
    if (stack2(i2)%value < 0) then
      err_str = 'FACTORIAL PROBLEM'
      return
    endif

  case (atan$)
    stack2(i2)%value = atan(stack2(i2)%value)

  case (atan2$)
    stack2(i2-1)%value = atan2(stack2(i2-1)%value, stack2(i2)%value)
    i2 = i2 - 1

  case (modulo$)
    stack2(i2-1)%value = modulo(stack2(i2-1)%value, stack2(i2)%value)
    i2 = i2 - 1

  case (sinh$)
    stack2(i2)%value = sinh(stack2(i2)%value)

  case (cosh$)
    stack2(i2)%value = cosh(stack2(i2)%value)

  case (tanh$)
    stack2(i2)%value = tanh(stack2(i2)%value)

  case (coth$)
    if (stack2(i2)%value == 0) then
      err_str = 'DIVIDE BY 0 IN EXPRESSION'
      return
    endif
    stack2(i2)%value = 1.0_rp / tanh(stack2(i2)%value)

  case (asinh$)
    stack2(i2)%value = asinh(stack2(i2)%value)

  case (acosh$)
    stack2(i2)%value = acosh(stack2(i2)%value)

  case (atanh$)
    stack2(i2)%value = atanh(stack2(i2)%value)

  case (acoth$)
    if (stack2(i2)%value == 0) then
      err_str = 'DIVIDE BY 0 IN EXPRESSION'
      return
    endif
    stack2(i2)%value = 1.0_rp / atanh(stack2(i2)%value)

  case (abs$)
    stack2(i2)%value = abs(stack2(i2)%value)

  case (sqrt$)
    if (stack2(i2)%value < 0) then
      err_str = 'SQRT ARGUMENT IS NEGATIVE '
      return
    endif
    stack2(i2)%value = sqrt(stack2(i2)%value)

  case (log$)
    if (stack2(i2)%value < 0) then
      err_str = 'LOG ARGUMENT IS NEGATIVE '
      return
    endif
    stack2(i2)%value = log(stack2(i2)%value)

  case (exp$)
    stack2(i2)%value = exp(stack2(i2)%value)

  case (int$)
    stack2(i2)%value = int(stack2(i2)%value)

  case (sign$)
    stack2(i2)%value = sign_of(stack2(i2)%value)

  case (nint$)
    stack2(i2)%value = nint(stack2(i2)%value)

  case (floor$)
    stack2(i2)%value = floor(stack2(i2)%value)

  case (ceiling$)
    stack2(i2)%value = ceiling(stack2(i2)%value)

  case (mass_of$)
    stack2(i2)%value = mass_of(nint(stack2(i2)%value))

  case (charge_of$)
    stack2(i2)%value = charge_of(nint(stack2(i2)%value))

  case (antiparticle$)
    stack2(i2)%value = antiparticle(nint(stack2(i2)%value))

  case (anomalous_moment_of$)
    stack2(i2)%value = anomalous_moment_of(nint(stack2(i2)%value))

  case (ran$)
    i2 = i2 + 1
    call ran_uniform(stack2(i2)%value)

  case (ran_gauss$)
    if (nint(stack(i-1)%value) == 0) then
      i2 = i2 + 1
      call ran_gauss(stack2(i2)%value)
    else
      call ran_gauss(value, sigma_cut = stack2(i2)%value)
      stack2(i2)%value = value
    endif

  case (species_const$)
    i2 = i2 + 1
    stack2(i2)%value = stack(i)%value

  case (species$)
    ! Nothing to do

  case default

    if (is_attribute(stack(i)%type, control_var$)) then
      if (.not. present(var)) then
        err_flag = .true.
        err_str = 'VAR ARGUMENT NOT PRESENT! GET HELP!'
        return
      endif
      ix = stack(i)%type - var_offset$
      i2 = i2 + 1
      if (logic_option(.false., use_old)) then
        stack2(i2)%value = var(ix)%old_value
      else
        stack2(i2)%value = var(ix)%value
      endif

    else
      err_str = 'INTERNAL ERROR #02: GET HELP'
      if (global_com%exit_on_error) then
        call out_io (s_fatal$, r_name, err_str)
        call err_exit
      endif
      return
    endif

  end select
enddo

if (i2 /= 1) then
  err_str = 'INTERNAL ERROR #03: GET HELP'
  if (global_com%exit_on_error) then
    call out_io (s_fatal$, r_name, err_str)
    call err_exit
  endif
  return
endif

value = stack2(1)%value
err_flag = .false.

end function expression_stack_value

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function expression_stack_to_string (stack, polish) result (str)
!
! Routine to convert an expression stack to a string
!
! Input:
!   stack(:)  -- expression_atom_struct: arithmetic expression
!   polish    -- logical, optional, Construct expression in reverse polish? Default is False.
!
! Output:
!   str       -- character(:), allocatable: Expression in string form.
!-

function expression_stack_to_string (stack, polish) result (str)

type (expression_atom_struct), target :: stack(:)
type (expression_atom_struct), pointer :: atom, atom2
type (expression_atom_struct) s2(size(stack))
type (var_length_string_struct) s2_name(size(stack))

character(:), allocatable :: str

integer i, i2, ix
logical, optional :: polish

!

allocate (character(1) :: str)
str = ''

! Polish notation

if (logic_option(.false., polish)) then
  do i = 1, size(stack)
    atom => stack(i)
    if (atom%type == end_stack$) exit
    if (atom%type >= 1 .and. atom%type <= size(expression_op_name)) then
      str = trim(str) // ', ' // expression_op_name(atom%type)    
    else
      str = trim(str) // ', ' // atom%name
    endif
  enddo

  str = str(3:)

! Standard notation

else

  i2 = 0
  i = 0

  do
    if (i+1 > size(stack)) exit
    if (stack(i+1)%type == end_stack$) exit

    i = i + 1
    atom => stack(i)

    if (is_attribute(atom%type, all_control_var$)) then
      i2 = i2 + 1
      s2(i2)%type = variable$

      if (atom%name == '') then
        s2_name(i2)%str = real_to_string(atom%value, 20, 14)
      else
        s2_name(i2)%str = trim(atom%name)
      endif
      cycle
    endif

    select case (atom%type)
    case (plus$, minus$, times$, divide$, power$)
      if (expression_eval_level(s2(i2-1)%type) < expression_eval_level(atom%type)) s2_name(i2-1)%str = '(' // trim(s2_name(i2-1)%str) // ')'
      if (atom%type == minus$ .or. atom%type == divide$) then
        if (expression_eval_level(s2(i2)%type) <= expression_eval_level(atom%type)) s2_name(i2)%str = '(' // trim(s2_name(i2)%str) // ')'
      else
        if (expression_eval_level(s2(i2)%type) < expression_eval_level(atom%type)) s2_name(i2)%str = '(' // trim(s2_name(i2)%str) // ')'
      endif
      s2_name(i2-1)%str = trim(s2_name(i2-1)%str) // trim(expression_op_name(atom%type)) // s2_name(i2)%str
      s2(i2-1)%type = atom%type
      i2 = i2 - 1

    case (numeric$, variable$, constant$, species_const$)
      i2 = i2 + 1
      s2(i2)%type = atom%type
      if (atom%name == '') then
        s2_name(i2)%str = real_to_string(atom%value, 20, 14)
      else
        s2_name(i2)%str = atom%name
      endif

    case (unary_minus$)
      if (expression_eval_level(s2(i2)%type) <= expression_eval_level(atom%type)) s2_name(i2)%str = '(' // trim(s2_name(i2)%str) // ')'
      s2_name(i2)%str = '-' // s2_name(i2)%str
 
    case (unary_plus$)
      if (expression_eval_level(s2(i2)%type) <= expression_eval_level(atom%type)) s2_name(i2)%str = '(' // trim(s2_name(i2)%str) // ')'
      s2_name(i2)%str = '+' // s2_name(i2)%str
 
    case (ran$)
      i2 = i2 + 1
      s2_name(i2)%str = trim(expression_op_name(atom%type)) // '()'
      s2%type = atom%type

    case (ran_gauss$)
      if (nint(stack(i-1)%value) == 0) then
        i2 = i2 + 1
        s2_name(i2)%str = trim(expression_op_name(atom%type)) // '()'
      else
        s2_name(i2)%str = trim(expression_op_name(atom%type)) // '(' // trim(s2_name(i2)%str) // ')'
      endif
      s2%type = atom%type

    case (atan2$, modulo$)
      i2 = i2 - 1
      s2_name(i2)%str = trim(expression_op_name(atom%type)) // '(' // trim(s2_name(i2)%str) // ',' // trim(s2_name(i2+1)%str) // ')'
      s2%type = atom%type

    case (factorial$)
      if (expression_eval_level(s2(i2)%type) <= expression_eval_level(atom%type)) s2_name(i2)%str = '(' // trim(s2_name(i2)%str) // ')'
      s2_name(i2)%str = trim(s2_name(i2)%str) // '!'
      s2%type = atom%type

    case (arg_count$)
      cycle

    case default ! Function
      s2_name(i2)%str = trim(expression_op_name(atom%type)) // '(' // trim(s2_name(i2)%str) // ')'
      s2%type = atom%type

    end select

  enddo

  str = s2_name(i2)%str
endif

end function expression_stack_to_string

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine split_expression_string (expr, width, indent, lines)
!
! Routine to break an expression into a number of lines for a nicer display.
! Used when printing expressions.
!
! Input:
!   expr      -- character(*): String containing the expression.
!   width     -- integer: Maximum width of split expression.
!   indent    -- integer: If positive: Number of spaces to indent for every line after the first.
!                         If negative: No indentation but first line is shortened by |indent|.
!   break_str -- character(*), optional: If present, only break lines at places where this string is.
!
! Output:
!   lines(:)  -- character(*), allocatable: Split expression.
!-

subroutine split_expression_string (expr, width, indent, lines, break_str)

integer width, indent
integer nn, n0, nl, ind, ww, i_split, large, i, j

real(rp), parameter :: weight(11) = [1.0_rp, 1.0_rp, 1.0_rp, 0.95_rp, 0.95_rp, 0.8_rp, 0.8_rp, 0.8_rp, 0.8_rp, 0.8_rp, 0.8_rp]
real(rp) score

character(*) expr
character(*), allocatable :: lines(:)
character(*), optional :: break_str
character(len(expr)) ex
character(width), allocatable :: li(:)
character(1) ch, ch0
character(11), parameter :: ch_split = '+-,*/)]}([{'

!

ind = max(0, -indent) ! Indent for first line.
nn = len_trim(expr)
n0 = 1 + 2*nn/(width-abs(indent))
call re_allocate (li, n0, .true., '')

ex = expr
nl = 0
do
  nl = nl + 1
  nn = len_trim(ex)
  ww = width - ind
  if (nn <= ww) then
    if (nl == 1) then
      li(nl) = ex
    else
      li(nl)(ind+1:) = ex
    endif
    exit
  endif

  i_split = ww
  if (present(break_str)) then
    nn = len(break_str)
    do j = ww, 2, -1
      if (ex(j-nn+1:j) /= break_str) cycle
      i_split = j
      exit
    enddo

  else
    score = 0
    do j = ww, 2, -1
      ch = ex(j:j)
      ch0 = ex(j-1:j-1)
      i = index(ch_split, ch)
      if (i == 0) cycle
      if ((ch == '-' .or. ch == '+') .and. (ch0 == 'e' .or. ch0 == 'E' .or. index(ch_split, ch0) /= 0)) cycle
      if (j * weight(i) < score) cycle
      i_split = j
      score = j * weight(i)
    enddo
  endif

  if (nl == 1) then
    li(nl) = ex(1:i_split)
  else
    li(nl)(ind+1:) = ex(1:i_split)
  endif
  ind = max(0, indent)  ! Indent for lines after first.

  ex = adjustl(ex(i_split+1:))
enddo

call re_allocate (lines, nl)
do i = 1, nl
  lines(i) = li(i)
enddo

end subroutine split_expression_string

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function linear_coef (stack, err_flag) result (coef)
!
! Routine to return the linear coefficient of a linear expression.
!
! Input:
!   stack(:) -- expression_atom_struct: Expression stack.
!
! Output:
!   err_flag  -- Logical: Set True if the expression is not linear
!   coef      -- real(rp): Linear coefficient.
!-

function linear_coef (stack, err_flag) result (coef)

type (expression_atom_struct) stack(:)

real(rp) coef
logical :: err_flag

integer i, i0, i1, n

character(40) err_str

!

n = size(stack)
err_flag = .true.
coef = 0

! If expression = "Var" then coef = 1.

if (n == 1 .and. is_attribute(stack(1)%type, control_var$)) then
  err_flag = .false.
  coef = 1
  return
endif

! Expression must have times at the end to be linear.

if (stack(n)%type /= times$) return

! To be linear, stack(1) or stack(n-1) must be variable and have no other variables.

if (is_attribute(stack(1)%type, control_var$)) then
  i0 = 2
  i1 = n-1
elseif (is_attribute(stack(n-1)%type, control_var$)) then
  i0 = 1
  i1 = n-2
else
  return
endif

coef = expression_stack_value(stack(i0:i1), err_flag, err_str)

end function linear_coef

end module
