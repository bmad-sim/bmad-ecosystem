! *******************************************************************
! *    Function to compare two strings over the first min_char      *
! *    number of characters and return .TRUE. if there is a match   *
! *******************************************************************
logical function string_comp(str1, str2, in_min_char)
  implicit none


  character*(*)::  str1
  character*(*)::  str2
  character,save:: string1*60, string2*60
  !   ----------------end_declare---------------------------
  integer :: in_min_char
  integer :: lenstr
  integer, save:: min_char, len1, len2, min_len
  string_comp = .FALSE.  !assume this until proven otherwise
  string1 = str1
  string2 = str2
  call str_upcase(string1, string1)
  call str_upcase(string2, string2)
  len1 = lenstr(string1)
  len2 = lenstr(string2)
  min_len = min(len1, len2)
  min_char = max(min_len, in_min_char) !compare all characters in short
  if( string1(1:min_char).eq.string2(1:min_char) ) string_comp = .TRUE.
  return
end function string_comp
