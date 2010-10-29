module column_mod

type column_struct 
  real row(15)
end type

integer, private :: private_dummy ! This is to suppress the ranlib "has no symbols" message

end module
