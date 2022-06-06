!input functions
MODULE infuncs
  USE globals
  IMPLICIT NONE
  PRIVATE
  PUBLIC get_cmdargs

CONTAINS

  SUBROUTINE get_cmdargs()
    CHARACTER(100) :: t_char
    CALL GET_COMMAND_ARGUMENT(1,t_char)
    READ(t_char,*)num_customers
  ENDSUBROUTINE get_cmdargs
END MODULE infuncs
