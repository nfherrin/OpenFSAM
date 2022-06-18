!input functions
MODULE infuncs
  USE globals
  IMPLICIT NONE
  PRIVATE
  PUBLIC get_cmdargs

CONTAINS

  SUBROUTINE get_cmdargs()
    CHARACTER(100) :: t_char
    INTEGER :: arg_count
    arg_count=COMMAND_ARGUMENT_COUNT()
    IF(arg_count .GE. 1)THEN
      CALL GET_COMMAND_ARGUMENT(1,t_char)
      READ(t_char,*)num_customers
    ELSE
      num_customers=10
    ENDIF
  ENDSUBROUTINE get_cmdargs
END MODULE infuncs
