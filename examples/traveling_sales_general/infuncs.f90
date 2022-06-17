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
    IF(arg_count .GE. 2)THEN
      CALL GET_COMMAND_ARGUMENT(1,t_char)
      READ(t_char,*)prob_dim
      CALL GET_COMMAND_ARGUMENT(2,t_char)
      READ(t_char,*)num_customers
    ELSEIF(arg_count .GE. 1)THEN
      CALL GET_COMMAND_ARGUMENT(1,t_char)
      READ(t_char,*)prob_dim
      num_customers=10
    ELSE
      prob_dim=1
      num_customers=30
    ENDIF
  ENDSUBROUTINE get_cmdargs
END MODULE infuncs
