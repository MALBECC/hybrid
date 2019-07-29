	module precision
!old program presicions
	integer, parameter :: sp = selected_real_kind(6,30)
	integer, parameter :: dp = selected_real_kind(14,100)
	integer, parameter :: qp = selected_real_kind(33, 4931)

!this will be the program presicion in future modifications
#ifdef SPRESS
	integer, parameter :: xp = selected_real_kind(6, 37)
#elif DPRESS
	integer, parameter :: xp = selected_real_kind(15, 307)
#elif QPRESS
	integer, parameter :: xp = selected_real_kind(33, 4931)
#endif
	end module precision

