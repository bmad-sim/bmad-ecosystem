        program test_find_first
	use sim_utils
        use str_find_first_substring_module
        implicit none
        logical match
        integer which,where
        character*40 str,sub1,sub2,sub3
        str=' Hello world testing s3 '
        sub1='s1' ; sub2='s2' ; sub3='s3'
        match=str_find_first_substring(str,where,which,sub1,sub2,sub3)
        print *,which,where
	stop
	end program
