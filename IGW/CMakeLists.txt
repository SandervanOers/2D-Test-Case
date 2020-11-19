include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)

#First declare the rules that will form the executable tests#
#############################################################


if(hpGEM_INCLUDE_PETSC_SUPPORT)

	configure_file(bCoeff.txt
			${PROJECT_BINARY_DIR}/bCoeff.txt @ONLY IMMEDIATE)
    include_directories(${PETSC_INCLUDES})

	add_executable(IGW.out
        main.cpp
        HEuler.cpp
        InitialConditions.cpp
	poly_eval.cpp
	airy.cpp
        )
	
		target_link_libraries(IGW.out Output ${PETSC_LIBRARIES})
else() 
	add_custom_target(IGW.out 
		COMMENT "\t ERROR:: You need PETSC to compile IGW.out. \n\t ERROR:: Please rerun cmake and enable PETSC support \n\t ERROR:: If using raw cmake add the flag -DhpGEM_INCLUDE_PETSC_SUPPORT='on'")
endif()

