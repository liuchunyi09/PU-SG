
	CCD++ algorithm for PU-SG model

=======================================================================================================================
	Programmer:	 Chun-Yi Liu

	Summary: 	 The program solves the PU-SG model for event recommendation in Event-Based Social Networks. Users can 
		input the user-event pairs, user-group pairs, users' location and events' location by four .txt file. And by 
		setting the parameter in main.app and Limits.h, users can test the performance in different parameters. The 
		program measures the four model i.e. PU-SG, PU-G, PU-S and PU-0 by AUC and compares our proposed model PU-SG with
		baseline methods i.e. Random, Most Popular, Location-Aware, User-KNN and Matrix Factorization by measuring P@3, 
		P@5 and MAP. Users will obtain a .txt files contain the results.

	Input:		 What is need are four .txt files.

			"Hou_ug.txt" : the file contains the user-group pairs. In this file, each line contain two integers separated 
				by space. The two integers are the user id and group id, e.g. 1 2 is the user with id '1' join the group
				with id '2'. NOTICE the user id number cannot be larger than the total user number, group is the same.

			"Hou_ue.txt" : this file is similar with "Hou_ug.txt" by replacing the group id to event id. NOTICE that the 
				event id number cannot be larger than the total event number, too.

			"Hou_ll_user.txt" : this file contains users' location information demonstrated by longitude and latitude. In
				this file, each line are two floating number separated by space. The location information corresponds to the
				user id, i.e, the first line records the longitude and latitude of the first user. NOTICE that we exchange 
				east longitude to positive value while the west longitude to negative value. And we exchange the south 
				latitude to '90 + south latitude' while the north latitude to '90 - latitude'. 

			"Hou_ll_event.txt" : this file is the location of event which plays the same role in "Hou_ll_user.txt".

	Output:		A .txt file named "Hou.txt" contains the results.
					The first four lines is the AUC of PU-SG, PU-G, PU-S, PU-0 in order.  The P@3, p@5 and MAP are clearly
					shown in the file.

	Parameter Control:
		In 'Limit.h', 
			'RATIO' is the ratio of test data in total data;
			'ACCU'	is the accuracy in stopping criterion;
			'ITMAX' is the max iteration times;

		In 'main.cpp',
			'alphax', 'lambda', 'rho_on' and 'rho_off' are corresponding to the parameter in PU-SG model.
			'disbest' is the dimension of latent factor.

