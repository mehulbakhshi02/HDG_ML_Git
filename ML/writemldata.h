/** @file
 * @brief Computes the L2 error of the solution or the error in target functional based on the solution that is there in the
 * global array ed. Should be modified to incorporate the solution at the current
 * state. Should also contain the details of the adjoint based error estimate and
 * error bounds. We decrese the order of the solution to work at the order p
*/

template<int D, int COMP, class Model>
void UnifyingFramework<D, COMP, Model>
::WriteMLData(const Solution & new_sol, LocalHeap & lh) {

	Solution adj;
	vector<double> error(ne);
	int its;
	if(AnisotropyData::adjoint_based)
  	{
		SolveAdjointSystem(new_sol, adj, error, its, lh);
	}

	//Writing error estimator csv file
	stringstream err_oss1(" ");
  	err_oss1 << "cell_wise_adj_err.csv";
  	ofstream err_log1(err_oss1.str().c_str(), ios::app);
  	err_log1.precision(16);
  	err_log1.setf(ios::scientific, ios::floatfield);
  
  	for(int i=0;i<ne;i++)
  	{
   		ML_ElementData<D, COMP> & ml_elidata = *ml_eldata[i];
    	// if(ml_elidata.nf!=ml_elidata.nf_bc) continue;

    	double normalized_error = log(error[i]);
		err_log1 << error[i] << "," << normalized_error << endl;
 	}

  	err_log1.close();
	cout<< "cell_wise_adj_err done" << endl;

	  // Writing the error monitors
	stringstream err_oss(" ");
    // oss << "solution-" << ne << "-" << p;
	err_oss << "mldata-" << ne << "-" << max_order;
    Model::GetFilename(err_oss);
  	err_oss << ".txt";
	ofstream err_log(err_oss.str().c_str(), ios::app);
	err_log.precision(16);
	int ne_int = 0;
	for(int i=0;i<ne;i++)
	{
		ML_ElementData<D, COMP> & ml_elidata = *ml_eldata[i];
		
		// KUNAL
		//if(ml_elidata.nf==ml_elidata.nf_bc)
				ne_int++;
	}
	
	Model::GetFilename(err_oss);
	Vec<Model::NumParam> param;
	Model::GetParameters(param);
	//err_log << "starting line 1" << endl;
	err_log << ne_int << "," << D << "," << COMP <<endl;
	
	// KUNAL
	/*
	for(int i=0;i<ne;i++)
	{
		// cout<<"Number of elements: "<<ne<<endl;
		// cout<<"ML error for the element: "<<i<<" : "<<ml_err[i]<<endl;

		ML_ElementData<D, COMP> & ml_elidata = *ml_eldata[i];

		for (int j = 0; j < ml_elidata.nf; ++j)
		{
		cout << "fcnr:" << ml_elidata.faces[j];
		}
		cout <<endl;

	}
	*/
	
	// KUNAL: I need the nip to know how many values are missing
	// To store the nip values at the faces
	int nip_universal;
	// A marker to check that the value is acquired
	int nip_acquired = 0; 
	// To debug the smallest 2 elements case
	// As a 6 faces are there in total for 2 triangular element
	int ctr = 0;
	for(int i = 0; i < ne && nip_acquired == 0; i++)
	{
		ML_ElementData<D, COMP> & ml_elidata = *ml_eldata[i];
		ElementData<D, COMP> & elidata = *eldata[i];
		for(int j = 0 ; j < ml_elidata.nf ; j++)
		{
			int fcnr = ml_elidata.faces[j];
			// vect[ctr] = fcnr;
			// int fcnr = elidata.faces[j];
			if(ml_elidata.nf==ml_elidata.nf_bc)
			{
				ML_FacetData<D, COMP> & ml_fd_1 = *ml_fadata[fcnr];
				nip_universal = ml_fd_1.nip;
				nip_acquired = 1; // Value acquired
			}
			// ctr = ctr + 1;
			// To fix the 2 elements case
			
		}
	}
	
	// To fix the 2 elements case
	if (ne == 2)
	{
		vector<int> vect(6);
		for(int i = 0; i < ne && nip_acquired == 0; i++)
		{
			ML_ElementData<D, COMP> & ml_elidata = *ml_eldata[i];
			for(int j = 0 ; j < ml_elidata.nf ; j++)
			{
				int fcnr = ml_elidata.faces[j];
				vect[ctr] = fcnr;
				ctr = ctr + 1;			
			}
		}
		// Check the first element
		for (int i = 0;i < 3;i++)
		{
			// Check the second element
			for (int j = 3;j < vect.size();j++)
			{
				if (vect[i] == vect[j])
				{
					int fcnr = vect[i];
					cout<< fcnr<<endl;
					ML_FacetData<D, COMP> & ml_fd_1 = *ml_fadata[fcnr];
					nip_universal = ml_fd_1.nip;
					nip_acquired = 1;
				}
			}
		}
	}
	
	// cout<<"nip_universal: "<<nip_universal<<endl;
	
	// KUNAL
	// cout<< "Total Elements: "<<ne<<endl;
	for(int i = 0; i < ne; i++)
	{
		// KUNAL
		// cout<< "Element Number: " << i <<endl;
		ML_ElementData<D, COMP> & ml_elidata = *ml_eldata[i];
		
		//
		ElementData<D, COMP> & elidata = *eldata[i];
		//
		
		// KUNAL: Commenting the next line
		// if(ml_elidata.nf!=ml_elidata.nf_bc) continue;
		
		
		// KUNAL: Trying to create a condition
		if(ml_elidata.nf!=ml_elidata.nf_bc)  // if it is a boundary element, then goes in
		{
			// For the boundary elements
			// cout<< "Boundary Elements Activate" << endl;
			
			int ndof_w = ml_elidata.ndof_w;
			// KUNAL
			//cout<< "ndof_w: "<< ndof_w<<endl;
			// Taking the quadrature nodes from the first face
			int fcnr = ml_elidata.faces[0];
			// int fcnr = ml_elidata.faces[0];
			ML_FacetData<D, COMP> & ml_fd_1 = *ml_fadata[fcnr];
			
			//  Original 
			// int nip_1 = ml_fd_1.nip;
			int nip_1 = 0;
			
			// We need to give the value of actual nip to use the original setup
			// We will use the nip_universal
			
			nip_1 = nip_universal;
						
			// cout << "nip_1: "<< nip_1<<endl;
			
			// Graph Neural Network
			// err_log << i <<","<< ml_elidata.nip << "," << ml_elidata.ndof_w*(D+1) + nip_1*(D+1)*ml_elidata.nf + (D*(D+1))/2 + Model::NumParam << "," << ml_elidata.nf_bc << endl;
			
			// Feedforward neural network
			//err_log << "starting line 2 (boundary element)" << endl;
			err_log << "1" << "," << (ml_elidata.ndof_w*(D+1) + nip_1*(D+1)*ml_elidata.nf)*COMP + (D*(D+1))/2 + Model::NumParam << "," << ml_elidata.nf_bc << endl;
			// cout<< "ml_elidata.ndof_w*(D+1): "<< ml_elidata.ndof_w*(D+1) <<endl;
			
			// This is the problem
			// cout<< "nip_1: "<< nip_1<<endl;
			// cout<< "ml_elidata.nf: "<< ml_elidata.nf <<endl;
			// cout<< "nip_1*(D+1)*ml_elidata.nf: "<< nip_1*(D+1)*ml_elidata.nf <<endl;
			
			
			// cout<< "(D*(D+1))/2: "<< (D*(D+1))/2 <<endl;
			// cout<< "Model::NumParam: "<< ml_elidata.ndof_w*(D+1) <<endl;
			//err_log << "volume, beta, theta (boundary element)" << endl;
			err_log << ml_elidata.vol << endl << ml_elidata.beta << endl << ml_elidata.theta << endl;
			
			//err_log << "parameters" << endl;
			for(int pr = 0; pr<Model::NumParam;pr++)
			{
				err_log << param(pr) <<endl;
	    	}
	    		
	    		
			// Generating the input data for solution and gradient
			for (int ll = 0; ll < COMP; ++ll)
			{
				//err_log << "Solution (boundary element)" << endl;
				for (int pp = 0; pp < ndof_w; ++pp)
				{
					err_log << ml_elidata.vecW[pp+ll*ndof_w] << endl; 
					//cout << "ml_elidata.vecW[pp+ll*ndof_w]: " << ml_elidata.vecW[pp+ll*ndof_w] << endl;
				}
			}

			// KUNAL
			// cout<< "Solution Printed" << endl;

			// Generating the input data for solution and gradient
			if (Model::Diffusion || Model::Source){
			for (int ll = 0; ll < COMP; ++ll)
			{
				//err_log << "Solution gradient (boundary element)" << endl;
				for (int pp = 0; pp < ndof_w; ++pp)
				{
					for(int dd = 0; dd < D; ++dd)
					{
						err_log << ml_elidata.vecQ[pp+ndof_w*(dd+D*ll)] << endl;
						//cout << "ml_elidata.vecQ[pp+ndof_w*(dd+D*ll)]: "<<ml_elidata.vecQ[pp+ndof_w*(dd+D*ll)] << endl;
					}
				}
			}}

			// KUNAL
			// cout<< "Solution gradient printed" << endl;

			// KUNAL
			// cout<<"ml_elidata.nf: "<<ml_elidata.nf<<endl;
			
			
			// KUNAL: Jump data starts from here
			// cout<<"Jump Data: "<< endl; 
			for(int j = 0 ; j < ml_elidata.nf ; j ++)
			{

				int fcnr = ml_elidata.faces[j];
				/*
				ML_FacetData<D, COMP> & ml_fd = *ml_fadata[fcnr];
				int nip = ml_fd.nip;
				*/
				
				if (!fadata[fcnr]) continue;
				FacetData<D, COMP> & fd = *fadata[fcnr];
				// KUNAL: Why they are doing this?
				// They have used a condition to only pass the internal elements
				// But without this the code is crashing
				// Original
				// if (fd.ndof_l == 0) continue; // skip boundary faces
				
				if (fd.bcnr == -2) continue; // zero-measure face
				
				// KUNAL
				if (fd.ndof_l == 0)
				{	
					// Boundary faces of the boundary element
					// Jump in w on quarature points
					// cout<< "Jump Data at the boundary: w: Printing" <<endl;
					// Original
					// for (int qp = 0; qp < nip; ++qp) 
					for (int qp = 0; qp < nip_universal; ++qp) 
					{
						// It does NOT exist at the boundary
						// Commenting
						// double *jump_w = &ml_fd.w_jump[qp*COMP];
						//err_log << "Solution Jump at boundary faces (boundary element)" << endl;
						for(int ll=0;ll<COMP;ll++)
						{
							
							// Testing  KUNAL
							FacetData<D, COMP> & fd = *fadata[fcnr];
							double *jump_w = &fd.lambda[qp*COMP];
							/*
							cout<< "Some w data accessed" << endl;
							cout << jump_w[ll] <<endl;
							// cout << bc_type << endl;// Giving 0 for
							if (jump_w[ll] != 0.0)
							{
								cout<< "Non-zero value"<<endl;
							}*/
							
							// Original
							err_log << jump_w[ll] <<endl;
							//cout<< "jump_w[ll]: "<< jump_w[ll] <<endl;
							
							// Marker
							// err_log << "Element Number: " << i<<", "<< "Face number: "<< fcnr << "Value required: w" <<endl;
						}
					}
					
					// Jump in q on quadrature points
					// cout<< "Jump Data at the boundary: q: Printing" <<endl;
					// Original
					// for (int qp = 0; qp < nip; ++qp) 
					if (Model::Diffusion || Model::Source){
					//err_log << "Solution gradient jump at boundary face (boundary element)" << endl;
					for (int qp = 0; qp < nip_universal; ++qp) 
					{
						// It does NOT exist at the boundary
						// Commenting
						// double *jump_q = &ml_fd.q_jump[qp * COMP * D];

						for(int ll=0;ll<COMP;ll++)
						{

							for(int dd=0;dd<D;dd++)
							{
								// Original
								// err_log << jump_q[dd+D*ll] << endl;
								err_log << 0.0 << endl;
								// Marker
								// err_log << "Element Number: " << i<<", "<< "Face number: "<< fcnr << ", "<< "Dimension Number: "<< dd << ", "<< "Value required: q" <<endl;
							}
						}

					}
					}
				}
				else
				{
					//nip_1 Internal face of the boundary element 
					ML_FacetData<D, COMP> & ml_fd = *ml_fadata[fcnr];
					int nip = ml_fd.nip;

					// Jump in w on quarature points
					// cout<< "Jump Data: w: Printing" <<endl;
					//err_log << "Solution Jump at internal face (boundary element)" << endl;
					for (int qp = 0; qp < nip; ++qp) 
					{
						double *jump_w = &ml_fd.w_jump[qp*COMP];
						for(int ll=0;ll<COMP;ll++)
						{
							err_log << jump_w[ll] <<endl;
							
							//cout<< "jump_w[ll]: "<< jump_w[ll] <<endl;
						}
					}
					
					// Jump in q on quadrature points
					// cout<< "Jump Data: q: Printing" <<endl;
					if (Model::Diffusion || Model::Source){
					//err_log << "Solution gradient Jump at internal face (boundary element)" << endl;
					for (int qp = 0; qp < nip; ++qp) 
					{
						double *jump_q = &ml_fd.q_jump[qp * COMP * D];

						for(int ll=0;ll<COMP;ll++)
						{

							for(int dd=0;dd<D;dd++)
							{
								err_log << jump_q[dd+D*ll] << endl;
								//cout << "jump_q[dd+D*ll]: "<<jump_q[dd+D*ll] << endl;
							}
						}

					}
				}
				}

			}
			
			// KUNAL
			// Boundary Element ML output
			//err_log << "ML Output Data (Boundary element)" << endl;
			//ML_OutputData<D, COMP> & ml_outdata = *ml_outputdata[i];
			//for (int j = 0; j < ml_outdata.nip; ++j) 
			//{
			// const double * pw_out  = &ml_outdata.w[j * COMP];
			//	for (int ll = 0; ll < COMP; ++ll)
			//	{
			//			err_log << pw_out[ll] << endl;  
						//cout<<"pw_out[ll]: Boundary Element: "<< pw_out[ll] <<endl;    
			//	}
			//}
			//err_log << adj_error.size() << endl;
			double normalized_error = log(error[i]);
			err_log << normalized_error << endl;

			
		}
		else
		{
			// For the internal elements
			// cout<< "Internal Elements Activate" << endl;
			
			int ndof_w = ml_elidata.ndof_w;
			// KUNAL
			// cout<< "ndof_w: "<< ndof_w<<endl;
			// Taking the quadrature nodes from the first face
			int fcnr = ml_elidata.faces[0];
			// int fcnr = ml_elidata.faces[0];
			ML_FacetData<D, COMP> & ml_fd_1 = *ml_fadata[fcnr];
			int nip_1 = ml_fd_1.nip;
			
			// cout << "nip_1: "<< nip_1<<endl;
			
			// Graph Neural Network
			// err_log << i <<","<< ml_elidata.nip << "," << ml_elidata.ndof_w*(D+1) + nip_1*(D+1)*ml_elidata.nf + (D*(D+1))/2 + Model::NumParam << "," << ml_elidata.nf_bc << endl;
			
			// Feedforward neural network
			//err_log << "starting line 2 (internal elements)" << endl;
			err_log << "1" << "," << (ml_elidata.ndof_w*(D+1) + nip_1*(D+1)*ml_elidata.nf)*COMP + (D*(D+1))/2 + Model::NumParam << "," << ml_elidata.nf_bc << endl;

			//err_log << "volume, beta, theta (internal elements)" << endl;
			err_log << ml_elidata.vol << endl << ml_elidata.beta << endl << ml_elidata.theta << endl;
			
			for(int pr = 0; pr<Model::NumParam;pr++)
			{
				err_log << param(pr) <<endl;
				//cout << "param(pr)" << param(pr) <<endl;
	    		}
	    		
	    		// Generating the input data for solution and gradient
			//err_log << "solution (internal elements)" << endl;
			for (int ll = 0; ll < COMP; ++ll)
			{
				for (int pp = 0; pp < ndof_w; ++pp)
				{
					err_log << ml_elidata.vecW[pp+ll*ndof_w] << endl;
					//cout << "ml_elidata.vecW[pp+ll*ndof_w]: " << ml_elidata.vecW[pp+ll*ndof_w] << endl;
				}
			}

			// KUNAL
			// cout<< "Solution Printed" << endl;

			// Generating the input data for solution and gradient
			if (Model::Diffusion || Model::Source){
			//err_log << "solution gradient (internal elements)" << endl;
			for (int ll = 0; ll < COMP; ++ll)
			{
				for (int pp = 0; pp < ndof_w; ++pp)
				{
					for(int dd = 0; dd < D; ++dd)
					{
						err_log << ml_elidata.vecQ[pp+ndof_w*(dd+D*ll)] << endl;
						//cout << "ml_elidata.vecQ[pp+ndof_w*(dd+D*ll)]: "<<ml_elidata.vecQ[pp+ndof_w*(dd+D*ll)] << endl;
					}
				}
			}
			}

			// KUNAL
			// cout<< "Solution gradient printed" << endl;

			// KUNAL
			// cout<<"ml_elidata.nf: "<<ml_elidata.nf<<endl;
			
			
			// KUNAL: Jump data starts from here
			// cout<<"Jump Data: "<< endl; 
			for(int j = 0 ; j < ml_elidata.nf ; j ++)
			{

				int fcnr = ml_elidata.faces[j];
				ML_FacetData<D, COMP> & ml_fd = *ml_fadata[fcnr];
				int nip = ml_fd.nip;

				if (!fadata[fcnr]) continue;
				FacetData<D, COMP> & fd = *fadata[fcnr];
				// KUNAL: Why they are doing this?
				// They have used a condition to only pass the internal elements
				// But without this the code is crashing
				if (fd.ndof_l == 0) continue; // skip boundary faces
				if (fd.bcnr == -2) continue; // zero-measure face

				// Jump in w on quarature points
				// cout<< "Jump Data: w: Printing" <<endl;
				//err_log << "solution jump (internal elements)" << endl;
				for (int qp = 0; qp < nip; ++qp) 
				{
					double *jump_w = &ml_fd.w_jump[qp*COMP];
					for(int ll=0;ll<COMP;ll++)
					{
						err_log << jump_w[ll] <<endl;
						//cout<< "jump_w[ll]: "<< jump_w[ll] <<endl;
					}
				}
				
				// Jump in q on quadrature points
				// cout<< "Jump Data: q: Printing" <<endl;
				if (Model::Diffusion || Model::Source){
				//err_log << "solution gradient jump (internal elements)" << endl;
				for (int qp = 0; qp < nip; ++qp) 
				{
					double *jump_q = &ml_fd.q_jump[qp * COMP * D];

					for(int ll=0;ll<COMP;ll++)
					{

						for(int dd=0;dd<D;dd++)
						{
							err_log << jump_q[dd+D*ll] << endl;
							//cout<< "jump_q[dd+D*ll]: "<< jump_q[dd+D*ll] <<endl;
						}
					}

				}
				}

			}
			
			// KUNAL
			// Internal Element ML output
			//err_log << "ML output data (internal element)" << endl;
			// ML_OutputData<D, COMP> & ml_outdata = *ml_outputdata[i];
			// for (int j = 0; j < ml_outdata.nip; ++j) 
			// {
			// 	const double * pw_out  = &ml_outdata.w[j * COMP];
			// 	for (int ll = 0; ll < COMP; ++ll)
			// 	{
			// 			err_log << pw_out[ll] << endl;  
			// 			//cout<<"pw_out[ll]: Internal Element: "<< pw_out[ll] <<endl;    
			// 	}
			// }
			double normalized_error = log(error[i]);
			err_log << normalized_error << endl;
		}

	}//end of loop over elements
	// err_log.close();
}
