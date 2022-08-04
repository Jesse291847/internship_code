data {
  int<lower=0> N; // Sample size
  int<lower=1> dim; // Number of dimensions
  vector[dim] x[N]; // Value for each sample on each dimension  
  cov_matrix[dim] S; // Psi matrix (parameter in inv.wishart)
  } 
  
parameters {
  vector[dim] mu; // mean vector
  cov_matrix[dim] sigma; // ***covariance matrix (what we are after)***
  }

model {
  for (n in 1:N) x[n] ~ multi_normal(mu, sigma);
  sigma ~ inv_wishart(dim, S);
  mu ~ normal(0, 1);
  }
  
generated quantities {
	matrix[dim, dim] prec;
	matrix[dim, dim] pcor;
	real cent1;
	real cent2;
	real cent3;
	real cent4;
	real cent5;
	real cent6;
	real cent7;
	real cent8;
	real cent9;
	real c1_1;
    	real c1_2;
    	real c1_3;
    	real c1_4;
    	real c1_5;
    	real c1_6;
    	real c1_7;
    	real c1_8;
    	real c1_9;
    	real c2_1;
    	real c2_2;
    	real c2_3;
   	real c2_4;
   	real c2_5;
    	real c2_6;
    	real c2_7;
    	real c2_8;
    	real c2_9;
    	real c3_1;
   	real c3_2;
   	real c3_3;
   	real c3_4;
   	real c3_5;
   	real c3_6;
   	real c3_7;
   	real c3_8;
   	real c3_9;
    	real c4_1;
    	real c4_2;
   	real c4_3;
    	real c4_4;
    	real c4_5;
   	real c4_6;
   	real c4_7;
   	real c4_8;
   	real c4_9;
	real c5_1;
    	real c5_2;
    	real c5_3;
    	real c5_4;
    	real c5_5;
    	real c5_6;
    	real c5_7;
    	real c5_8;
    	real c5_9;
    	real c6_1;
    	real c6_2;
    	real c6_3;
    	real c6_4;
    	real c6_5;
    	real c6_6;
    	real c6_7;
    	real c6_8;
    	real c6_9;
    	real c7_1;
    	real c7_2;
    	real c7_3;
    	real c7_4;
    	real c7_5;
    	real c7_6;
    	real c7_7;
    	real c7_8;
    	real c7_9;
    	real c8_1;
    	real c8_2;
    	real c8_3;
    	real c8_4;
    	real c8_5;
    	real c8_6;
    	real c8_7;
    	real c8_8;
    	real c8_9;
    	real c9_1;
    	real c9_2;
    	real c9_3;
    	real c9_4;
    	real c9_5;
    	real c9_6;
    	real c9_7;
    	real c9_8;
    	real c9_9;
    	
	
	for (n in 1:N) 
		// Precision matrix
		prec = inverse(sigma);
		
		// Partial correlation matrix (GGM)
		for(i in 1:dim) {
		  for(j in 1:dim) {
		    pcor[i,j] = prec[i,j] * (-1) / (sqrt(prec[i,i]) * sqrt(prec[j,j]));
		  }
		}
		
		// Centrality (expected influence)
		cent1 = sum(pcor[1,]) - pcor[1,1];
    		cent2 = sum(pcor[2,]) - pcor[2,2];
   		cent3 = sum(pcor[3,]) - pcor[3,3];
    		cent4 = sum(pcor[4,]) - pcor[4,4];
    		cent5 = sum(pcor[5,]) - pcor[5,5];
    		cent6 = sum(pcor[6,]) - pcor[6,6];
    		cent7 = sum(pcor[7,]) - pcor[7,7];
    		cent8 = sum(pcor[8,]) - pcor[8,7];
    		cent9 = sum(pcor[9,]) - pcor[9,7];
    		
    		// Centrality difference distributions
    		c1_1 = cent1 - cent1;
    		c1_2 = cent1 - cent2;
    		c1_3 = cent1 - cent3;
    		c1_4 = cent1 - cent4;
    		c1_5 = cent1 - cent5;
    		c1_6 = cent1 - cent6;
    		c1_7 = cent1 - cent7;
    		c1_8 = cent1 - cent8;
    		c1_9 = cent1 - cent9;
    		c2_1 = cent2 - cent1;
    		c2_2 = cent2 - cent2;
    		c2_3 = cent2 - cent3;
   	 	c2_4 = cent2 - cent4;
   	 	c2_5 = cent2 - cent5;
    		c2_6 = cent2 - cent6;
    		c2_7 = cent2 - cent7;
    		c2_8 = cent2 - cent8;
    		c2_9 = cent2 - cent9;
    		c3_1 = cent3 - cent1;
   	 	c3_2 = cent3 - cent2;
   	 	c3_3 = cent3 - cent3;
   	 	c3_4 = cent3 - cent4;
   	 	c3_5 = cent3 - cent5;
   	 	c3_6 = cent3 - cent6;
   	 	c3_7 = cent3 - cent7;
   	 	c3_8 = cent3 - cent8;
   	 	c3_9 = cent3 - cent9;
    		c4_1 = cent4 - cent1;
    		c4_2 = cent4 - cent2;
   	 	c4_3 = cent4 - cent3;
    		c4_4 = cent4 - cent4;
    		c4_5 = cent4 - cent5;
   	 	c4_6 = cent4 - cent6;
   	 	c4_7 = cent4 - cent7;
   	 	c4_8 = cent4 - cent8;
   	 	c4_9 = cent4 - cent9;
		c5_1 = cent5 - cent1;
    		c5_2 = cent5 - cent2;
    		c5_3 = cent5 - cent3;
    		c5_4 = cent5 - cent4;
    		c5_5 = cent5 - cent5;
    		c5_6 = cent5 - cent6;
    		c5_7 = cent5 - cent7;
    		c5_8 = cent5 - cent8;
    		c5_9 = cent5 - cent9;
    		c6_1 = cent6 - cent1;
    		c6_2 = cent6 - cent2;
    		c6_3 = cent6 - cent3;
    		c6_4 = cent6 - cent4;
    		c6_5 = cent6 - cent5;
    		c6_6 = cent6 - cent6;
    		c6_7 = cent6 - cent7;
    		c6_8 = cent6 - cent8;
    		c6_9 = cent6 - cent9;
    		c7_1 = cent7 - cent1;
    		c7_2 = cent7 - cent2;
    		c7_3 = cent7 - cent3;
    		c7_4 = cent7 - cent4;
    		c7_5 = cent7 - cent5;
    		c7_6 = cent7 - cent6;
    		c7_7 = cent7 - cent7;
    		c7_8 = cent7 - cent8;
    		c7_9 = cent7 - cent9;
    		c8_1 = cent8 - cent1;
    		c8_2 = cent8 - cent2;
    		c8_3 = cent8 - cent3;
    		c8_4 = cent8 - cent4;
    		c8_5 = cent8 - cent5;
    		c8_6 = cent8 - cent6;
    		c8_7 = cent8 - cent7;
    		c8_8 = cent8 - cent8;
    		c8_9 = cent8 - cent9;
    		c9_1 = cent9 - cent1;
    		c9_2 = cent9 - cent2;
    		c9_3 = cent9 - cent3;
    		c9_4 = cent9 - cent4;
    		c9_5 = cent9 - cent5;
    		c9_6 = cent9 - cent6;
    		c9_7 = cent9 - cent7;
    		c9_8 = cent9 - cent8;
    		c9_9 = cent9 - cent9;
  }
