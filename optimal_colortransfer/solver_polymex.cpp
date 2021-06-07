#include <Eigen/Dense>
#include "mex.h"




using namespace Eigen;




MatrixXcd polysolver(double *d)
{
	// Compute coefficients
    //const double* d = data.data();
    VectorXd coeffs(30);
    coeffs[0] = std::pow(d[0],6) - 9*std::pow(d[0],4)*d[6];
    coeffs[1] = 2*std::pow(d[0],5) - 12*std::pow(d[0],3)*d[6];
    coeffs[2] = std::pow(d[0],4) - 4*std::pow(d[0],2)*d[6];
    coeffs[3] = 2*std::pow(d[0],4) - 6*std::pow(d[0],2)*d[6];
    coeffs[4] = 2*std::pow(d[0],3) - 4*d[0]*d[6];
    coeffs[5] = std::pow(d[0],2) - d[6];
    coeffs[6] = -2*std::pow(d[0],3)*d[1];
    coeffs[7] = -2*std::pow(d[0],2)*d[1];
    coeffs[8] = -2*d[0]*d[1];
    coeffs[9] = std::pow(d[1],2) - d[6];
    coeffs[10] = std::pow(d[2],6) - 9*std::pow(d[2],4)*d[6];
    coeffs[11] = 2*std::pow(d[2],5) - 12*std::pow(d[2],3)*d[6];
    coeffs[12] = std::pow(d[2],4) - 4*std::pow(d[2],2)*d[6];
    coeffs[13] = 2*std::pow(d[2],4) - 6*std::pow(d[2],2)*d[6];
    coeffs[14] = 2*std::pow(d[2],3) - 4*d[2]*d[6];
    coeffs[15] = std::pow(d[2],2) - d[6];
    coeffs[16] = -2*std::pow(d[2],3)*d[3];
    coeffs[17] = -2*std::pow(d[2],2)*d[3];
    coeffs[18] = -2*d[2]*d[3];
    coeffs[19] = std::pow(d[3],2) - d[6];
    coeffs[20] = std::pow(d[4],6) - 9*std::pow(d[4],4)*d[6];
    coeffs[21] = 2*std::pow(d[4],5) - 12*std::pow(d[4],3)*d[6];
    coeffs[22] = std::pow(d[4],4) - 4*std::pow(d[4],2)*d[6];
    coeffs[23] = 2*std::pow(d[4],4) - 6*std::pow(d[4],2)*d[6];
    coeffs[24] = 2*std::pow(d[4],3) - 4*d[4]*d[6];
    coeffs[25] = std::pow(d[4],2) - d[6];
    coeffs[26] = -2*std::pow(d[4],3)*d[5];
    coeffs[27] = -2*std::pow(d[4],2)*d[5];
    coeffs[28] = -2*d[4]*d[5];
    coeffs[29] = std::pow(d[5],2) - d[6];



	// Setup elimination template
	static const int coeffs0_ind[] = { 0,10,1,0,10,11,20,2,1,11,12,21,2,12,22,0,10,20,3,13,1,0,10,11,20,21,4,3,13,14,2,1,11,12,21,23,22,4,14,2,12,22,24,3,13,10,0,20,23,5,15,4,3,13,14,11,1,21,23,24,5,15,4,14,12,2,22,24,25,5,15,13,3,23,25,5,15,14,4,24,25,0,10,20,6,16,1,11,0,20,10,21,7,6,16,17,26,2,12,1,21,11,22,7,17,27,2,22,12,6,16,3,13,20,10,0,23,26,8,18,7,6,16,17,26,4,14,3,23,21,11,1,13,24,27,8,18,7,17,27,28,4,24,22,12,2,14,6,16,20,10,0,26,9,19,7,17,6,26,21,11,1,16,27,9,19,29,7,27,22,12,2,17,8,18,16,6,26,5,15,23,13,3,25,28,8,18,17,7,27,28,5,25,24,14,4,15,15,5,25 };
	static const int coeffs1_ind[] = { 29,19,9,9,19,26,16,6,29,9,19,8,18,26,16,6,23,13,3,28,29,9,29,27,17,7,19,9,19,29,8,28,27,17,7,24,14,4,18,29,19,9,28,18,8,19,9,29,28,18,8,25,15,5,18,8,28,25,15,5 };
	static const int C0_ind[] = { 0,3,26,27,28,29,38,52,53,54,55,64,79,80,90,108,111,129,130,133,134,135,136,137,141,155,156,157,158,159,160,161,162,163,167,168,181,183,184,187,188,193,194,212,215,216,217,218,233,234,237,238,239,240,241,242,243,244,245,259,261,262,265,266,268,269,270,271,272,290,293,294,295,296,311,317,318,320,321,322,323,351,352,362,364,367,377,378,379,380,387,388,390,391,392,393,402,403,404,405,406,413,414,417,418,428,431,432,439,446,449,455,456,459,460,461,466,467,468,471,472,473,474,475,479,481,482,483,484,485,486,487,491,492,493,495,496,499,500,505,506,509,510,511,512,513,517,533,534,540,541,542,544,546,549,559,560,561,562,566,567,568,569,570,573,574,584,587,588,592,593,594,595,602,605,606,607,608,611,612,615,616,617,622,623,629,630,632,633,634,635,639,640,641,642,643,647,658,659,660 } ;
	static const int C1_ind[] = { 20,21,22,39,40,46,47,48,50,56,59,65,66,69,70,71,72,73,74,76,77,93,94,98,99,100,101,109,110,115,119,120,121,122,123,124,125,126,127,147,148,149,150,151,152,164,165,166,173,174,175,176,177,178,190,191,192,199,200,201 };

	Matrix<double,26,26> C0; C0.setZero();
	Matrix<double,26,8> C1; C1.setZero();
	for (int i = 0; i < 200; i++) { C0(C0_ind[i]) = coeffs(coeffs0_ind[i]); }
	for (int i = 0; i < 60; i++) { C1(C1_ind[i]) = coeffs(coeffs1_ind[i]); } 

	Matrix<double,26,8> C12 = C0.partialPivLu().solve(C1);




	// Setup action matrix
	Matrix<double,11, 8> RR;
	RR << -C12.bottomRows(3), Matrix<double,8,8>::Identity(8, 8);

	static const int AM_ind[] = { 8,5,0,7,1,9,10,2 };
	Matrix<double, 8, 8> AM;
	for (int i = 0; i < 8; i++) {
		AM.row(i) = RR.row(AM_ind[i]);
	}

	Matrix<std::complex<double>, 3, 8> sols;
	sols.setZero();

	// Solve eigenvalue problem
	EigenSolver<Matrix<double, 8, 8> > es(AM);
	ArrayXcd D = es.eigenvalues();	
	ArrayXXcd V = es.eigenvectors();

V = (V / V.row(0).array().replicate(8, 1)).eval();


    sols.row(0) = V.row(1).array();
    sols.row(1) = V.row(3).array();
    sols.row(2) = D.transpose().array();



	return sols;
}
// Action =  z
// Quotient ring basis (V) = 1,x,x*z,y,y*z,z,z^2,z^3,
// Available monomials (RR*V) = x*z^2,y*z^2,z^4,1,x,x*z,y,y*z,z,z^2,z^3,

double distapprox_poly_full(double a,double b,double c,double x,double y){
    
    double k = 3*a*x*x+2*b*x+c;
    double m = (a*x*x*x+b*x*x+c*x)-k*x;

    return (k*x-y+m)*(k*x-y+m)/(k*k+1);
    }


void optpoly(int nn, double *xx, double *yy, double ep, double *poly, mxLogical *inliers)
{
    double input[7];
    //double polye[9];
    MatrixXcd sols;
    //const VectorXd data = Map<const VectorXd>(input, 21);
    const double lille = 1E-9;              
    //double bq2,bq3,bq4,q2,q3,q4,qn2;
    int nrins,bestins = 0;
    double s1,s2,d,e;
    
    input[6] = ep;
    for (int i1 = 0;i1 < nn;i1++){
        input[0] = xx[i1];
        input[1] = yy[i1];
        for (int i2 = i1+1;i2 < nn;i2++){
                input[2] = xx[i2];
                input[3] = yy[i2];
                for (int i3 = i2+1; i3 < nn;i3++){
                    input[4] = xx[i3];
                    input[5] = yy[i3];
                    sols = polysolver(input);
                    //int count = 0;
                    for (Index i = 0; i < sols.size(); i = i+3) {
                        if ((abs(sols(i).imag())+abs(sols(i+1).imag())+abs(sols(i+2).imag()))<lille){                   
                            d = 2*sols(i+1).real()/(3*sols(i).real());
                            e = sols(i+2).real()/(3*sols(i).real());
                            if ((d*d)>=(4*e)){
                                s1 = -d/2-sqrt(d*d/4-e);
                                s2 = -d/2+sqrt(d*d/4-e);
                                
                                if  ((s1<0 || s1>1) && (s2<0 || s2>1)) {
                                
                                    nrins = 0;
                                    for (int j = 0; j < nn;j++){
                                    nrins+=(distapprox_poly_full(sols(i).real(),sols(i+1).real(),sols(i+2).real(),xx[j],yy[j])<(ep+lille));
                                    }
                                    if (nrins>bestins){
                                        bestins = nrins;
                                        poly[0] = sols(i).real();
                                        poly[1] = sols(i+1).real();
                                        poly[2] = sols(i+2).real();
                                        
                                    }
                                    if (bestins==nn){
                                        i1 = nn;
                                        i2 = nn;
                                        i3 = nn;
                                        }
                                }
                            
                            }
                        }
                    }
                    
                }
            }
       }
    
    
    
    
    
    for (int j = 0; j < nn;j++)
        inliers[j] = (distapprox_poly_full(poly[0],poly[1],poly[2],xx[j],yy[j])<(ep+lille));
        
 
    }
    



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 3) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:testpoly:nrhs", "Three input required.");
	}
	if (nlhs != 2) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:testpoly:nlhs", "Two outputs required.");
	}    
	if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:testpoly:notDouble", "Input data must be type double.");
	}
		if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:test:notDouble", "Input data must be type double.");
	}
    if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:test:notDouble", "Input data must be type double.");
	}

       if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:test:notDouble", "Input data must be type double.");
	}

    double *xx = mxGetPr(prhs[0]);
    double *yy = mxGetPr(prhs[1]);
    double ep = mxGetScalar(prhs[2]);
    int nn = mxGetNumberOfElements(prhs[0]);
    plhs[0] = mxCreateDoubleMatrix(3,1,mxREAL);
    plhs[1] = mxCreateLogicalMatrix(nn,1);
    double *poly = mxGetPr(plhs[0]);
    mxLogical *inliers = mxGetLogicals(plhs[1]);
    //mexPrintf("hepp %d \n",nn);
    
    optpoly(nn,xx,yy,ep,poly,inliers);

    
	
	
}


