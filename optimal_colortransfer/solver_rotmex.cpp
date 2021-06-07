#include <Eigen/Dense>
#include "mex.h"


using namespace Eigen;


MatrixXcd rotsolver(double *d)
{
	// Compute coefficients
   // const double* d = data.data();
    
    //for (int i = 0;i<20;i++)
    //    mexPrintf("%5.2f\n",d[i]);
    
    VectorXd coeffs(30);
    coeffs[0] = d[0]*d[9] - d[1]*d[10] - d[2]*d[11] + d[18];
    coeffs[1] = 2*d[1]*d[9] + 2*d[0]*d[10];
    coeffs[2] = -d[0]*d[9] + d[1]*d[10] - d[2]*d[11] + d[18];
    coeffs[3] = 2*d[2]*d[9] + 2*d[0]*d[11];
    coeffs[4] = 2*d[2]*d[10] + 2*d[1]*d[11];
    coeffs[5] = -d[0]*d[9] - d[1]*d[10] + d[2]*d[11] + d[18];
    coeffs[6] = 2*d[2]*d[10] - 2*d[1]*d[11];
    coeffs[7] = -2*d[2]*d[9] + 2*d[0]*d[11];
    coeffs[8] = 2*d[1]*d[9] - 2*d[0]*d[10];
    coeffs[9] = d[0]*d[9] + d[1]*d[10] + d[2]*d[11] + d[18];
    coeffs[10] = d[3]*d[12] - d[4]*d[13] - d[5]*d[14] + d[19];
    coeffs[11] = 2*d[4]*d[12] + 2*d[3]*d[13];
    coeffs[12] = -d[3]*d[12] + d[4]*d[13] - d[5]*d[14] + d[19];
    coeffs[13] = 2*d[5]*d[12] + 2*d[3]*d[14];
    coeffs[14] = 2*d[5]*d[13] + 2*d[4]*d[14];
    coeffs[15] = -d[3]*d[12] - d[4]*d[13] + d[5]*d[14] + d[19];
    coeffs[16] = 2*d[5]*d[13] - 2*d[4]*d[14];
    coeffs[17] = -2*d[5]*d[12] + 2*d[3]*d[14];
    coeffs[18] = 2*d[4]*d[12] - 2*d[3]*d[13];
    coeffs[19] = d[3]*d[12] + d[4]*d[13] + d[5]*d[14] + d[19];
    coeffs[20] = d[6]*d[15] - d[7]*d[16] - d[8]*d[17] + d[20];
    coeffs[21] = 2*d[7]*d[15] + 2*d[6]*d[16];
    coeffs[22] = -d[6]*d[15] + d[7]*d[16] - d[8]*d[17] + d[20];
    coeffs[23] = 2*d[8]*d[15] + 2*d[6]*d[17];
    coeffs[24] = 2*d[8]*d[16] + 2*d[7]*d[17];
    coeffs[25] = -d[6]*d[15] - d[7]*d[16] + d[8]*d[17] + d[20];
    coeffs[26] = 2*d[8]*d[16] - 2*d[7]*d[17];
    coeffs[27] = -2*d[8]*d[15] + 2*d[6]*d[17];
    coeffs[28] = 2*d[7]*d[15] - 2*d[6]*d[16];
    coeffs[29] = d[6]*d[15] + d[7]*d[16] + d[8]*d[17] + d[20];



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

void myquat2rot(double q2,double q3, double q4, double *Re){
    double q1 = 1.0;
    Re[0] = q1*q1+q2*q2-q3*q3-q4*q4;
    Re[1] = 2*q1*q4+2*q2*q3;
    Re[2] = 2*q2*q4-2*q1*q3;
    Re[3] = 2*q2*q3-2*q1*q4;
    Re[4] = q1*q1-q2*q2+q3*q3-q4*q4;
    Re[5] = 2*q1*q2+2*q3*q4;
    Re[6] = 2*q1*q3+2*q2*q4;
    Re[7] = 2*q3*q4-2*q1*q2;
    Re[8] = q1*q1-q2*q2-q3*q3+q4*q4;   
    }

double myprod(double *yx,double *R){
    double p = 0;
    for (int i = 0;i<9;i++)
        p+=yx[i]*R[i];
    return p;
    }


void optrot(int nn, double *xx, double *yy, double *ss,double *yx, double *R, mxLogical *inliers)
{
    double input[21];
    double Re[9];
    MatrixXcd sols;
    //const VectorXd data = Map<const VectorXd>(input, 21);
    const double lille = 1E-9;              
    double bq2,bq3,bq4,q2,q3,q4,qn2;
    int nrins,bestins = 0;
    
    for (int i1 = 0;i1 < nn;i1++){
        input[0] = -2*yy[3*i1];
        input[1] = -2*yy[3*i1+1];
        input[2] = -2*yy[3*i1+2];
        input[9] = xx[3*i1];
        input[10] = xx[3*i1+1];
        input[11] = xx[3*i1+2];
        input[18] = ss[i1];     
        for (int i2 = i1+1;i2 < nn;i2++){
                input[3] = -2*yy[3*i2];
                input[4] = -2*yy[3*i2+1];
                input[5] = -2*yy[3*i2+2];
                input[12] = xx[3*i2];
                input[13] = xx[3*i2+1];
                input[14] = xx[3*i2+2];
                input[19] = ss[i2];                 
                for (int i3 = i2+1; i3 < nn;i3++){
                    input[6] = -2*yy[3*i3];
                    input[7] = -2*yy[3*i3+1];
                    input[8] = -2*yy[3*i3+2];
                    input[15] = xx[3*i3];
                    input[16] = xx[3*i3+1];
                    input[17] = xx[3*i3+2];
                    input[20] = ss[i3];                 
                    sols = rotsolver(input);
                    //int count = 0;
                    for (Index i = 0; i < sols.size(); i = i+3) {
                        if ((abs(sols(i).imag())+abs(sols(i+1).imag())+abs(sols(i+2).imag()))<lille){
                           // count++;
                            q2 = sols(i).real();
                            q3 = sols(i+1).real();
                            q4 = sols(i+2).real();
                            qn2 = 1+q2*q2+q3*q3+q4*q4;
                            //mexPrintf("%5.4f %5.4f %5.4f \n",q2,q3,q4);
                            nrins = 0;
                            for (int j = 0; j < nn;j++){
                                myquat2rot(q2,q3,q4,Re);
                                nrins += ((myprod(&yx[9*j],Re)/qn2+ss[j])<lille);
                                //mexPrintf("%15.12f\n",(myprod(&yx[9*j],Re)/qn2+ss[j]));
                                }
                            if (nrins>bestins){
                                bestins = nrins;
                                bq2 = q2;
                                bq3 = q3;
                                bq4 = q4;
                               // mexPrintf("%d %d %d %d %d\n",i1+1,i2+1,i3+1,count,bestins);
                                }
                            }
                        }
                
		
                
    }
                }
        
            }
    //mexPrintf("%d\n",bestins);
    myquat2rot(bq2,bq3,bq4,R);
    qn2 = 1+bq2*bq2+bq3*bq3+bq4*bq4;
    for (int j = 0; j < nn;j++){
        inliers[j] = ((myprod(&yx[9*j],R)/qn2+ss[j])<0);
        }
    for (int i = 0;i<9;i++)
        R[i] = R[i]/qn2;
    }
    


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 4) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:test:nrhs", "Four input required.");
	}
	if (nlhs != 2) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:test:nlhs", "Two outputs required.");
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
    if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3])) {
		mexErrMsgIdAndTxt("automatic_generator_cvpr:test:notDouble", "Input data must be type double.");
	}
    
	double *xx = mxGetPr(prhs[0]);
    double *yy = mxGetPr(prhs[1]);
    double *ss = mxGetPr(prhs[2]);
    double *yx = mxGetPr(prhs[3]);
    int nn = mxGetNumberOfElements(prhs[2]);
    plhs[0] = mxCreateDoubleMatrix(3,3,mxREAL);
    plhs[1] = mxCreateLogicalMatrix(nn,1);
    double *R = mxGetPr(plhs[0]);
    mxLogical *inliers = mxGetLogicals(plhs[1]);
    //mexPrintf("hepp %d \n",nn);
    
    optrot(nn,xx,yy,ss,yx,R,inliers);

}


