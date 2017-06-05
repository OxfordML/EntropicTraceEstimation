#include <iostream>
#include <cmath>
#include "mex.h"

using namespace std;

mxArray * blockHankel(const mxArray *mxY, const int dimsY[], int nBlockRows)
{
    double *Y = mxGetPr(mxY);
    int nc = dimsY[1] - nBlockRows + 1;
    int nr = nBlockRows * dimsY[0];
    mxArray *mxHy = mxCreateDoubleMatrix(nr, nc, mxREAL);
    double *Hy = mxGetPr(mxHy);
    for(size_t k = 0; k < nc; ++k)
        for(size_t i = 0; i < nBlockRows; ++i)
            for(size_t j = 0; j < dimsY[0]; ++j)
            {
                Hy[k*nBlockRows*dimsY[0]+i*dimsY[0]+j] = Y[(i+k)*dimsY[0]+j];
            }
    return mxHy;
}

void normalizeGram(double *Gy, size_t m)
{
    size_t nElements = m * m;
    double sum = 0.0;
    for(size_t i = 0; i < nElements; ++i)
    {
        sum += Gy[i] * Gy[i];
    }
    sum = sqrt(sum);
    for(size_t i = 0; i < nElements; ++i)
    {
        Gy[i] /= sum;
    }
}

void regularizeGram(double *Gy, size_t m, double sigma)
{
    for(size_t i = 0; i < m; ++i)
        for(size_t j = 0; j < m; ++j)
        {
            if(i == j)
                Gy[i*m+j] += sigma;
        }
}

mxArray * gram(const mxArray *mxHy, bool normalize, double sigma)
{
    const int *dimsHy = mxGetDimensions(mxHy);
    int m = dimsHy[0];
    int n = dimsHy[1];
    double *Hy = mxGetPr(mxHy);
    mxArray *mxGy = mxCreateDoubleMatrix(m, m, mxREAL);
    double *Gy = mxGetPr(mxGy);
    double sum;
    
    for(size_t j = 0; j < m; ++j)
        for(size_t i = 0; i < m; ++i)
        {
            if(i < j)
            {
                Gy[j*m+i] = Gy[i*m+j];
                continue;
            }
            sum = 0.0;
            for(size_t k = 0; k < n; ++k)
            {
                sum += Hy[k*m+i] * Hy[k*m+j];
            }
            Gy[j*m+i] = sum;
        }
    
    if(normalize)
    {
        normalizeGram(Gy, m);
        regularizeGram(Gy, m, sigma);
    }
    return mxGy;
}

mxArray * getGram(const mxArray *mxY, const mxArray *mxNr)
{
    if (mxGetNumberOfDimensions(mxY) != 2)
        mexErrMsgTxt("Input dimension must be 2.\n");
    if (!mxIsScalar(mxNr))
        mexErrMsgTxt("Second input must be a scalar.\n");
    
    const int *dimsY = mxGetDimensions(mxY);
    if (dimsY[0] == 0 || dimsY[1] == 0)
        mexErrMsgTxt("Input should not be empty");
    
    int nBlockRows = (int) mxGetScalar(mxNr);
    if (nBlockRows > dimsY[1])
        mexErrMsgTxt("number of rows should not be larger than sequence length");
    
    mxArray *mxHy = blockHankel(mxY, dimsY, nBlockRows);
    
    bool normalize = true;
    double sigma = 1e-4;
    mxArray *mxGy = gram(mxHy, normalize, sigma);
    
    mxDestroyArray(mxHy);
    
    return mxGy;
}

mxArray * getGram2(const mxArray *mxY, const mxArray *mxNr)
{
    if (mxGetNumberOfDimensions(mxY) != 2)
        mexErrMsgTxt("Input dimension must be 2.\n");
    if (!mxIsScalar(mxNr))
        mexErrMsgTxt("Second input must be a scalar.\n");
    
    const int *dimsY = mxGetDimensions(mxY);
    if (dimsY[0] == 0 || dimsY[1] == 0)
        mexErrMsgTxt("Input should not be empty");
    
    int nBlockRows = (int) mxGetScalar(mxNr);
    if (nBlockRows > dimsY[1])
        mexErrMsgTxt("number of rows should not be larger than sequence length");
    
    double *Y = mxGetPr(mxY);
    int nc = dimsY[1] - nBlockRows + 1;
    int nr = nBlockRows * dimsY[0];
    mxArray *mxGy = mxCreateDoubleMatrix(nr, nr, mxREAL);
    double *Gy = mxGetPr(mxGy);
    
    double sum;
    int nrY = dimsY[0];
    for(size_t j = 0; j < nr; ++j)
        for(size_t i = 0; i < nr; ++i)
        {
            if(i < j)
            {
                Gy[i+j*nr] = Gy[j+i*nr];
                continue;
            }
            if(j > nrY)
            {
                Gy[i+j*nr] = Gy[i-nrY+(j-nrY)*nr] - Y[i-nrY] * Y[j-nrY] + Y[i+(nc-1)*nrY] * Y[j+(nc-1)*nrY];
                continue;
            }
            sum = 0.0;
            for(size_t k = 0; k < nc; ++k)
            {
                sum += Y[i+k*nrY] * Y[j+k*nrY];
            }
            Gy[i+j*nr] = sum;
        }
    
//    bool normalize = true;
    double sigma = 1e-4;
    normalizeGram(Gy, nr);
    regularizeGram(Gy, nr, sigma);
    
    return mxGy;
}

mxArray * process(const mxArray *cellArrayPtr, const mxArray *mxNr)
{
    if(!mxIsCell(cellArrayPtr))
    {
        mexPrintf("the first input must be a cell");
    }
    const mxArray *cellPtr;
    mwSize nCells = mxGetNumberOfElements(cellArrayPtr);
    mwSize ndim = 2;
    mwSize dims[2];
    dims[0] = 1;
    dims[1] = nCells;
    mxArray *outputCell = mxCreateCellArray(ndim, dims);
    mxArray *Gy;
    for(size_t i = 0; i < nCells; ++i)
    {
        cellPtr = mxGetCell(cellArrayPtr, i);
        if(cellPtr == NULL)
        {
            mexPrintf("\tEmpty Cell\n");
        }
        Gy = getGram2(cellPtr, mxNr);
        mxSetCell(outputCell, i, Gy);
    }
    return outputCell;
}

// matlab entry point
// Hy = mex_getHH(y, opt)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs != 2)
        mexErrMsgTxt("Wrong number of inputs\n");
    if (nlhs != 1)
        mexErrMsgTxt("Wrong number of outputs\n");
    plhs[0] = process(prhs[0],prhs[1]);
}