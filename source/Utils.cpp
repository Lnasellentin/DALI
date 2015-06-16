/*
 Some ultilities for Dali.
 Ln(a) Sellentin
 Universität Heidelberg
 */

#include "Utils.h"

void printxyz(double x, double y, double z)
{
    printf("%.3f %.3f %.3f \n",x,y,z);
}



double vec1Matvec2(gsl_vector* vec1, gsl_matrix* mat, gsl_vector* vec2, int dim)
{
    gsl_vector * hilf = gsl_vector_alloc(dim);
    double scalarResult;
    gsl_blas_dgemv(CblasNoTrans, 1.0, mat, vec2, 0.0, hilf);  //Mat*vec2
    gsl_blas_ddot(hilf,vec1 , &scalarResult);

    gsl_vector_free(hilf);
    return scalarResult;
}



void printmatrix(gsl_matrix* mat, int zeilen, int spalten)
{
    for(int i = 0; i < zeilen; i ++)
    {
        for(int j = 0; j < spalten; j++)
        {
            printf("%.2e ",gsl_matrix_get(mat, i, j));
        }
        printf("\n");
    }
}


void invmatrix(gsl_matrix* tobeinverted, unsigned int dim, gsl_matrix* inverted)
{
    /*Takes the matrix "tobeinverted", which is of dimension "dim", and inverts it. The result is stored in an externally provided gsl_matrix "inverted". */

    /*Do not allocate space for the resulting matrix within this function!
     The code would compile, and store the matrix here - but you want to use it outside this function. However, the outside allocated matrix will then not be the one you just calculated!*/
    gsl_matrix* hilf = gsl_matrix_alloc(dim,dim);

    /*use memcopy, such that the content of the matrix is copied.
     * Else, the pointers would be set equal.*/
    gsl_matrix_memcpy(hilf,tobeinverted);
    int signum;

    gsl_permutation * perm = gsl_permutation_alloc(dim);
    gsl_linalg_LU_decomp(hilf, perm, &signum);
    gsl_linalg_LU_invert(hilf,perm,inverted);

    gsl_matrix_free(hilf);
    gsl_permutation_free(perm);
    /* cout << "Inverted Matrix:" << endl;
     printmatrix(inverted,dim,dim);*/

}



double determinant(gsl_matrix* f, unsigned int dim)
{
    /*returns the determinant of the square matrix f, when dim is the dimension of f*/

    gsl_matrix* hilf = gsl_matrix_alloc(dim,dim);
    gsl_matrix_memcpy(hilf,f); //use memcopy, to not delete the external matrix f
    int signum;

    gsl_permutation * perm = gsl_permutation_alloc(dim);
    gsl_linalg_LU_decomp(hilf, perm, &signum);
    double det = gsl_linalg_LU_det(hilf,signum);

    gsl_permutation_free(perm);
    gsl_matrix_free(hilf);
    return det;
}







/*scalar product without gsl_vector*/
/*assuming M is square with size dim times dim*/
double vMv(vector<double> vec1,gsl_matrix* M, vector<double> vec2,int dim )
{
    double res = 0;

    for(int i = 0; i < dim; i++)
    {
        for(int j = 0; j < dim; j++)
        {
            res+=vec1[i]*gsl_matrix_get(M, i, j)*vec2[j];
        }
    }


    return res;
}


double trace(gsl_matrix* m, int dim )
{
    double tr = 0.0;

    for(int i = 0; i < dim; i++)
    {
        tr += gsl_matrix_get(m,i,i);
    }

    return tr;
}

/*multiplies four matrices*/
void mult4mat(gsl_matrix* m1, gsl_matrix* m2, gsl_matrix* m3, gsl_matrix* m4, int dim, gsl_matrix* result)
{

    gsl_matrix* intermed1 = gsl_matrix_calloc(dim,dim);
    gsl_matrix* intermed2 = gsl_matrix_calloc(dim,dim);

    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,m3,m4,0.0,intermed1);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,m2,intermed1,0.0,intermed2);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,m1,intermed2,0.0,result);

    gsl_matrix_free(intermed1);
    gsl_matrix_free(intermed2);
}



void mult3mat(gsl_matrix* m1, gsl_matrix* m2, gsl_matrix* m3, int dim, gsl_matrix* result)
{

    gsl_matrix* intermed1 = gsl_matrix_calloc(dim,dim);

    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,m2,m3,0.0,intermed1);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,m1,intermed1,0.0,result);

    gsl_matrix_free(intermed1);
}

/*multiplying two matrices, and storing the result in result. The matrix result must be *externally* allocated.*/
void mult2mat(gsl_matrix* m1, gsl_matrix* m2, int dim, gsl_matrix* result)
{
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,m1,m2,0.0,result);

}

/*flag alldiag shall be true if m1 and m2 are diagonal matrices*/
void mult2mat(gsl_matrix* m1, gsl_matrix* m2, int dim, gsl_matrix* result, bool alldiag)
{
  if(!alldiag)
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,m1,m2,0.0,result);
  
  
  if(alldiag)
  {
    gsl_matrix_set_identity(result);
    for(int i = 0; i < dim; i++)
    {
      gsl_matrix_set(result,i,i,gsl_matrix_get(m1,i,i)*gsl_matrix_get(m2,i,i));
    }
  }

}


/*subfish must be an externally allocated 2 by 2 matrix*/
void retSubFish(int p1, int p2, int dimbigfish, gsl_matrix* bigfish, gsl_matrix* subfish)
{
    /*the parameters p1 and p2 get selected, the others are cast away*/

    gsl_vector* transporig = gsl_vector_alloc(dimbigfish);
    gsl_vector* transp2 = gsl_vector_alloc(2);

    gsl_matrix * rectangle = gsl_matrix_alloc(2,dimbigfish); //zeilen,spalten

    gsl_matrix_get_row(transporig, bigfish,p1);
    gsl_matrix_set_row(rectangle,0,transporig);
    gsl_matrix_get_row(transporig,bigfish,p2);
    gsl_matrix_set_row(rectangle,1,transporig);

    gsl_matrix_get_col(transp2,rectangle,p1);
    gsl_matrix_set_col(subfish,0,transp2);
    gsl_matrix_get_col(transp2,rectangle,p2);
    gsl_matrix_set_col(subfish,1,transp2);

    gsl_vector_free(transporig);
    gsl_vector_free(transp2);
    gsl_matrix_free(rectangle);

}




void Gnuplotconfidencecontours(string specifier, vector<double> testprob)
{

    vector<double> pdf(testprob);

    double sum = 0;


    double sumsig1 = 0;
    double sumsig2 = 0;
    double sumsig3 = 0;

    double contsig1 = 0;
    double contsig2 = 0;
    double contsig3 = 0;

    for(int i= 0; i < pdf.size(); i++)
    {
        sum += pdf[i];
    }
    // cout << "Sum: " << sum << endl;
    sort(pdf.begin(), pdf.end());
    // cout << pdf.size() << endl;

    /*.size-1 weil ein array ab null zählt*/

    for(int i = pdf.size() - 1; i >= 0; i--)
    {
        sumsig1 += pdf[i];
        contsig1 = pdf[i];

        if (sumsig1 > sum * 0.68)break; //0.6823

    }


    for(int i = pdf.size() - 1; i >= 0; i--)
    {
        sumsig2 += pdf[i];
        contsig2 = pdf[i];

        if (sumsig2 > sum * 0.90)break; //0.954

    }


    for(int i = pdf.size() - 1; i >= 0; i--)
    {
        sumsig3 += pdf[i];
        contsig3 = pdf[i];

        if (sumsig3 > sum * 0.95)break; //0.9995
        /*Wenn man hier stattdessen nur 95 Prozent wie bei Wolz nimmt, sind die Konturen viel schmaler als beim SNIa-paper. */

    }


    cout << specifier << "set cntrparam levels discrete " << contsig1 << ", " << contsig2 << ", " << contsig3 << endl;
    cout << "  the largest probability encountered, is: " << pdf[pdf.size()-1] << endl;

}





void pythoncontours(vector<double> testprob, int xcount, int ycount, string executablename,string datafile, string xaxis, string yaxis)
{

    vector<double> pdf(testprob);

    double sum = 0;

    double sumsig1 = 0;
    double sumsig2 = 0;
    double sumsig3 = 0;

    double contsig1 = 0;
    double contsig2 = 0;
    double contsig3 = 0;

    for(int i= 0; i < pdf.size(); i++)
    {
        sum += pdf[i];
    }
    // cout << "Sum: " << sum << endl;
    sort(pdf.begin(), pdf.end());


    /*.size-1 weil ein array ab null zählt*/

    for(int i = pdf.size() - 1; i >= 0; i--)
    {
        sumsig1 += pdf[i];
        contsig1 = pdf[i];

        if (sumsig1 > sum * 0.68)break; //0.6823

    }


    for(int i = pdf.size() - 1; i >= 0; i--)
    {
        sumsig2 += pdf[i];
        contsig2 = pdf[i];

        if (sumsig2 > sum * 0.90)break; //0.954

    }


    for(int i = pdf.size() - 1; i >= 0; i--)
    {
        sumsig3 += pdf[i];
        contsig3 = pdf[i];

        if (sumsig3 > sum * 0.95)break; //0.9995

    }


    FILE * pythonex;
    pythonex = fopen(executablename.c_str(), "w");

    fprintf(pythonex, "import numpy as np\nimport matplotlib.pyplot as plt\nfrom matplotlib import rc\n\n" );
    fprintf(pythonex, "font ={'size': 20}\nplt.rc('font',**font)\nplt.rc('text',usetex=True)\norigin = 'lower'\n\n");

    fprintf(pythonex, "N2x=%.d\nN2y=%.d\n",xcount,ycount);
    fprintf(pythonex, "x2,y2,z2 = np.genfromtxt(r'%s', unpack=True)\n",datafile.c_str());

    fprintf(pythonex, "xi2 = np.reshape(x2,(-1,N2x))\nyi2 = np.reshape(y2,(-1,N2y))\nzi2 = np.reshape(z2,(-1,N2x))\n");
    fprintf(pythonex, "print \"plotting...\"\n");
    fprintf(pythonex, "levels2 = [%.3e,%.3e,%.3e]\n",contsig1, contsig2, contsig3);
    fprintf(pythonex, "CS2 = plt.contour(xi2, yi2, zi2, levels2, colors = ('navy','steelblue'),linewidths = (2),origin = origin)\n");
    fprintf(pythonex, "plt.xlabel(r'%s')\nplt.ylabel(r'%s')\n",xaxis.c_str(),yaxis.c_str());
    fprintf(pythonex, "\n");


    string PDFName = datafile.c_str() +(string)".pdf";
    fprintf(pythonex, "plt.savefig('%s')\n",PDFName.c_str());

    fclose(pythonex);

}



void pythoncontours_w(vector<double> testprob, int xcount, int ycount, string executablename,string datafile, string xaxis, string yaxis)
{

    vector<double> pdf(testprob);

    double sum = 0;

    double sumsig1 = 0;
    double sumsig2 = 0;
    double sumsig3 = 0;

    double contsig1 = 0;
    double contsig2 = 0;
    double contsig3 = 0;

    for(int i= 0; i < pdf.size(); i++)
    {
        sum += pdf[i];
    }
    // cout << "Sum: " << sum << endl;
    sort(pdf.begin(), pdf.end());


    /*.size-1 weil ein array ab null zählt*/

    for(int i = pdf.size() - 1; i >= 0; i--)
    {
        sumsig1 += pdf[i];
        contsig1 = pdf[i];

        if (sumsig1 > sum * 0.68)break; //0.6823

    }


    for(int i = pdf.size() - 1; i >= 0; i--)
    {
        sumsig2 += pdf[i];
        contsig2 = pdf[i];

        if (sumsig2 > sum * 0.90)break; //0.954

    }


    for(int i = pdf.size() - 1; i >= 0; i--)
    {
        sumsig3 += pdf[i];
        contsig3 = pdf[i];

        if (sumsig3 > sum * 0.95)break; //0.9995

    }


    FILE * pythonex;
    pythonex = fopen(executablename.c_str(), "w");

    fprintf(pythonex, "import numpy as np\nimport matplotlib.pyplot as plt\nfrom matplotlib import rc\n\n" );
    fprintf(pythonex, "font ={'size': 20}\nplt.rc('font',**font)\nplt.rc('text',usetex=True)\norigin = 'lower'\n\n");

    fprintf(pythonex, "N2x=%.d\nN2y=%.d\n",xcount,ycount);
    fprintf(pythonex, "x2,y2,z2 = np.genfromtxt(r'%s', unpack=True)\n",datafile.c_str());

    fprintf(pythonex, "xi2 = np.reshape(x2,(-1,N2x))\nyi2 = np.reshape(y2,(-1,N2y))\nzi2 = np.reshape(z2,(-1,N2x))\n");
    fprintf(pythonex, "print \"plotting...\"\n");
    fprintf(pythonex, "levels2 = [%.3e,%.3e,%.3e]\n",contsig1, contsig2, contsig3);
    fprintf(pythonex, "CS2 = plt.contour(xi2, yi2, zi2, levels2, colors = ('lightgrey','WhiteSmoke'),linewidths = (2),origin = origin)\n");
    fprintf(pythonex, "plt.xlabel(r'%s')\nplt.ylabel(r'%s')\n",xaxis.c_str(),yaxis.c_str());
    fprintf(pythonex, "\n");


    string PDFName = datafile.c_str() +(string)".pdf";
    fprintf(pythonex, "plt.savefig('%s')\n",PDFName.c_str());

    fclose(pythonex);

}







