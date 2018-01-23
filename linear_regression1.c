<<<<<<< HEAD


=======
//
>>>>>>> 67a88c6330de743a3cc3020c79b555decd5795d4
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include <string.h>

#define true 1
#define false 0
#define THRESHOLD 2/sqrt(150)
#define NUM_OBSERVATION 150 /* read 1024 bytes at a time */
#define NUM_OBSERVATIONS 150
#define TERMS_TO_PREDICT 18
//#define FILE_TO_READ "series.csv"
//#define FILE_TO_WRITE "series_out.csv"

double calcMean(double arr[],int length)
{
	double sum=0.0;
	for(int i=0;i<length;i++)
	{
		sum+=arr[i];
	}
	return (sum/length);
}

int calcArrMin(int a[], int n) {
  int c, min, index;
 
  min = a[0];
  index = 0;
 
  for (c = 1; c < n; c++) {
    if (a[c] < min) {
       index = c;
       min = a[c];
    }
  }
 
  return index;
}

void Duplicate(double arr_x[],double arr_y[],int length){
	for(int i=0;i<length;i++){
		arr_y[i]=arr_x[i];
	}
}

double * squareArr(double arr[],int size)
{
	double newArr[size];
	for(int i=0;i<size;i++)
        {
                newArr[i]= arr[i]*arr[i];
        }
        return newArr;
}

double calcArrSum(double arr[],int length)
{
	double sum=0;
	for(int i=0;i<length;i++)
	{
		sum+=arr[i];
	}
	return sum;
}

void calcArrDiff(double arr1[],double  arr2[],double arr_diff[],int size)
{
	for(int i=0;i<size;i++)
	{
		arr_diff[i]= arr1[i]-arr2[i];
	}
	return;
}

void calcMDArrDiff(int m,int n,double arr1[][n],double arr2[][n],double arr_diff[][n]){
		for(int i=0;i<m;i++){
			for(int j=0;j<n;j++){
				arr_diff[i][j]=arr1[i][j]-arr2[i][j];
			}
		}
	return;
}

double * calcArrProd(double arr1[],double arr2[],int size)
{
	double newArr[size];
	for(int i=0;i<size;i++)
	{
		newArr[i]=arr1[i]*arr2[i];
	}
	return newArr;
}

double calcSD(double arr[],int length)
{
	double avg=calcMean(arr,length);
	double newArr[length],*squaredArr;
	for(int i=0;i<length;i++)
	{
		newArr[i]=arr[i]-avg;
	}
	squaredArr=squareArr(newArr,length);
	double sum=calcArrSum(squaredArr,length); 
	double variance=sum/(length-1);
	
	return sqrt(variance);
}
void calcArrLead(double arr[],double arr_lead[],int length,int lag){
	
	for(int i=0;i<(length-lag);i++)
	{
		arr_lead[i]=arr[i+lag];
	}
	return;
}

void calcArrLag(double arr[],double arr_lag[],int length,int lag){
	for(int i=0;i<(length-lag);i++)
	{
		arr_lag[i]=arr[i];
		/*printf("%lf\n",arr_lag[i])*/;
	}
	return;
}

void calcArrLeadLag(double arr[],double arr_lead[],int length,int lead,int lag){

        for(int i=0;i<(length-lead-lag);i++)
        {
                arr_lead[i]=arr[i+lag];
        }
        return;
}



void calc2DArrMean(int m,int n,double arr[][n],double mean_arr[][n]){
        double x;
        for(int i=0;i<n;i++){
                double tmp_arr[m];
                for(int j=0;j<m;j++){
                        tmp_arr[j]=arr[j][i];
                }
                mean_arr[0][i]=calcMean(tmp_arr,m);

        }
        return;
}

void modArrNorm(int m,int n,double arr[][n],double norm_arr[][n]){
        double mean_arr[1][n];
        calc2DArrMean(m,n,arr,mean_arr);
        for(int i=0;i<m;i++){
                for(int j=0;j<n;j++){
                        norm_arr[i][j]=(arr[i][j]-mean_arr[0][j]);
                }
        }
        return;
}

void Cbind(int m,int n,double arr1[],double arr2[],double arr3[],double cbind_arr[][n]){
		
		for(int i=0;i<m;i++){
			cbind_arr[i][0]=arr1[i];
		}
		
		for(int i=0;i<m;i++){
			cbind_arr[i][1]=arr2[i];
		}
		
		for(int i=0;i<m;i++){
			cbind_arr[i][2]=arr3[i];
		}
		return;
}

void Cbind4(int m,int n,double arr1[],double arr2[],double arr3[],double arr4[],double cbind_arr[][n]){

                for(int i=0;i<m;i++){
                        cbind_arr[i][0]=arr1[i];
                }

                for(int i=0;i<m;i++){
                        cbind_arr[i][1]=arr2[i];
                }

                for(int i=0;i<m;i++){
                        cbind_arr[i][2]=arr3[i];
                }
		
		for(int i=0;i<m;i++){
                        cbind_arr[i][3]=arr4[i];
                }
                return;
}

void D2MD(int m,double arr[],double arr_2d[][1]){
		for(int i=0;i<m;i++){
				arr_2d[i][0]=arr[i];
		}
		return;
}

void MD2D(int m,double arr[][1],double arr_1d[]){
		for(int i=0;i<m;i++){
				arr_1d[i]=arr[i][0];
		}
		return;
}


void inverse(int n,double a[][n],double inv_arr[][n],int size){

  double determinant=0;
  for(int i=0;i<n;i++){
      determinant = determinant + (a[0][i]*(a[1][(i+1)%3]*a[2][(i+2)%3] - a[1][(i+2)%3]*a[2][(i+1)%3]));
  }
   /*printf("\nInverse of matrix is: \n\n");*/
   for(int i=0;i<n;i++){
      for(int j=0;j<n;j++){
        inv_arr[j][i]=((a[(i+1)%3][(j+1)%3] * a[(i+2)%3][(j+2)%3]) - (a[(i+1)%3][(j+2)%3]*a[(i+2)%3][(j+1)%3]))/ determinant;
        }
   }


   return;
}

void inverse4(double m[4][4],double inv[4][4]){

    double det;
    int i;

    inv[0][0] = m[1][1]  * m[2][2] * m[3][3] -
             m[1][1]  * m[2][3] * m[3][2] -
             m[2][1]  * m[1][2]  * m[3][3] +
             m[2][1]  * m[1][3]  * m[3][2] +
             m[3][1] * m[1][2]  * m[2][3] -
             m[3][1] * m[1][3]  * m[2][2];

    inv[1][0] = -m[1][0]  * m[2][2] * m[3][3] +
              m[1][0]  * m[2][3] * m[3][2] +
              m[2][0]  * m[1][2]  * m[3][3] -
              m[2][0]  * m[1][3]  * m[3][2] -
              m[3][0] * m[1][2]  * m[2][3] +
              m[3][0] * m[1][3]  * m[2][2];

    inv[2][0] = m[1][0]  * m[2][1] * m[3][3] -
             m[1][0]  * m[2][3] * m[3][1] -
             m[2][0]  * m[1][1] * m[3][3] +
             m[2][0]  * m[1][3] * m[3][1] +
             m[3][0] * m[1][1] * m[2][3] -
             m[3][0] * m[1][3] * m[2][1];

    inv[3][0] = -m[1][0]  * m[2][1] * m[3][2] +
               m[1][0]  * m[2][2] * m[3][1] +
               m[2][0]  * m[1][1] * m[3][2] -
               m[2][0]  * m[1][2] * m[3][1] -
               m[3][0] * m[1][1] * m[2][2] +
               m[3][0] * m[1][2] * m[2][1];

    inv[0][1] = -m[0][1]  * m[2][2] * m[3][3] +
              m[0][1]  * m[2][3] * m[3][2] +
              m[2][1]  * m[0][2] * m[3][3] -
              m[2][1]  * m[0][3] * m[3][2] -
              m[3][1] * m[0][2] * m[2][3] +
              m[3][1] * m[0][3] * m[2][2];

    inv[1][1] = m[0][0]  * m[2][2] * m[3][3] -
             m[0][0]  * m[2][3] * m[3][2] -
             m[2][0]  * m[0][2] * m[3][3] +
             m[2][0]  * m[0][3] * m[3][2] +
             m[3][0] * m[0][2] * m[2][3] -
             m[3][0] * m[0][3] * m[2][2];

    inv[2][1] = -m[0][0]  * m[2][1] * m[3][3] +
              m[0][0]  * m[2][3] * m[3][1] +
              m[2][0]  * m[0][1] * m[3][3] -
              m[2][0]  * m[0][3] * m[3][1] -
              m[3][0] * m[0][1] * m[2][3] +
              m[3][0] * m[0][3] * m[2][1];

    inv[3][1] = m[0][0]  * m[2][1] * m[3][2] -
              m[0][0]  * m[2][2] * m[3][1] -
              m[2][0]  * m[0][1] * m[3][2] +
              m[2][0]  * m[0][2] * m[3][1] +
              m[3][0] * m[0][1] * m[2][2] -
              m[3][0] * m[0][2] * m[2][1];

    inv[0][2] = m[0][1]  * m[1][2] * m[3][3] -
             m[0][1]  * m[1][3] * m[3][2] -
             m[1][1]  * m[0][2] * m[3][3] +
             m[1][1]  * m[0][3] * m[3][2] +
             m[3][1] * m[0][2] * m[1][3] -
             m[3][1] * m[0][3] * m[1][2];

    inv[1][2] = -m[0][0]  * m[1][2] * m[3][3] +
              m[0][0]  * m[1][3] * m[3][2] +
              m[1][0]  * m[0][2] * m[3][3] -
              m[1][0]  * m[0][3] * m[3][2] -
              m[3][0] * m[0][2] * m[1][3] +
              m[3][0] * m[0][3] * m[1][2];

    inv[2][2] = m[0][0]  * m[1][1] * m[3][3] -
              m[0][0]  * m[1][3] * m[3][1] -
              m[1][0]  * m[0][1] * m[3][3] +
              m[1][0]  * m[0][3] * m[3][1] +
              m[3][0] * m[0][1] * m[1][3] -
              m[3][0] * m[0][3] * m[1][1];

    inv[3][2] = -m[0][0]  * m[1][1] * m[3][2] +
               m[0][0]  * m[1][2] * m[3][1] +
               m[1][0]  * m[0][1] * m[3][2] -
               m[1][0]  * m[0][2] * m[3][1] -
               m[3][0] * m[0][1] * m[1][2] +
               m[3][0] * m[0][2] * m[1][1];

    inv[0][3] = -m[0][1] * m[1][2] * m[2][3] +
              m[0][1] * m[1][3] * m[2][2] +
              m[1][1] * m[0][2] * m[2][3] -
              m[1][1] * m[0][3] * m[2][2] -
              m[2][1] * m[0][2] * m[1][3] +
              m[2][1] * m[0][3] * m[1][2];

    inv[1][3] = m[0][0] * m[1][2] * m[2][3] -
             m[0][0] * m[1][3] * m[2][2] -
             m[1][0] * m[0][2] * m[2][3] +
             m[1][0] * m[0][3] * m[2][2] +
             m[2][0] * m[0][2] * m[1][3] -
             m[2][0] * m[0][3] * m[1][2];

    inv[2][3] = -m[0][0] * m[1][1] * m[2][3] +
               m[0][0] * m[1][3] * m[2][1] +
               m[1][0] * m[0][1] * m[2][3] -
               m[1][0] * m[0][3] * m[2][1] -
               m[2][0] * m[0][1] * m[1][3] +
               m[2][0] * m[0][3] * m[1][1];
			   
    inv[3][3] = m[0][0] * m[1][1] * m[2][2] -
              m[0][0] * m[1][2] * m[2][1] -
              m[1][0] * m[0][1] * m[2][2] +
              m[1][0] * m[0][2] * m[2][1] +
              m[2][0] * m[0][1] * m[1][2] -
              m[2][0] * m[0][2] * m[1][1];

    det = m[0][0] * inv[0][0] + m[0][1] * inv[1][0] + m[0][2] * inv[2][0] + m[0][3] * inv[3][0];


    det = 1.0 / det;
        for(int i=0;i<4;i++){
            for (int j = 0; j < 4; j++){
                inv[i][j] = inv[i][j] * det;
                /*printf("%lf\t",inv[i][j]);*/
                }
        /*printf("\n")*/;
        }

    return;
}

void transpose(int m,int n, double arr[][n],double t_arr[][m],int size){

        for(int i=0;i<n;i++){
                for(int j=0;j<m;j++){
                        t_arr[i][j]=arr[j][i];
                }
        }
        return;
}


void product3(int m1,int n1,int m2,int n2,double arr1[][n1],double arr2[][n2],double p_arr[][n2],int size){
        double sum=0.0;
        for(int i=0;i<m1;i++){
                for(int j=0;j<n2;j++){
                        for(int k=0;k<n1;k++){
                                sum=sum+arr1[i][k]*arr2[k][j];
                        }
                p_arr[i][j]=sum;
                sum=0.0;
                }
        }

        return;
}


double *  LR1(double arr_x[],double arr_y[],int length){
	double arr_x_diff[length];
	double arr_y_diff[length];
	double mean_x,mean_y;
	double Sx,Sy,Corr;
	double *estimates = malloc(sizeof(double)*2);
	double beta,intercept;

	mean_x=calcMean(arr_x,length);
	mean_y=calcMean(arr_y,length);
	for(int i=0;i<length;i++)
	{
		arr_x_diff[i]=arr_x[i]-mean_x;
		arr_y_diff[i]=arr_y[i]-mean_y;
	}
	
	Sx=calcSD(arr_x,length);
	Sy=calcSD(arr_y,length);
	Corr=calcArrSum(calcArrProd(arr_x_diff,arr_y_diff,length),length)/((length-1)*Sx*Sy);
	beta=Corr*Sy/Sx;
	intercept=mean_y-(beta*mean_x);
	estimates[0]=beta;
	estimates[1]=intercept;
	printf("Mean_x: %lf\tMean_y: %lf\tCorr: %lf\tSx: %lf\tSy: %lf",mean_x,mean_y,Corr,Sx,Sy);
	return estimates;
}

double * LR2(double arr_x1[],double arr_x2[],double arr_y[],int length){
	double arr_x1_diff[length];
	double arr_x2_diff[length];
	double arr_y_diff[length];
	double beta1,beta2,intercept;
	double mean_x1,mean_x2,mean_y;
	double Sx1,Sx2,Sx1y,Sx2y,Sx1x2;
	double SumX1,SumX2,SumY,SumX1Y,SumX2Y,SumX1X2;
	double *estimates = malloc(sizeof(double)*3);;
	
	mean_x1=calcMean(arr_x1,length);
	mean_x2=calcMean(arr_x2,length);
	mean_y=calcMean(arr_y,length);
	
	for(int i=0;i<length;i++)
	{
		arr_x1_diff[i]=arr_x1[i]-mean_x1;
		arr_x2_diff[i]=arr_x2[i]-mean_x2;
		arr_y_diff[i]=arr_y[i]-mean_y;
	}

	SumX1=calcArrSum(arr_x1,length);
	SumX2=calcArrSum(arr_x2,length);
	SumY=calcArrSum(arr_y,length);
	SumX1Y=calcArrSum(calcArrProd(arr_x1,arr_y,length),length);
	SumX2Y=calcArrSum(calcArrProd(arr_x2,arr_y,length),length);
	SumX1X2=calcArrSum(calcArrProd(arr_x1,arr_x2,length),length);

	Sx1=calcArrSum(squareArr(arr_x1_diff,length),length);
	Sx2=calcArrSum(squareArr(arr_x2_diff,length),length);
	Sx1y=SumX1Y-((SumX1*SumY)/length);
	Sx2y=SumX2Y-((SumX2*SumY)/length);
	Sx1x2=SumX1X2-((SumX1*SumX2)/length);

	
	beta1=((Sx2*Sx1y)-(Sx1x2*Sx2y))/((Sx1*Sx2)-(Sx1x2*Sx1x2));
	beta2=((Sx1*Sx2y)-(Sx1x2*Sx1y))/((Sx1*Sx2)-(Sx1x2*Sx1x2));
	intercept=mean_y-(beta1*mean_x1)-(beta2*mean_x2);
	estimates[0]=beta1;
	estimates[1]=beta2;
	estimates[2]=intercept;
	return estimates;
}

double *  LRM(int m,int n,double X[][n],double Y[][1]){
	int size=m;
	double Xn[m][n],Yn[m][1],X_bar[1][n],Y_bar[1][1],Xt[n][m],XX[n][n],XXi[n][n],XXX[n][m],XY[m][1];
        double beta[n][1],intercept,sum;
	double *estimates = malloc(sizeof(double)*(n+1));
	modArrNorm(m,n,X,Xn);
        modArrNorm(m,1,Y,Yn);
        transpose(m,n,Xn,Xt,size);
        product3(n,m,m,n,Xt,Xn,XX,size);
        if(n==3){
		inverse(3,XX,XXi,size);
	}else{
		inverse4(XX,XXi);
	}
        product3(n,n,n,m,XXi,Xt,XXX,size);
        product3(n,m,m,1,XXX,Yn,beta,size);
	for(int i=0;i<n;i++){
		/*printf("%lf\n",beta[i][0]);*/
		estimates[i]=beta[i][0];
	}

	/*Intercept Calculation*/
	calc2DArrMean(m,n,X,X_bar);
	calc2DArrMean(m,1,Y,Y_bar);
	sum=0.0;
	for(int i=0;i<n;i++){
		sum=sum+(estimates[i]*X_bar[0][i]);
	}
	intercept=Y_bar[0][0]-sum;
	estimates[n]=intercept;
	/*printf("%f\n",intercept);*/
	return estimates;

}
void LR1_Pred(double arr[],double arr_pred[],double beta,double intercept,int length){
	for(int i=0;i<length;i++)
	{
		arr_pred[i]=(arr[i]*beta)+intercept;
	}
	return;
}


void LR2_Pred(double arr_x1[],double arr_x2[],double arr_pred[],double beta1,double beta2,double intercept,int length){
	for(int i=0;i<length;i++)
	{
		arr_pred[i]=(arr_x1[i]*beta1)+(arr_x2[i]*beta2)+intercept;
	}
	return;
}

void LRM_Pred(int m,int n,double X[][n],double Y_pred[][1],double estimates[n]){
	
	for(int i=0;i<m;i++){
		double sum=0.0;
		for(int j=0;j<n;j++){
			
               		sum=sum+(X[i][j]*estimates[j]);
			
		}
	Y_pred[i][0]=sum+estimates[n];
	}
       	return;
}	


double Corr(double arr_x[],double arr_y[],int length){
		double arr_x_diff[length];
        double arr_y_diff[length];
        double mean_x,mean_y;
        double Sx,Sy,Corr;

        mean_x=calcMean(arr_x,length);
        mean_y=calcMean(arr_y,length);
        for(int i=0;i<length;i++)
        {
                arr_x_diff[i]=arr_x[i]-mean_x;
                arr_y_diff[i]=arr_y[i]-mean_y;
        }

        Sx=calcSD(arr_x,length);
        Sy=calcSD(arr_y,length);
        Corr=calcArrSum(calcArrProd(arr_x_diff,arr_y_diff,length),length)/((length-1)*Sx*Sy);
        return Corr;
}

int DFTest(double arr[],double arr_recov[],int length){
        int d=0;
        int dd=0;
        double * estimates;
        double Y_t[length-1],Y_t_1[length-1];

        do{

               calcArrLead(arr,Y_t,length,1);
               calcArrLag(arr,Y_t_1,length,1);
	       /*printf("\nLength_DFT: %d\n",length);*/
               estimates=LR1(Y_t_1,Y_t,length);
			
		if(d>0){
			for(int i=0;i<(length-dd);i++){
					arr[i]=(dd*arr[i+dd])-arr[i];
					/*printf("Value:%lf\n",arr[i]);*/
			}
		length=length-1;
		}
	      arr_recov[d]=arr[length-1];
              d=d+1;
              dd=1;
	      /*printf("DFT: %lf\n",estimates[0]);*/
}while(!(estimates[0]<=1.001 && estimates[0]>=-1.001));
        free(estimates);
        return (d-1);
}


void Drift(double arr[],double arr_stry[],int length,int d){

	for(int i=0;i<d;i++){
		for(int j=0;i<(length-d);i++){
			arr_stry[j]=arr_stry[j+d]-arr_stry[j];
		}
	}
}


void EAFMatrix(double arr[],double arr_eaf[][3],int length){
        	/*p=0*/
		
		double Y_t[length-3],Y_t_1[length-3],Y_t_2[length-3],Y_t_3[length-3],Y_t_cap[length-3],Y_error[length-3];
		double phi,c,phi1,phi2;
		double E_t[length-6],E_t_1[length-6],E_t_2[length-6],E_t_3[length-6];
		double *AR1_estimates,*AR2_estimates;
		
		calcArrLeadLag(arr,Y_t,length,0,3);
        	calcArrLeadLag(arr,Y_t_1,length,1,2);
        	calcArrLeadLag(arr,Y_t_2,length,2,1);
		calcArrLeadLag(arr,Y_t_3,length,3,0);
		length=length-3;
		
		arr_eaf[0][0]=Corr(Y_t,Y_t_1,length);
		arr_eaf[0][1]=Corr(Y_t,Y_t_2,length);
		arr_eaf[0][2]=Corr(Y_t,Y_t_3,length);
		
		/*P=1*/
		/*printf("\nEAF_Length: %d\n",length);*/
		AR1_estimates=LR1(Y_t_1,Y_t,length);
        	phi=AR1_estimates[0];
        	c=AR1_estimates[1];
       
	        LR1_Pred(Y_t_1,Y_t_cap,phi,c,length);
        	calcArrDiff(Y_t,Y_t_cap,Y_error,length);
		
		calcArrLeadLag(Y_error,E_t,length,0,3);
	        calcArrLeadLag(Y_error,E_t_1,length,1,2);
        	calcArrLeadLag(Y_error,E_t_2,length,2,1);
		calcArrLeadLag(Y_error,E_t_3,length,3,0);
		
		/*Change Length here*/
		arr_eaf[1][0]= Corr(E_t,E_t_1,length-3);
		arr_eaf[1][1]= Corr(E_t,E_t_2,length-3);
		arr_eaf[1][2]= Corr(E_t,E_t_3,length-3);
		
		/*P=2*/
		AR2_estimates=LR2(Y_t_1,Y_t_2,Y_t,length);
	        phi1=AR2_estimates[0];
		phi2=AR2_estimates[1];
        	c=AR2_estimates[2];
      		
	        LR2_Pred(Y_t_1,Y_t_2,Y_t_cap,phi1,phi2,c,length);
        	calcArrDiff(Y_t,Y_t_cap,Y_error,length);
		
		calcArrLeadLag(Y_error,E_t,length,0,3);
        	calcArrLeadLag(Y_error,E_t_1,length,1,2);
	        calcArrLeadLag(Y_error,E_t_2,length,2,1);
		calcArrLeadLag(Y_error,E_t_3,length,3,0);
		
		/*Change Length here*/
		arr_eaf[2][0]=Corr(E_t,E_t_1,length-3);
		arr_eaf[2][1]=Corr(E_t,E_t_2,length-3);
		arr_eaf[2][2]=Corr(E_t,E_t_3,length-3);
		/*printf("\nEAF Matrix\n");
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
				printf("%lf\t",arr_eaf[i][j]);
			}
			printf("\n");
		}*/
		printf("\n\n");
		return;
		
}

double MAPE(double arr_y[],double arr_y_cap[],int length){
	double arr_error_per[length];
	for(int i=0;i<length;i++){
		arr_error_per[i]=fabs((arr_y[i]-arr_y_cap[i])/arr_y[i]);
		}
		
		return (calcArrSum(arr_error_per,length)/length)*100;
}

double MAE(double arr_y[],double arr_y_cap[],int length){
        double arr_error[length];
        for(int i=0;i<length;i++){
                arr_error[i]=fabs((arr_y[i]-arr_y_cap[i]));
                }

                return (calcArrSum(arr_error,length)/length)*100;
}



double * AR1(double arr[],int length){
	double *estimates;
	double *forecast = malloc(sizeof(double)*18);
	double Y_n;
	double Y_t[length-1],Y_t_1[length-1],Y_cap[length-1];
	double mape;
	
	calcArrLead(arr,Y_t,length,1);
	calcArrLag(arr,Y_t_1,length,1);
	length=length-1;
	printf("\nThis is it\n");	
	estimates=LR1(Y_t_1,Y_t,length);

	for(int i=0;i<length;i++){
		Y_cap[i]=(Y_t_1[i]*estimates[0])+estimates[1];	
	}
	
	printf("Phi: %lf\tConstant:%lf\n",estimates[0],estimates[1]);
	
	mape=MAPE(Y_t,Y_cap,length);
	
	Y_n=Y_t[length-1];
	for(int i=0;i<16;i++){
		forecast[i]=(estimates[0]*Y_n)+estimates[1];
		Y_n=forecast[i];
	}
	free(estimates);
	printf("MAPE: %lf\t",mape);
	forecast[16]=mape;
	return forecast;
	
}

double *  AR2(double arr[],int length){
	double *estimates;double Y_n,Y_n_1,mape;
	double Y_t[length-2],Y_t_1[length-2],Y_t_2[length-2],Y_cap[length-2];
	double *forecast=malloc(sizeof(double)*18);

	calcArrLead(arr,Y_t,length,2);
	calcArrLag(arr,Y_t_1,length,1);
	calcArrLead(Y_t_1,Y_t_1,length,1);
	calcArrLag(arr,Y_t_2,length,2);
	length=length-2;

	estimates=LR2(Y_t_1,Y_t_2,Y_t,length-2);

	for(int i=0;i<length;i++){
		Y_cap[i]=(estimates[0]*Y_t_1[i])+(estimates[1]*Y_t_2[i])+estimates[2];
	}
	
	
	mape=MAPE(Y_t,Y_cap,length);
	printf("Phi: %lf\tPhi2: %lf\tConstant:%lf\n",estimates[0],estimates[1],estimates[2]);
	printf("MAPE: %lf\t",mape);	
	Y_n=Y_t[length-1],Y_n_1=Y_t[length-2];
	for(int i=0;i<16;i++){
                forecast[i]=(estimates[0]*Y_n)+(estimates[1]*Y_n_1)+estimates[2];
                Y_n_1=Y_n;
		Y_n=forecast[i];
	}
	free(estimates);
	forecast[16]=mape;
	return forecast;
}


double * MA1(double arr[],int length){
	double phi1,phi2,c1, * AR2_estimates,* MA2_estimates,E_n,mape;
        double theta,c2;
        double theta_new,c2_new;
        double theta_int,c2_int;
	double DY_t[length-2],DY_t_1[length-2],DY_t_1_tmp[length-1],DY_t_2[length-2],DY_t_cap[length-2],DY_error[length-2];
	double Y_t[length-3],Y_t_int[length-3],E_t_1[length-3],Y_t_cap[length-3],Y_t_cap_int[length-3],Y_error[length-3];
        double *forecast=malloc(sizeof(double)*3);;
	int length_int;

        calcArrLead(arr,DY_t,length,2);
        calcArrLead(arr,DY_t_1_tmp,length,1);
        calcArrLag(DY_t_1_tmp,DY_t_1,length,1);
        calcArrLag(arr,DY_t_2,length,2);

        length=length-2;
        AR2_estimates=LR2(DY_t_1,DY_t_2,DY_t,length);
        phi1=AR2_estimates[0];
        phi2=AR2_estimates[1];
        c1=AR2_estimates[2];
	free(AR2_estimates);
        /*printf("phi1:%lf\nphi2:%lf\ncon:%lf\n",phi1,phi2,c1);*/
        LR2_Pred(DY_t_1,DY_t_2,DY_t_cap,phi1,phi2,c1,length);
        calcArrDiff(DY_t,DY_t_cap,DY_error,length);
        /*for(int i=0;i<20;i++){
                printf("D:%lf\n",DY_error[i]);
        }*/
        calcArrLead(DY_t,Y_t,length,1);
        calcArrLag(DY_error,E_t_1,length,1);

        length=length-1;
        MA2_estimates=LR1(E_t_1,Y_t,length);
        theta=MA2_estimates[0];
        c2=MA2_estimates[1];
        /*printf("phi:%lf\ntheta:%lf\ncons:%lf\n",phi_cap,theta,c2);*/
 	free(MA2_estimates);
        theta_new=theta;
        c2_new=c2;
	
	theta_int=theta;
        c2_int=c2;
	length_int=length;
	Duplicate(Y_t_cap,Y_t_cap_int,length);
	Duplicate(Y_t,Y_t_int,length);
	
        int counter=0;
 do{

                
                theta=theta_new;
                c2=c2_new;
                /*printf("phi_cap%lf\ntheta:%lf\nc:%lf\n",phi_cap,theta,c2);*/
                LR1_Pred(E_t_1,Y_t_cap,theta,c2,length);
                calcArrDiff(Y_t,Y_t_cap,Y_error,length);

                /*for(int i=0;i<length;i++){
                        printf("Diff:%lf\n",Y_error[i]);
                }*/
                calcArrLead(Y_t,Y_t,length,1);
                calcArrLag(Y_error,E_t_1,length,1);

                /*for(int i=0;i<20;i++){
                        printf("Y:%lf Y_1:%lf E:%lf\n",Y_t[i],Y_t_1[i],E_t_1[i]);
                }*/
                length=length-1;
                /*printf("Length:%d\n",length);*/
                MA2_estimates=LR1(E_t_1,Y_t,length);

                
                theta_new=MA2_estimates[0];
                c2_new=MA2_estimates[1];
		free(MA2_estimates);
                counter=counter+1;
}while(!(fabs(theta-theta_new)>0.01 && counter>30));

        printf("NIteration: %d\nTheta: %lf\nConstant: %lf\n",counter,theta_new,c2_new);
	
	if(counter>30){
		theta_new=theta_int;
		c2_new=c2_int;
		length=length_int;
	
		mape=MAPE(Y_t_int,Y_t_cap_int,length);
		};
	
	mape=MAPE(Y_t,Y_t_cap,length);
	printf("MAPE: %lf\t",mape);
	
	E_n=E_t_1[length-1];
        forecast[0]=(theta_new*E_n)+c2_new;
        forecast[1]=mape;
        return forecast;
}

double *  MA2(double arr[],int length){
        double phi1,phi2,phi3,c1,E_n,E_n_1,E_n_2;
        int n=3,length_int;
	double *forecast=malloc(sizeof(double)*3);
        double theta1,theta2,c2,mape;
        double theta1_new,theta2_new,c2_new;
	double theta1_int,theta2_int,c2_int;
        double *AR3_estimates,*MA2_estimates,*MA2_estimates_new;

        double DY_t[length-3],DY_t_1[length-3],DY_t_2[length-3],DY_t_3[length-3],DY_t_cap[length-2],DY_t_cap_int[length-2],DY_error[length-2];
        double Y_t[length-5],Y_t_int[length-5],Y_t_1[length-5],Y_t_2[length-5],Y_t_cap[length-5],Y_t_cap_int[length-5],Y_error[length-5];

	calcArrLeadLag(arr,DY_t,length,0,3);
        calcArrLeadLag(arr,DY_t_1,length,1,2);
        calcArrLeadLag(arr,DY_t_2,length,2,1);
        calcArrLeadLag(arr,DY_t_3,length,3,0);


        length=length-3;

        double DY_X[length][3];
        double DY_Y[length][1];
        double DY_Y_cap[length][1];
        double DY_Y_error[length][1];
        double E_t[length],E_t_1[length-2],E_t_2[length-2];

        D2MD(length,DY_t,DY_Y);
        Cbind(length,n,DY_t_1,DY_t_2,DY_t_3,DY_X);

        AR3_estimates=LRM(length,n,DY_X,DY_Y);

        phi1=AR3_estimates[0];
        phi2=AR3_estimates[1];
        phi3=AR3_estimates[2];
        c1=AR3_estimates[3];
	free(AR3_estimates);
        
        LRM_Pred(length,n,DY_X,DY_Y_cap,AR3_estimates);
        calcMDArrDiff(length,1,DY_Y,DY_Y_cap,DY_Y_error);

        MD2D(length,DY_Y_error,E_t);
        
	calcArrLeadLag(DY_t,Y_t,length,0,2);
        calcArrLeadLag(E_t,E_t_1,length,1,1);
        calcArrLeadLag(E_t,E_t_2,length,2,0);

        length=length-2;

	MA2_estimates=LR2(E_t_1,E_t_2,Y_t,length);

        theta1=MA2_estimates[0];
        theta2=MA2_estimates[1];
        c2=MA2_estimates[2];
	free(MA2_estimates);
        /*printf("phi:%lf\ntheta:%lf\ncons:%lf\n",phi_cap,theta,c2);*/
        
        theta1_new=theta1;
        theta2_new=theta2;
        c2_new=c2;
        

	theta1_int=theta1;
        theta2_int=theta2;
        c2_int=c2;
	length_int=length;
        LR2_Pred(E_t_1,E_t_2,Y_t_cap_int,theta1,theta2,c2,length);
	int counter=0;
        do{
                
                theta1=theta1_new;
                theta2=theta2_new;
                c2=c2_new;
                
		LR2_Pred(E_t_1,E_t_2,Y_t_cap,theta1,theta2,c2,length);
                calcArrDiff(Y_t,Y_t_cap,Y_error,length);
				
		calcArrLeadLag(Y_t,Y_t,length,0,2);
                calcArrLeadLag(Y_error,E_t_1,length,1,1);
		calcArrLeadLag(Y_error,E_t_2,length,2,0);

				length=length-2;
                /*printf("Length:%d\n",length);*/
                MA2_estimates=LR2(E_t_1,E_t_2,Y_t,length);

                theta1_new=MA2_estimates[0];
                theta2_new=MA2_estimates[1];
                c2_new=MA2_estimates[2];
		free(MA2_estimates);
                counter=counter+1;
		/*printf("%lf\t%lf\t%lf\n",theta1_new,theta2_new,c2_new);*/
}while(!(fabs(theta1-theta1_new) <0.01 && fabs(theta2-theta2_new) < 0.01 && counter>30));
        printf("%d",counter);
	if(counter>30){
		theta1_new=theta1_int;
		theta2_new=theta2_int;
		c2_new=c2_int;
		length=length_int;
                mape=MAPE(Y_t_int,Y_t_cap_int,length);
	};
	printf("Theta1: %lf\tTheta2: %lf\tConstant: %lf\n",theta1_new,theta2_new,c2_new);
	mape=MAPE(Y_t,Y_t_cap,length);
	printf("MAPE: \n%lf\n",mape);
	E_n=E_t_1[length-1],E_n_1=E_t_1[length-2];
        forecast[0]=((theta1_new*E_n)+(theta2_new*E_n_1)+c2_new);
	forecast[1]=mape;
        return forecast;
}

double *  AR1MA1(double arr[],int length){
	double phi1,phi2,c1, * AR2_estimates,* MA2_estimates,Y_n,Y_n_1,E_n;
	double *forecast=malloc(sizeof(double)*18);
	double phi_cap,theta,c2;
	double phi_cap_new,theta_new,c2_new;
	double phi_cap_int,theta_int,c2_int,mape;
	double Y_t[length-3],Y_t_int[length-3],Y_t_1[length-3],Y_t_1_tmp[length-2],E_t_1[length-3],Y_t_cap[length-3],Y_t_cap_int[length-3],Y_error[length-3];
	double DY_t[length-2],DY_t_1[length-2],DY_t_1_tmp[length-1],DY_t_2[length-2],DY_t_cap[length-2],DY_error[length-2];
	int length_int;
	
	calcArrLead(arr,DY_t,length,2);
	calcArrLead(arr,DY_t_1_tmp,length,1);
	calcArrLag(DY_t_1_tmp,DY_t_1,length,1);
	calcArrLag(arr,DY_t_2,length,2);

	length=length-2;	
	AR2_estimates=LR2(DY_t_1,DY_t_2,DY_t,length);
	phi1=AR2_estimates[0];
	phi2=AR2_estimates[1];
	c1=AR2_estimates[2];
	free(AR2_estimates);
	/*printf("phi1:%lf\nphi2:%lf\ncon:%lf\n",phi1,phi2,c1);*/
	LR2_Pred(DY_t_1,DY_t_2,DY_t_cap,phi1,phi2,c1,length);
	calcArrDiff(DY_t,DY_t_cap,DY_error,length);
	/*for(int i=0;i<20;i++){
		printf("D:%lf\n",DY_error[i]);
	}*/
	calcArrLead(DY_t,Y_t,length,1);
	calcArrLag(DY_t,Y_t_1,length,1);
	calcArrLag(DY_error,E_t_1,length,1);
 	
	length=length-1;	
    	MA2_estimates=LR2(Y_t_1,E_t_1,Y_t,length);
	phi_cap=MA2_estimates[0];
	theta=MA2_estimates[1];
	c2=MA2_estimates[2];
	free(MA2_estimates);
	/*printf("phi:%lf\ntheta:%lf\ncons:%lf\n",phi_cap,theta,c2);*/
	phi_cap_new=phi_cap;
	theta_new=theta;
	c2_new=c2;

	phi_cap_int=phi_cap;
        theta_int=theta;
        c2_int=c2;
	length_int=length;
	Duplicate(Y_t,Y_t_int,length);
	Duplicate(Y_t_cap,Y_t_cap_int,length);
	int counter=0;
	
	do{
		
		phi_cap=phi_cap_new;
		theta=theta_new;
		c2=c2_new;
		/*printf("phi_cap%lf\ntheta:%lf\nc:%lf\n",phi_cap,theta,c2);*/	
		LR2_Pred(Y_t_1,E_t_1,Y_t_cap,phi_cap,theta,c2,length);
		calcArrDiff(Y_t,Y_t_cap,Y_error,length);
		
		/*for(int i=0;i<length;i++){
			printf("Diff:%lf\n",Y_error[i]);
		}*/
		calcArrLead(Y_t,Y_t,length,1);
		calcArrLead(Y_t_1,Y_t_1,length,1);
		calcArrLag(Y_error,E_t_1,length,1);
		
		/*for(int i=0;i<20;i++){
			printf("Y:%lf Y_1:%lf E:%lf\n",Y_t[i],Y_t_1[i],E_t_1[i]);
		}*/
		length=length-1;
		/*printf("Length:%d\n",length);*/
		MA2_estimates=LR2(Y_t_1,E_t_1,Y_t,length);	
		
		phi_cap_new=MA2_estimates[0];
		theta_new=MA2_estimates[1];
		c2_new=MA2_estimates[2];
		free(MA2_estimates);
		counter=counter+1;
}while(!(fabs(phi_cap-phi_cap_new)<0.01 && fabs(theta-theta_new)<0.01 && counter>30));

	/*printf("NIteration:%d\nPhi:%lf\nTheta:%lf\nConstant:%lf\n",counter,phi_cap_new,theta_new,c2_new);*/
	Y_n=Y_t[length-1],E_n=E_t_1[length-1];
	printf("%lf\t%lf",Y_n,E_n);

	if(counter>30){
		phi_cap_new=phi_cap_int;
		theta_new=theta_int;
		c2_new=c2_int;
		length=length_int;
		mape=MAPE(Y_t_int,Y_t_cap_int,length);
	};
	printf("NIteration: %d\nPhi: %lf\nTheta: %lf\nConstant: %lf\n",counter,phi_cap_new,theta_new,c2_new);
	mape=MAPE(Y_t,Y_t_cap,length);
	printf("\nMAPE:%lf\n",mape);
        for(int i=0;i<16;i++){
                forecast[i]=(phi_cap_new*Y_n)+(theta_new*E_n)+c2_new;
                E_n=0;
                Y_n=forecast[i];
        }
	forecast[16]=mape;
        return forecast;	 
}

double * AR1MA2(double arr[],int length){
        double phi1,phi2,phi3,c1,Y_n,Y_n_1,E_n,E_n_1;
	int n=3,length_int;
	double *forecast=malloc(sizeof(double)*18);
	double phi_cap,theta1,theta2,c2;
        double phi_cap_new,theta1_new,theta2_new,c2_new;
	double phi_cap_int,theta1_int,theta2_int,c2_int,mape;
	double *AR3_estimates,*MA3_estimates,*MA3_estimates_new;
        
	double DY_t[length-3],DY_t_1[length-3],DY_t_2[length-3],DY_t_3[length-3],DY_t_cap[length-2],DY_error[length-2];
	double Y_t[length-5],Y_t_int[length-5],Y_t_1[length-5],Y_t_2[length-5],Y_t_cap[length-5],Y_t_cap_int[length-5],Y_error[length-5];
        

		
        calcArrLeadLag(arr,DY_t,length,0,3);
        calcArrLeadLag(arr,DY_t_1,length,1,2);
	calcArrLeadLag(arr,DY_t_2,length,2,1);
        calcArrLeadLag(arr,DY_t_3,length,3,0);
	
	
        length=length-3;
		
	double DY_X[length][3];
	double DY_Y[length][1];
	double DY_Y_cap[length][1];
	double DY_Y_error[length][1];
	double E_t[length],E_t_1[length-2],E_t_2[length-2];
		
	D2MD(length,DY_t,DY_Y);
	Cbind(length,n,DY_t_1,DY_t_2,DY_t_3,DY_X);
	
	AR3_estimates=LRM(length,n,DY_X,DY_Y);
		
        phi1=AR3_estimates[0];
        phi2=AR3_estimates[1];
	phi3=AR3_estimates[2];
	c1=AR3_estimates[3];
	free(AR3_estimates);
        
        /*printf("phi1:%lf\nphi2:%lf\ncon:%lf\n",phi1,phi2,c1);*/
        LRM_Pred(length,n,DY_X,DY_Y_cap,AR3_estimates);
        calcMDArrDiff(length,1,DY_Y,DY_Y_cap,DY_Y_error);
        /*for(int i=0;i<20;i++){
                printf("D:%lf\n",DY_error[i]);
        }*/
        MD2D(length,DY_Y_error,E_t);
	calcArrLeadLag(DY_t,Y_t,length,0,2);
        calcArrLeadLag(DY_t,Y_t_1,length,1,1);
	calcArrLeadLag(E_t,E_t_1,length,1,1);
        calcArrLeadLag(E_t,E_t_2,length,2,0);
        
	length=length-2;
		
	double Y_X[length][n];
	double Y_Y[length][1];
	double Y_Y_cap[length][1];
	double Y_Y_error[length][1];
		
	D2MD(length,Y_t,Y_Y);
	Cbind(length,3,Y_t_1,E_t_1,E_t_2,Y_X);
        MA3_estimates=LRM(length,n,Y_X,Y_Y);
        
	phi_cap=MA3_estimates[0];
        theta1=MA3_estimates[1];
	theta2=MA3_estimates[2];
        c2=MA3_estimates[3];
	free(MA3_estimates);
        /*printf("phi:%lf\ntheta:%lf\ncons:%lf\n",phi_cap,theta,c2);*/
        phi_cap_new=phi_cap;
        theta1_new=theta1;
	theta2_new=theta2;
        c2_new=c2;

	phi_cap_int=phi_cap;
        theta1_int=theta1;
        theta2_int=theta2;
        c2_int=c2;
	length_int=length;
	Duplicate(Y_t,Y_t_int,length);
	Duplicate(Y_t_cap_int,Y_t_cap_int,length);
	/*MA3_estimates_new=MA3_estimates;*/
        int counter=0;
	double Pred_estimates[4];	
	 do{
                phi_cap=phi_cap_new;
               	theta1=theta1_new;
		theta2=theta2_new;
               	c2=c2_new;
		/*MA3_estimates=MA3_estimates_new;*/
               	/*printf("phi_cap%lf\ntheta:%lf\nc:%lf\n",phi_cap,theta,c2);*/
		Pred_estimates[0]=phi_cap;
		Pred_estimates[1]=theta1;
		Pred_estimates[2]=theta2;
		Pred_estimates[3]=c2;

	        LRM_Pred(length,3,Y_X,Y_Y_cap,Pred_estimates);
		calcMDArrDiff(length,1,Y_Y,Y_Y_cap,Y_Y_error);

                /*for(int i=0;i<length;i++){
                        printf("Diff:%lf\n",Y_error[i]);
                }*/
                MD2D(length,Y_Y_error,E_t);
		
		calcArrLeadLag(Y_t,Y_t,length,0,2);
		calcArrLeadLag(Y_t_1,Y_t_1,length,0,2);
		calcArrLeadLag(E_t,E_t_1,length,1,1);
		calcArrLeadLag(E_t,E_t_2,length,2,0);
                length=length-2;
		
		D2MD(length,Y_t,Y_Y);
		Cbind(length,3,Y_t_1,E_t_1,E_t_2,Y_X);
		MA3_estimates_new=LRM(length,n,Y_X,Y_Y);
		
		phi_cap_new=MA3_estimates_new[0];
        	theta1_new=MA3_estimates_new[1];
        	theta2_new=MA3_estimates_new[2];
        	c2_new=MA3_estimates_new[3];
		free(MA3_estimates_new);
		counter=counter+1;
		printf("%lf\t%lf\t%lf\t%lf\n",phi_cap_new,theta1_new,theta2_new,c2_new);
}while(!(fabs(phi_cap-phi_cap_new)<0.01 && fabs(theta1-theta1_new) <0.01 && fabs(theta2-theta2_new) < 0.01 && counter>30));
	printf("%d",counter);
	
	if(counter>30){
		phi_cap_new=phi_cap_int;
		theta1_new=theta1_int;
		theta2_new=theta2_int;
		length=length_int;
		mape=MAPE(Y_t_int,Y_t_cap_int,length);
	};
	
	printf("Phi: %lf\tTheta1: %lf\tTheta2: %lf\tConstant: %lf\n",phi_cap_new,theta1_new,theta2_new,c2_new);
	mape=MAPE(Y_t,Y_t_cap,length);
	printf("MAPE: %lf\n",mape);
	Y_n=Y_t[length-1],E_n=E_t_1[length-1],E_n_1=E_t_2[length-2];
        for(int i=0;i<16;i++){
                forecast[i]=(phi_cap_new*Y_n)+(theta1_new*E_n)+(theta2_new*E_n_1)+c2_new;
                E_n=0;
		E_n_1=0;	
                Y_n=forecast[i];
        }
	forecast[16]=mape;
	return forecast;
}

double * AR2MA1(double arr[],int length){
        double phi1,phi2,phi3,c1,Y_n,Y_n_1,E_n;
        int n=3,length_int;
	double *forecast=malloc(sizeof(double)*18);
        double phi_cap1,phi_cap2,theta,c2;
        double phi_cap1_new,phi_cap2_new,theta_new,c2_new;
	double phi_cap1_int,phi_cap2_int,theta_int,c2_int,mape;
        double *AR3_estimates,*MA3_estimates,*MA3_estimates_new;

        double DY_t[length-3],DY_t_1[length-3],DY_t_2[length-3],DY_t_3[length-3],DY_t_cap[length-2],DY_error[length-2];
        double Y_t[length-5],Y_t_int[length-5],Y_t_1[length-5],Y_t_2[length-5],Y_t_cap[length-5],Y_t_cap_int[length-5],Y_error[length-5];



        calcArrLeadLag(arr,DY_t,length,0,3);
        calcArrLeadLag(arr,DY_t_1,length,1,2);
        calcArrLeadLag(arr,DY_t_2,length,2,1);
        calcArrLeadLag(arr,DY_t_3,length,3,0);


        length=length-3;

        double DY_X[length][3];
        double DY_Y[length][1];
        double DY_Y_cap[length][1];
        double DY_Y_error[length][1];
        double E_t[length],E_t_1[length-1];

        D2MD(length,DY_t,DY_Y);
        Cbind(length,n,DY_t_1,DY_t_2,DY_t_3,DY_X);

        AR3_estimates=LRM(length,n,DY_X,DY_Y);

        phi1=AR3_estimates[0];
        phi2=AR3_estimates[1];
        phi3=AR3_estimates[2];
        c1=AR3_estimates[3];
	free(AR3_estimates);
        printf("phi1:%lf\nphi2:%lf\nphi3:%lf\ncon:%lf\n",AR3_estimates[0],AR3_estimates[1],AR3_estimates[2],AR3_estimates[3]);
        LRM_Pred(length,n,DY_X,DY_Y_cap,AR3_estimates);
        calcMDArrDiff(length,1,DY_Y,DY_Y_cap,DY_Y_error);
        /*for(int i=0;i<20;i++){
                printf("D:%lf\n",DY_error[i]);
        }*/
		MD2D(length,DY_Y_error,E_t);
        
		calcArrLeadLag(DY_t,Y_t,length,0,1);
        calcArrLeadLag(DY_t_1,Y_t_1,length,0,1);
        calcArrLeadLag(DY_t_2,Y_t_2,length,0,1);
        calcArrLeadLag(E_t,E_t_1,length,1,0);

        length=length-1;

        double Y_X[length][n];
        double Y_Y[length][1];
        double Y_Y_cap[length][1];
        double Y_Y_error[length][1];

        D2MD(length,Y_t,Y_Y);
        Cbind(length,3,Y_t_1,Y_t_2,E_t_1,Y_X);
        MA3_estimates=LRM(length,n,Y_X,Y_Y);

        phi_cap1=MA3_estimates[0];
        phi_cap2=MA3_estimates[1];
        theta=MA3_estimates[2];
        c2=MA3_estimates[3];
	free(MA3_estimates);
        /*printf("phi:%lf\ntheta:%lf\ncons:%lf\n",phi_cap,theta,c2);*/
        phi_cap1_new=phi_cap1;
        phi_cap2_new=phi_cap2;
	theta_new=theta;
        c2_new=c2;
        
	phi_cap1_int=phi_cap1;
        phi_cap2_int=phi_cap2;
        theta_int=theta;
        c2_int=c2;
	length_int=length;
	Duplicate(Y_t,Y_t_int,length);
	Duplicate(Y_t_cap,Y_t_cap_int,length);
	/*MA3_estimates_new=MA3_estimates;*/
        int counter=0;
        double Pred_estimates[4];
         do{
                phi_cap1=phi_cap1_new;
                phi_cap2=phi_cap2_new;
                theta=theta_new;
                c2=c2_new;
                /*MA3_estimates=MA3_estimates_new;*/
                /*printf("phi_cap%lf\ntheta:%lf\nc:%lf\n",phi_cap,theta,c2);*/
                Pred_estimates[0]=phi_cap1;
                Pred_estimates[1]=phi_cap2;
                Pred_estimates[2]=theta;
                Pred_estimates[3]=c2;

                LRM_Pred(length,3,Y_X,Y_Y_cap,Pred_estimates);
                calcMDArrDiff(length,1,Y_Y,Y_Y_cap,Y_Y_error);

                /*for(int i=0;i<length;i++){
                        printf("Diff:%lf\n",Y_error[i]);
                }*/
                MD2D(length,Y_Y_error,E_t);

                calcArrLeadLag(Y_t,Y_t,length,0,1);
                calcArrLeadLag(Y_t_1,Y_t_1,length,0,1);
                calcArrLeadLag(Y_t_2,Y_t_2,length,0,1);
                calcArrLeadLag(E_t,E_t_1,length,1,0);
                length=length-1;

                D2MD(length,Y_t,Y_Y);
                Cbind(length,3,Y_t_1,Y_t_2,E_t_1,Y_X);
                MA3_estimates_new=LRM(length,n,Y_X,Y_Y);

                phi_cap1_new=MA3_estimates_new[0];
                phi_cap2_new=MA3_estimates_new[1];
                theta_new=MA3_estimates_new[2];
                c2_new=MA3_estimates_new[3];
		free(MA3_estimates_new);
                /*printf("%lf\t%lf\t%lf\t%lf\n",phi_cap1_new,phi_cap2_new,theta_new,c2_new);*/
		counter=counter+1;
}while(!(fabs(phi_cap1-phi_cap1_new)<0.01 && fabs(phi_cap2-phi_cap2_new) < 0.01 && fabs(theta-theta_new) < 0.01 && counter>30));
        printf("%d",counter);
	
	Y_n=Y_t[length-1],Y_n_1=Y_t[length-2],E_n=E_t_1[length-1];

	if(counter>30){
		phi_cap1_new=phi_cap1_int;
		phi_cap2_new=phi_cap1_int;
		theta_new=theta_int;
		c2_new=c2_int;
		length=length_int;
		mape=MAPE(Y_t_int,Y_t_cap_int,length);
	};	
	
	printf("Phi1: %lf\tPhi2: %lf\tTheta: %lf\tConstant: %lf\n",phi_cap1_new,phi_cap2_new,theta_new,c2_new);
	mape=MAPE(Y_t,Y_t_cap,length);
	printf("MAPE: %lf\n",mape);
        for(int i=0;i<16;i++){
                forecast[i]=(phi_cap1_new*Y_n)+(phi_cap2_new*Y_n_1)+(theta_new*E_n)+c2_new;
                E_n=0;
		Y_n_1=Y_n;	
                Y_n=forecast[i];
        }
	forecast[16]=mape;
        return forecast;

}

double *  AR2MA2(double arr[],int length){
        double phi1,phi2,phi3,c1,Y_n,Y_n_1,E_n,E_n_1;
        int n=4,length_int;
	double *forecast=malloc(sizeof(double)*18);
        double phi_cap1,phi_cap2,theta1,theta2,c2;
        double phi_cap1_new,phi_cap2_new,theta1_new,theta2_new,c2_new;
	double phi_cap1_int,phi_cap2_int,theta1_int,theta2_int,c2_int,mape;
        double *AR3_estimates,*MA3_estimates,*MA3_estimates_new;

        double DY_t[length-3],DY_t_1[length-3],DY_t_2[length-3],DY_t_3[length-3],DY_t_cap[length-2],DY_error[length-2];
        double Y_t[length-5],Y_t_int[length-5],Y_t_1[length-5],Y_t_2[length-5],Y_t_cap[length-5],Y_t_cap_int[length-5],Y_error[length-5];



        calcArrLeadLag(arr,DY_t,length,0,3);
        calcArrLeadLag(arr,DY_t_1,length,1,2);
        calcArrLeadLag(arr,DY_t_2,length,2,1);
        calcArrLeadLag(arr,DY_t_3,length,3,0);
	
        length=length-3;

        double DY_X[length][3];
        double DY_Y[length][1];
        double DY_Y_cap[length][1];
        double DY_Y_error[length][1];
        double E_t[length],E_t_1[length-2],E_t_2[length-2];
	D2MD(length,DY_t,DY_Y);
        Cbind(length,3,DY_t_1,DY_t_2,DY_t_3,DY_X);
        
	AR3_estimates=LRM(length,3,DY_X,DY_Y);

        phi1=AR3_estimates[0];
        phi2=AR3_estimates[1];
        phi3=AR3_estimates[2];
        c1=AR3_estimates[3];
	free(AR3_estimates);
        /*printf("%lf\t%lf\t%lf\t%lf\t%d\n",phi1,phi2,phi3,theta2,length);*/
	LRM_Pred(length,3,DY_X,DY_Y_cap,AR3_estimates);
        calcMDArrDiff(length,1,DY_Y,DY_Y_cap,DY_Y_error);
	
        MD2D(length,DY_Y_error,E_t);

        calcArrLeadLag(DY_t,Y_t,length,0,2);
        calcArrLeadLag(DY_t_1,Y_t_1,length,0,2);
        calcArrLeadLag(DY_t_2,Y_t_2,length,0,2);
        calcArrLeadLag(E_t,E_t_1,length,1,1);
	calcArrLeadLag(E_t,E_t_2,length,2,0);
	
	 /*for(int i=0;i<length-2;i++){
                        printf("%lf,%lf,%lf,%lf,%lf\n",Y_t[i],Y_t_1[i],Y_t_2[i],E_t_1[i],E_t_2[i]);
                }*/	
        length=length-2;

        double Y_X[length][n];
        double Y_Y[length][1];
        double Y_Y_cap[length][1];
        double Y_Y_error[length][1];

        D2MD(length,Y_t,Y_Y);
        Cbind4(length,4,Y_t_1,Y_t_2,E_t_1,E_t_2,Y_X);
        
	MA3_estimates=LRM(length,n,Y_X,Y_Y);

        phi_cap1=MA3_estimates[0];
        phi_cap2=MA3_estimates[1];
        theta1=MA3_estimates[2];
	theta2=MA3_estimates[3];
        c2=MA3_estimates[4];
	free(MA3_estimates);
	/*printf("%lf\t%lf\t%lf\t%lf\n",phi1,phi2,phi3,c1);*/
        /*printf("phi:%lf\ntheta:%lf\ncons:%lf\n",phi_cap,theta,c2);*/
        phi_cap1_new=phi_cap1;
        phi_cap2_new=phi_cap2;
	theta1_new=theta1;
	theta2_new=theta2;
        c2_new=c2;

	phi_cap1_int=phi_cap1;
        phi_cap2_int=phi_cap2;
        theta1_int=theta1;
        theta2_int=theta2;
        c2_int=c2;
	length_int=length;
	Duplicate(Y_t,Y_t_int,length);
	Duplicate(Y_t_cap,Y_t_cap_int,length);
	/*printf("%lf\t%lf\t%lf\t%lf\t%lf\t%d\n",phi_cap1_new,phi_cap2_new,theta1_new,theta2_new,c2_new,length);*/
        int counter=0;
        double Pred_estimates[5];
         do{
                phi_cap1=phi_cap1_new;
                phi_cap2=phi_cap2_new;
                theta1=theta1_new;
		theta2=theta2_new;
                c2=c2_new;
                /*MA3_estimates=MA3_estimates_new;*/
                /*printf("phi_cap%lf\ntheta:%lf\nc:%lf\n",phi_cap,theta,c2);*/
                Pred_estimates[0]=phi_cap1;
                Pred_estimates[1]=phi_cap2;
                Pred_estimates[2]=theta1;
		Pred_estimates[3]=theta2;
                Pred_estimates[4]=c2;

                LRM_Pred(length,4,Y_X,Y_Y_cap,Pred_estimates);
                calcMDArrDiff(length,1,Y_Y,Y_Y_cap,Y_Y_error);
		/*for(int i=0;i<length;i++){
                        printf("%lf\n",Y_Y_cap[i][0]);
                }*/
                
		MD2D(length,Y_Y_error,E_t);

                calcArrLeadLag(Y_t,Y_t,length,0,2);
                calcArrLeadLag(Y_t_1,Y_t_1,length,0,2);
                calcArrLeadLag(Y_t_2,Y_t_2,length,0,2);
                calcArrLeadLag(E_t,E_t_1,length,1,1);
		calcArrLeadLag(E_t,E_t_2,length,2,0);
                length=length-2;
		
		/*for(int i=0;i<length;i++){
                        printf("%lf\t%lf\t%lf\t%lf\t%lf\n",Y_t[i],Y_t_1[i],Y_t_2[i],E_t_1[i],E_t_2[i]);
                }*/
                D2MD(length,Y_t,Y_Y);
                Cbind4(length,4,Y_t_1,Y_t_2,E_t_1,E_t_2,Y_X);
                MA3_estimates_new=LRM(length,n,Y_X,Y_Y);

                phi_cap1_new=MA3_estimates_new[0];
                phi_cap2_new=MA3_estimates_new[1];
                theta1_new=MA3_estimates_new[2];
		theta2_new=MA3_estimates_new[3];
                c2_new=MA3_estimates_new[4];
		free(MA3_estimates_new);
		printf("%lf\t%lf\t%lf\t%lf\t%lf\t%d\n",phi_cap1_new,phi_cap2_new,theta1_new,theta2_new,c2_new,length);
                /*printf("%lf\t%lf\t%lf\n",fabs(phi_cap1-phi_cap1_new),fabs(phi_cap2-phi_cap2_new),fabs(theta1-theta1_new));*/
                counter=counter+1;
}while(!(fabs(phi_cap1-phi_cap1_new)<0.01 && fabs(phi_cap2-phi_cap2_new)<0.01 && fabs(theta1-theta1_new) <0.01 && fabs(theta2-theta2_new) < 0.01 && counter>30));

        /*printf("%d",counter);*/
	
	Y_n=Y_t[length-1],Y_n_1=Y_t[length-2],E_n=E_t_1[length-1],E_n_1=E_t_2[length-2];

	if(counter>30){
		phi_cap1_new=phi_cap1_int;
		phi_cap2_new=phi_cap2_int;
		theta1_new=theta1_int;
		theta2_new=theta2_int;
		c2_new=c2_int;
		length=length_int;
		mape=MAPE(Y_t_int,Y_t_cap_int,length);
	};
	
	printf("Phi1: %lf\tPhi2: %lf\tTheta1: %lf\tTheta2: %lf\tConstant: %lf\n",phi_cap1_new,phi_cap2_new,theta1_new,theta2_new,c2_new);
	mape=MAPE(Y_t,Y_t_cap,length);
	printf("MAPE: %lf\n",mape);
        for(int i=0;i<16;i++){
                forecast[i]=(phi_cap1_new*Y_n)+(phi_cap2_new*Y_n_1)+(theta1_new*E_n)+(theta2_new*E_n_1)+c2_new;
                E_n_1=0;
		E_n=0;	
		Y_n_1=Y_n;
                Y_n=forecast[i];
        }
	forecast[16]=mape;
        return forecast;
}

void Forecast_Func(double arr_forecast[],double arr_recov[],double arr_forecast_final[],int length,int d){
        double arr_forecast_final_tmp[length],X;
        Duplicate(arr_forecast,arr_forecast_final,length);
        for(int i=0;i<d;i++){
                X=arr_recov[d-i-1];
                for(int j=0;j<length;j++){
                        arr_forecast_final_tmp[j]=arr_forecast_final[j]+X;
                        X=arr_forecast_final_tmp[j];
                        arr_forecast_final[j]=arr_forecast_final_tmp[j];
                        }
        }
}


int main(int argc, char **argv)
{	
	if(argc < 2)
        {
                printf("Usage : command <<source_file.csv>>  <<destination_file.csv>>\n");
                return(0);
        }
	
	/** read the input and output file names from the command line  **/
        char *FILE_TO_READ=argv[1];
        char *FILE_TO_WRITE=argv[2];

	FILE *fp=fopen(FILE_TO_READ,"r+");
        FILE *fp2=fopen(FILE_TO_WRITE,"w+");
        
	if(!fp){
                printf("File does not exist\n");
                return(0);
        }

	char str[2000];
        //double dl[NUM_OBSERVATIONS];
	double arr[NUM_OBSERVATIONS];
	double actual[16];
        //double *ddl;
  	int arr_d[1374];
 	int counter=-1; 
	while(fgets(str,1500,fp)!= NULL){
		const char s[2] = ",";
        	char *token;
               	char *ch;
                char c;
		counter=counter+1;
		//double k[]={111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126};		
		fputs(str,fp2);
		int columnCount=0;
		/* get the first token */
                token = strtok(str, s);

		/* walk through other tokens */
                while( token != NULL )
                {
		        /* scan the read token into required format*/
                        if(columnCount >0 && columnCount <= NUM_OBSERVATIONS)
                        //if(columnCount >0)
                        {
                                //sscanf(token,"%lf",ddl);
				sscanf(token,"%lf",&arr[columnCount-1]);
                                //printf("%lf\n",ddl[0]);
				//printf("inside while %lf\n",arr[columnCount]);
                        }else if(columnCount >  NUM_OBSERVATIONS && columnCount <= (NUM_OBSERVATIONS+16)){
				sscanf(token,"%lf",&actual[columnCount-151]);
			}
                        token = strtok(NULL, s);
                        columnCount++;
                }
		
		double *Estimates;
		double *Forecast;
		double arr_forecast[18],arr_forecast_final[18],mape_os;

		int m,n,size,length,d;
		double EAF[3][3],threshold;
		int priority_base[3][3]={{1,3,6},{2,5,8},{4,7,9}};
		int priority[9]={99,99,99,99,99,99,99,99,99},model_index,points;
		

		//double arr[150]={8,4,2,4,3,2,4,2,3,3,3,1,3,1,2,2,4,5,4,1,6,0,4,4,6,2,4,1,2,3,5,9,4,4,3,2,5,4,4,3,5,9,2,5,2,2,4,5,2,7,4,1,4,3,10,7,6,2,7,10,5,3,1,5,2,5,7,7,4,6,6,3,2,4,3,6,3,3,6,5,4,5,3,3,4,9,8,2,1,10,4,4,4,6,6,3,5,5,1,6,3,9,6,6,1,4,5,4,3,3,6,4,4,4,8,4,4,4,5,2,3,8,3,8,4,4,4,4,5,4,4,4,4,10,10,13,2,7,8,4,6,9,6,8,8,7,1,5,4,7};
	
	/*double actual[16]={8,4,2,4,3,2,4,2,3,3,3,1,3,1,2,2};*/
	//printf("done!\n");
	length=sizeof(arr)/sizeof(arr[0]);
	//printf("done done\n");
	double inarr[length];
	double arr_recov[100];
	
	d=DFTest(arr,arr_recov,length);
	arr_d[counter]=d;
	
	//printf("\nD value:%d\n\n",d);
	Duplicate(arr,inarr,length);
	Drift(arr,inarr,length,d);
	length=sizeof(inarr)/sizeof(inarr[0]);
	
	/*for(int i=0;i<length;i++){
		printf("\nValue:%lf\n",inarr[i]);
	}*/
	
	m=length;
	size=length;
	threshold=1.96/(sqrt(length));
	//printf("\nNew Length: %d\n",length);	
	EAFMatrix(inarr,EAF,length);
	
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			if(EAF[i][j] < (-1*threshold) || EAF[i][j] > threshold){
				priority_base[i][j]=0;
			}
		/*printf("%d\t",priority_base[i][j]);*/
		}
		/*printf("\n");*/
	}
	

	int index=-1;
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			if(priority_base[i][j]!=0){
				index=index+1;
				priority[index]=priority_base[i][j];

				}
		}
	}
	
	model_index=priority[calcArrMin(priority,9)];
	/*model_index=9;*/
	
	if(model_index==1){
		printf("\n\nUse Average as Forecast\n");	
	}else if(model_index==2){
		Forecast=AR1(inarr,length);
		printf("\n\nUsing AR1 Model for forecasting\n");	
	}else if(model_index==3){
		Forecast=MA1(inarr,length);
		 printf("\n\nUsing MA1 Model for forecasting\n");
	}else if(model_index==4){
		Forecast=AR2(inarr,length);
		 printf("\n\nUsing AR2 Model for forecasting\n");
	}else if(model_index==5){
		Forecast=AR1MA1(inarr,length);
		 printf("\n\nUsing AR1MA1 Model for forecasting\n");
	}else if(model_index==6){
		Forecast=MA2(inarr,length);
		 printf("\n\nUsing MA2 Model for forecasting\n");
	}else if(model_index==7){
		Forecast=AR2MA1(inarr,length);
		 printf("\n\nUsing AR2MA1 Model for forecasting\n");
	}else if(model_index==8){
		Forecast=AR1MA2(inarr,length);
		 printf("\n\nUsing AR1MA2 Model for forecasting\n");
	}else if(model_index==9){
		Forecast=AR2MA2(inarr,length);
		 printf("\n\nUsing AR2MA2 Model for forecasting\n");	
	}
	
	if(model_index>1){
	for(int i=0;i<16;i++){
		/*printf("%lf\n",Forecast[i]);*/
		arr_forecast[i]=Forecast[i];
		}
		printf("\nMAPE: %lf\n",Forecast[16]);
		arr_forecast[16]=Forecast[16];
		
		/*mape_os=MAPE(actual,arr_forecast,16);
        	arr_forecast[17]=mape_os;
        	printf("\nMAPE_OS: %lf\n",mape_os);*/

	}
	
	if(d>0){
		Forecast_Func(arr_forecast,arr_recov,arr_forecast_final,16,d);
	}else{
		Duplicate(arr_forecast,arr_forecast_final,16);
	}


 	 mape_os=MAPE(actual,arr_forecast_final,16);
         arr_forecast[17]=mape_os;
         printf("\nMAPE_OS: %lf\n",mape_os);
	
	/** get the current position in the file. **/
        int len=ftell(fp);
        c = fgetc(fp); // moves offset by 1
        fseek(fp,len,SEEK_SET); // Resets the pointer to the place before the previous command
	fseek(fp2,ftell(fp2)-2,SEEK_SET);
	
	for(int i=0;i<TERMS_TO_PREDICT;i++)
        {
        	/** append , and the predicted value **/
                fputs(",",fp2);
               	//fprintf(fp2,"%lf",k);
                fprintf(fp2,"%lf",arr_forecast[i]);
        }

	/** append \n to mark end of line **/
        if(c!=EOF)
        {
        	fputs("\n",fp2);
        }


	if(model_index>1){
		free(Forecast);
	}
	
	}
	for(int i=0;i<counter;i++){
                printf("D: %d",arr_d[i]);
                }

	fclose(fp);
        fclose(fp2);
	printf("completed!\n");
	return(0);
}
