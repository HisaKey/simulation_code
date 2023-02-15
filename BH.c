#include<stdio.h>
#include<math.h>


double fv (double u);

#define dphai 0.01

int main()
{

	int i, i_max, phai_max;
	double u, u_new,
		   phai;
	
	
	u = 1.347; phai = 0.0;/*+の方程式を扱う場合は u==限りなく0に近いやつ　φ=0.0を代入-の方程式のときは+のやつを実行したら途中でnanが出るからその時のuφを代入*/
	phai_max =3; i_max = phai_max / dphai;
	
	for (i=0; i <=i_max; i++)
		{
		
		 phai =dphai * i;

		 u_new = u + fv(u) * dphai;
		 u = u_new;
		 
		 
		 printf("%f\t %f\t %f\t %f\n",phai, 6/u ,6*cos(phai+2.39)/u,6*sin(phai+2.39)/u);/*+のときは何も足さない。ーのときはφの初期値を足す*/
		 }
	return 0;
}

		double fv( double u)
		{
			return  -sqrt(1 - u*u + 2*u*u*u/6) ;/*コメントで触れている+-はここの符号のこと*/
		}