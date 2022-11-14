#include<stdio.h>
#include<math.h>
  double fx(double v);
  double fv(double x, double v, double t);
  

#define t_max 5.0 /*時間の上限*/
#define dt 0.001/*時間の刻みはば*/

int main()
{
	int i, i_max,m;
	double x,v,t,x_new,v_new,n;
	FILE *fp;
	/*パラメータを表示*/
	/*scanf("%lf", &n);*/
	printf("#t\t\t x\t\tv \t\tE#\n");
	
	x=1; v=-1/3000;  i_max = t_max/dt;
			fp = fopen("Emden3.txt","w");
	for ( i=0; i<i_max; i++)
	{
		t = dt*i+0.001;
		v_new = v + fv(x,v,t)*dt;
		x_new = x + fx(v)*dt;
		x = x_new;   v = v_new;
		fprintf (fp,"%f\t %lf\t %lf \n",t,x,v);
	}
	fclose(fp);
return 0;
}
	
	double fx (double v)
	{
		return v;
	}
	
	double fv(double x, double v, double t)
	{
		return -pow(x,3)-2*v/t;
	}
