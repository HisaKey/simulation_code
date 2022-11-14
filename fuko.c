#include<stdio.h>
#include<math.h>

double fx(double xv);
double fxv(double x, double xv, double yv);

double fy(double yv);
double fyv(double y, double yv, double xv);

#define g 9.80665      /*重力加速度*/
#define l 100000        /*糸の長さ*/
#define pi 3.141592
#define rad (pi/180.0)
#define omega 0.0000727 
#define m 30
#define T 10000
#define dt 0.1
#define t_max 86400.0

int main()
{
        FILE *fp1, *fp2, *fp3;
        int i,i_max; int a;
        double x,xv,x_new,xv_new,
               y,yv,y_new,yv_new,t;/*各種変数*/

        x=0.5; y=0.0;
        xv=0.0; yv=0;/*各種初期値*/
        i_max = t_max/dt;

        fp1 = fopen("1.dat","w") ;
        fp2 = fopen("2.dat","w") ;
        fp3 = fopen("3.dat","w") ;
              for ( i=0; i<=i_max;i++)
	for ( i=0; i<=100000;i++)
	
	{

                t = dt*i;

                x_new = x + fx(xv)*dt;
                xv_new = xv + fxv(x,xv,yv)*dt;

                 y_new = y + fy(yv)*dt;
                yv_new = yv + fyv(y,yv,xv)*dt;

                x=x_new; y=y_new; xv=xv_new; yv=yv_new;
				printf("%lf",x);
                fprintf(fp1,"%lf\t %lf \n",x,y);
                fprintf(fp2,"%lf\t %lf \n",t,x);
                fprintf(fp3,"%lf\t %lf \n",t,y);
				

	}
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
return 0 ;
}

double fx(double xv)
{
        return xv;
}
             
double fxv(double x,double xv,double yv)
{
        return -1*(g/l)*x+2*omega*sin(60*rad)*yv;
}

double fy(double yv)
{
        return yv;
}

double fyv(double y,double yv,double xv)
{
        return -1*(g/l)*y-2*omega*sin(60*rad)*xv;
}