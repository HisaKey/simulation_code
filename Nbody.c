#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#define NMAX 16384

#define VMAX 1.25

void calc_force(int n, double m[], double x[][3], double a[][3], double eps2)//粒子間の万有引力を計算している関数
{
	double r = 0.0;
	double r3inv = 0.0;
	for (int i =0;i<n;i++){
		for(int k=0;k<3;k++){
			a[i][k] = 0.0;
		}
	}
	
	double rr[3] ={0.0, 0.0, 0.0};
	for (int i =0;i<n;i++){
		for(int j=i+1;j<n;j++){
			rr[0] = x[j][0]-x[i][0];//粒子iとjのx軸方向の距離
			rr[1] = x[j][1]-x[i][1];//y方向の距離
			rr[2] = x[j][2]-x[i][2];//z方向の距離
			r = sqrt(rr[0]*rr[0] + rr[1]*rr[1] + rr[2]*rr[2]);//粒子ij間の距離
			r3inv = fabs(1.0/(r*r*r));//rの三乗の逆数。すなわちr^(-3)
			for(int k=0;k<3;k++){
				a[i][k] -= -m[j]*rr[k]*r3inv;
				a[j][k] -= m[i]*rr[k]*r3inv;
			}
		}
	}
}

void leap_frog(
			   int n, double m[], double x[][3], double v[][3], double a[][3], double dt, double eps2)//時間発展(蛙飛び法)
{
 
	double v_half[n][3];
	for(int i=0;i<n;i++){
		for(int k=0;k<3;k++){
			v_half[i][k] = v[i][k] + 0.5*a[i][k]*dt;
			x[i][k] += v_half[i][k]*dt;
		}
	}

	
	calc_force(n, m, x, a, eps2);
		
	for(int i=0;i<n;i++){
		for(int k=0;k<3;k++){
			v[i][k] = v_half[i][k] + 0.5*a[i][k]*dt;
		}
	}
}

double calc_energy(int n, double m[], double x[][3], double v[][3],double eps2)//エネルギーを計算する関数
{
	/*運動エネルギー(K)の計算*/
	double K =0.0;
	for(int i=0; i<n; i++){
		K += m[i]*(v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2])/2.0;
	}

	/*ポテンシャルエネルギー(W)の計算*/
	double W = 0.0;
	double r2 = 0.0;
	double sigma = 0.0;
	for(int i=0; i<n-1;i++){
		for(int j=i+1; j<n;j++){
			r2 = (x[j][0]-x[i][0])*(x[j][0]-x[i][0]) + (x[j][1]-x[i][1])*(x[j][1]-x[i][1])+ (x[j][2]-x[i][2])*(x[j][2]-x[i][2]);
			W += m[i]*m[j]/sqrt(r2 + eps2*eps2);
		}
	}
	return(K+W);
}

/*
double particle_bunpu(void)
{
	//粒子分布を入れる
	double x, y, z;
	do{
		x = 2.0*drand48() - 1.0;
		y = 2.0*drand48() - 1.0;
		z = 2.0*drand48() - 1.0;
		r = sqrt(x*x+y*y+z*z);
	}while(r >= 1.0);
	return (x,y,z);
}
*************************************/;
double gaussian(void) //平均0,分散1のガウス分布を導出する関数。Box-Muller法を用いている
{
	double x, y, r2;
	double z;
	do{
		x = 2.0*drand48() - 1.0;
		y = 2.0*drand48() - 1.0;// x,y共に-1から１までの少数の乱数を出力
		r2 = x*x + y*y;
	}while(r2 >= 1.0 || r2 == 0.0);//x,y共にガウス分布の範囲内であるときにzを出力
	z = sqrt(-2.0*log(r2)/r2)*x;
	return (z);
}

void make_spherical(int n,double m[],double x[][3],double v[][3],double r_v,double eps2)//粒子の位置、速度分布を出すための関数
{
	int i =0;
	double xx =0;
	double yy =0;
	double zz =0;
	double r =0;
	double mm = 1.0/n;
	double W = 0.0;
	double r2 = 0.0;
	double sigma = 0.0;
	
	//粒子分布を作る
	while(i<n){
		xx = 2.0*drand48() - 1.0;
		yy = 2.0*drand48() - 1.0;
		zz = 2.0*drand48() - 1.0;
		r = sqrt(xx*xx+yy*yy+zz*zz);
		if (r <= 1.0 && r >= 0.25)//条件式を満たすxx,yy,zzを粒子分布として採用。条件式の中身を変えることで分布も変えることができる
		{
			x[i][0] = xx;
			x[i][1] = yy;
			x[i][2] = zz;
			m[i] = mm;
			i += 1;
		}
	}
		for(int i=0; i<n-1;i++){
			for(int j=i+1; j<n;j++){
				r2 = (x[j][0]-x[i][0])*(x[j][0]-x[i][0]) + (x[j][1]-x[i][1])*(x[j][1]-x[i][1])+ (x[j][2]-x[i][2])*(x[j][2]-x[i][2]);
				W -= m[i]*m[j]/sqrt(r2 + eps2*eps2);//速度分散を出すために重力エネルギーを出している
			}
		}
		sigma = sqrt(2.0*r_v*fabs(W)/3.0);//速度分散
	for(int i =0; i<n;i++){
		for(int k=0;k<3;k++){
			v[i][k] = sigma*gaussian();//速度分布を出力
			}
		}
}

int main(void)//メイン間数
{
	/* 粒子数 */
	int n;
	/* 粒子の物理量(質量、位置、速度、加速度) */
	static double m[NMAX], x[NMAX][3], v[NMAX][3], a[NMAX][3];
	/* 全エネルギー、 全エネルギー(初期値)、 ビリアル比*/
	double e, e_ini, r_v;
	/* 時刻、 時間刻み幅、 終了時刻, データの出力間隔 */ 
	double t, dt, t_end, t_out;
	/* ソフトニング長 */
	double eps2;
	t_out = 1.0;
//	n =4096;
//	r_v = 0.1;
	eps2 = 0.03125;
	dt = 0.0078125;
	t_end = 1.0;
	//パラメータの設定 粒子数とビリアル比を標準入力で表現できるようにしている
	printf("n= ");
	scanf("%d",&n);
	printf("n= %d\n", n);
	
	printf("r_v= ");
	scanf("%lf", &r_v);
	printf("r_v= %lf\n", r_v);
	//必要な変数の設定
	for (int i =0; i<n;i++){
		m[i] = 1.0/n;//質量の設定ここでは全ての粒子を等質量にしている
	}
	make_spherical(n, m, x, v,r_v,eps2);//初期分布を設定

	e_ini = calc_energy(n, m, x, v, eps2);//初期エネルギーを求めている
	//初期の各粒子の加速度を導出 
	calc_force(n, m, x, a, eps2);
	
	FILE *data,*gp;
	char *data_file;
	gp = popen("gnuplot","w");//gnuplotを開く
	while(t<t_end){
		data_file = "out.dat";//全粒子の位置を書き込むためのファイル
		data = fopen(data_file,"w");
		
		for (int i =0 ;i<n; i++){
			fprintf(data,"%f %f %f\n",x[i][0],x[i][1],x[i][2]);
		}//このfor文で、ある時間の全粒子の位置を全てファイルに書き込んでいる
		fprintf(data,"\n");
		fclose(data);
		//gnuplotの設定を変更、出力する範囲を[-1.5:1.5]*[-1.5:1.5]*[-1.5:1.5]に限定している
		fprintf(gp, "set xrange[-1.0:1.0]\n");
		fprintf(gp, "set yrange[-1.0:1.0]\n");
		fprintf(gp, "set zrange[-1.0:1.0]\n"); 
		fprintf(gp, "set size ratio 1\n");

		fprintf(gp, "plot 'out.dat'\n");//実際に出力
		fflush(gp);
		sleep(t_out);//スリープ
		leap_frog(n, m, x, v, a, dt, eps2);
		t += dt;
	}
	pclose(gp);
	
//	e = calc_energy(n, m, x, v, r_v);
//	printf("%e %e\n", dt, fabs((e - e_ini)/e_ini)); 数値計算のごさを見るためのもの。時間刻みはばを小さくするほど、誤差が小さくなる(=精度が良くなる)

	return(0);
}
