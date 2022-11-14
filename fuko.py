import numpy as np
import math
def fx(xv):
	return xv

def fxv(x,xv,yv):
	return -1*(g/l)*x + 2*omega*np.sin(math.radians(60))*yv;

def fy(yv):
	return yv

def fyv(y,yv,xv):
	return -1*(g/l)*y -2*omega*np.sin(math.radians(60))*xv

g = 9.80655 #地表での重力加速度
l = 100000 #振り子の紐の長さ
omega = 0.0000727 #地球の角速度
dt = 0.1 #時間刻みはば
t_max = 86400.0 #計算時間の上限値(1日の秒数に対応)

if __name__ =="__main__":

	i_max = t_max/dt
	x = 0.5
	y = 0.0
	xv = 0.0
	yv = 0.0

	fp1 = open("x-y_fuko.txt","w")
	fp2 = open("t-x_fuko.txt","w")
	fp3 = open("t-y_fuko.txt","w")
	
	for i in range(int(i_max)):
		t = dt*i
		
		x_new = x + fx(xv)*dt
		xv_new = xv + fxv(x,xv,yv) * dt
		
		y_new = y + fy(yv)*dt
		yv_new = yv + fyv(y,yv,xv) * dt
		
		x = x_new
		y = y_new
		xv = xv_new
		yv = yv_new
		
		xy = str(x)+"\t"+str(y)
		tx = str(t)+"\t"+str(x)
		ty = str(t)+"\t"+str(y)
		
		fp1.write(xy+"\n")
		fp2.write(tx+"\n")
		fp3.write(ty+"\n")
	fp1.close()
	fp2.close()
	fp3.close()