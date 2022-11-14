def fx(v):
	return v

def fv(x,v,t):
	return -x**3 -2*v/t

if __name__ =="__main__":
	t_max =5.0
	dt = 0.001
	x=1
	v = -1/3000
	i_max = t_max = t_max/dt
	patn_w = "Emden3_test.txt"
	fp = open("Emden3_test.txt","w")
	for i in range(int(i_max)):
		t = dt*i+0.001
		v_new = v + fv(x,v,t)*dt;
		x_new = x + fx(v)*dt;
		x = x_new;   v = v_new;
		
		t_str = str(t)
		x_str = str(x)
		v_str = str(v)
		
		txv = t_str+"\t"+x_str+"\t"+v_str
		fp.write(txv+"\n")
	fp.close()