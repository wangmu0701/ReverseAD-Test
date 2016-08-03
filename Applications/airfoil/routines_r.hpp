#include <fstream>
#include <cmath>
#include <cstdlib>
#include <complex>


template <typename T>
int time_cell(T* x1, T* x2, T* x3, T* x4, T* q,
	      T& adt)
	      
{	      

//-----------------------compute area/timestep for an individual cell-----------------------------------------	      

  int i; 
  T ri, u, v, p, c;
  T dx[4];
  T dy[4];

	dx[0] = x2[0] - x1[0];
	dx[1] = x3[0] - x2[0];
	dx[2] = x4[0] - x3[0];
	dx[3] = x1[0] - x4[0];

	dy[0] = x2[1] - x2[1];
	dy[1] = x3[1] - x2[1];
	dy[2] = x4[1] - x3[1];
	dy[3] = x1[1] - x4[1];
	
	ri = 1.0/q[0];
        //ri = q[0];
	u = ri*q[1];
	v =  ri*q[2];
	p = gm1*(q[3] - 0.5*ri*(q[1]*q[1] + q[2]*q[2]));
	c = sqrt(gam*ri*p);
	//c = gam*ri*p;
	
	adt = 0;

	for (i = 0; i < 4; i++) 
	{
	  adt = adt+fabs(u*dy[i]-v*dx[i])+c*sqrt(dx[i]*dx[i]+dy[i]*dy[i]); 
	  //adt = adt+fabs(u*dy[i]-v*dx[i])+c*(dx[i]*dx[i]+dy[i]*dy[i]); 
	} 

	adt = adt/cfl;

return 0;

} //end


/*********************************************************************************************************/

template <typename T>
int flux_face(T* x1, T* x2, T* q1, T* q2, T adt1, 
	      T adt2, T* res1, T* res2)
	      
{	      
	      
//------------------------compute flux through an individual face-----------------------------------------	      


int n, r; 
T dx, dy, mu, r1i, u1, v1, p1, vol1, r2i, u2, v2, p2, vol2;
T f[4];

	dx = x1[0] - x2[0];
	dy = x1[1] - x2[1];
	mu = 0.5*(adt1 + adt2)*eps;

	r1i = 1.0/q1[0];
	//r1i = q1[0];
	u1 = r1i*q1[1];
	v1 = r1i*q1[2];
	p1 = gm1*(q1[3] - 0.5*r1i*(q1[1]*q1[1] + q1[2]*q1[2]));
	vol1 = u1*dy - v1*dx;
	
	r2i = 1.0/q2[0];
        //r2i = q2[0];
	u2 = r2i*q2[1];
	v2 = r2i*q2[2];
	p2 = gm1*(q2[3] - 0.5*r2i*(q2[2]*q2[2] + q2[3]*q2[3]));
	vol2 = u2*dy - v2*dx;
	
	f[0] = 0.5*(vol1*q1[0]         + vol2*q2[0]         );
	f[1] = 0.5*(vol1*q1[1] + p1*dy + vol2*q2[1] + p2*dy );
	f[2] = 0.5*(vol1*q1[2] - p1*dx + vol2*q2[2] - p2*dx );
	f[3] = 0.5*(vol1*(q1[3] + p1)  + vol2*(q2[3] + p2)  );


	for (n = 0; n < 4; n++)
	{
		f[n] = f[n] + mu*(q1[n]-q2[n]);
		res1[n] = res1[n] + f[n];
		res2[n] = res2[n] - f[n];
	} 


return 0;

} //end


/*********************************************************************************************************/

template <typename T>
int flux_wall(T* x1, T* x2, T* q, T* res) 

{

//--------------------compute momentum from an individual wall face------------------------------------------
	       

int n; 
T dx, dy, ri, u,v, p;

	dx = x1[0] - x2[0];
	dy = x1[1] - x2[1];

	ri = 1.0/q[0];
	//ri = q[0];
	u = ri*q[1];
	v = ri*q[2];
	p = gm1*(q[3] - 0.5*ri*(q[1]*q[1] + q[2]*q[2]));
	
	res[1] = res[1] - p*dy;
	res[2] = res[2] + p*dx;
		  
        return 0;
       
} //end


/********************************************************************************************************/

template <typename T>
int lift_wall(T* x1, T* x2, T* q, T &lift) {

T dx, ri, u, v, p; // COMPUTE MOMENTUM FLUX FROM AN INDIVIDUAL WALL FACE

	dx = x1[0] - x2[0];
	
	ri = 1.0/q[0];
	//ri = q[0];
	u = ri*q[1];
	v = ri*q[2];
	p = gm1*(q[3] - 0.5*ri*(q[1]*q[1] + q[2]*q[2]));

	lift = lift + p*dx;
	
	return 0;
       
} //end

