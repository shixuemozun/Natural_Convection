#include<iostream>
# include<cmath>
# include<cstdlib>
# include<iomanip>
# include<fstream>
# include<sstream>
# include<string>
using namespace std;
const int NX = 120;
const int NY = 120; 

double u[NX+1][NY+1],u1[NX+1][NY+1],v[NX+1][NY+1],v1[NX+1][NY+1],T[NX+1][NY+1],T1[NX+1][NY+1],p[NX+1][NY+1],p1[NX+1][NY+1],pg[NX+1][NY+1];
double U1[NX+1][NY+1];
double Pr,Gr;
double A,B,C;
double dx,dy,dt;
int i,j;
double a,b,c,d; 
double Lx,Ly;

int main()
{
	
	dx = 1;
	dy = 1;
    dt = 0.001;
    Pr = 1;
    Gr = 0.05;
    
    Lx = dx*double(NX);
	Ly = dy*double(NY);

//��ʼ�� 

  for(i=0;i<=NX;i++)
  {
  	
  	for(j=0;j<=NY;j++)
  	{
  		u[i][j] = 0;
  		u1[i][j] = 0; 
  		v[i][j] = 0;
  		v1[i][j] = 0;
  		U1[i][j] = 0;
		  
		T[i][j] = 0;
  		T1[i][j] = 0;
  		
  		T[0][j] = 0.5;
  		T[NX][j] = -0.5;
  		
  		T1[0][j] = 0.5;
  		T1[NX][j] = -0.5;
  		
  		p[i][j] = 0;
		p1[i][j] = 0;   
	  
	  }
  	
  	
  }

/*	
// x ����������
        
	  A =  Pr*((u[i+1][j]-2*u[i][j]+u[i-1][j])/(dx*dx)+(u[i][j+1]-2*u[i][j]+u[i][j-1])/(dy*dy));
	  u1[i][j] = u[i][j]+dt*(-1*u[i][j]*((u[i+1][j]-u[i-1][j])/(2*dx))-1*v[i][j]*((u[i][j+1]-u[i][j-1])/(2*dy))-((p[i+1][j]-p[i-1][j])/(2*dx))+A);
	
//y����������

 	B = Pr*((v[i+1][j]-2*v[i][j]+v[i-1][j])/(dx*dx)+(v[i][j+1]-2*v[i][j]+v[i][j-1])/(dy*dy))+Gr*Pr*Pr*T[i][j];	
	v1[i][j] = v[i][j]+dt*(-1*u[i][j]*((v[i+1][j]-v[i-1][j])/(2*dx))-1*v[i][j]*((v[i][j+1]-v[i][j-1])/(2*dy))-((p[i][j+1]-p[i][j-1])/(2*dy))+B);
	
	
//��������
   C = ((T[i+1][j]-2*T[i][j]+T[i-1][j])/(dx*dx))+((T[i][j+1]-2*T[i][j]+T[i][j-1])/(dy*dy)); 
   T1[i][j] = T[i][j]+dt*(-1*u[i][j]*((T[i+1][j]-T[i-1][j])/(2*dx))-1*v[i][j]*((T[i][j+1]-T[i][j-1])/(2*dy))+C);		
*/	
	
for(int t=0;t<2500000;t++)
{
cout << t <<endl;

	
//Ux

for(i=1;i<NX;i++)
{
	
	for(j=1;j<NY;j++)
	{ 
		
	  A =  Pr*((u[i+1][j]-2*u[i][j]+u[i-1][j])/(dx*dx)+(u[i][j+1]-2*u[i][j]+u[i][j-1])/(dy*dy));
	  u1[i][j] = u[i][j]+dt*(-1*u[i][j]*((u[i+1][j]-u[i-1][j])/(2*dx))-1*v[i][j]*((u[i][j+1]-u[i][j-1])/(2*dy))-((p[i+1][j]-p[i][j])/(dx))+A);	
		
	}
	

	}	
	
//Ux�߽�����
//���ұ߽� 
for(j=1;j<NY;j++) 
{
	u1[0][j] = 0;
	u1[NX][j] = 0;
	
}

//���±߽�

for(i=0;i<=NX;i++)
{
	u1[i][0] = 0;
	u1[i][NY] = 0;
	
	
}
 
//���� u ֵ

for(i=0;i<=NX;i++)
{
	for(j=0;j<=NY;j++)
	{
		u[i][j] = u1[i][j];
		
		
	}
	
	
 } 
 
 
 
 

//Uy

for(i=1;i<NX;i++)
{
	for(j=1;j<NY;j++)
	{
		
	B = Pr*((v[i+1][j]-2*v[i][j]+v[i-1][j])/(dx*dx)+(v[i][j+1]-2*v[i][j]+v[i][j-1])/(dy*dy))+Gr*Pr*Pr*T[i][j];	
	v1[i][j] = v[i][j]+dt*(-1*u[i][j]*((v[i+1][j]-v[i-1][j])/(2*dx))-1*v[i][j]*((v[i][j+1]-v[i][j-1])/(2*dy))-((p[i][j+1]-p[i][j])/(dy))+B);	
		
		
	}
	
	
}

//Uy�߽�
//���ұ߽�

for(j=1;j<NY;j++)
{
	
	v1[0][j] = 0;
	v1[NX][j] = 0;
	
 } 

//���±߽�

for(i=0;i<=NX;i++)
{
	v1[i][0] = 0;
	v1[i][NY] = 0;
	
	
 } 
 
 
 //���� Uy ֵ 
 for(i=0;i<=NX;i++)
{
	for(j=0;j<=NY;j++)
	{
		v[i][j] = v1[i][j];
		
		
	}
	
	
 } 
 
for(int l=0;l<2;l++)
{
  

//ѹ������ 

for(i=0;i<=NX;i++)
{
	
	for(j=0;j<=NY;j++)
	{
		
		pg[i][j] = p[i][j];
		
	 } 
	
 } 

for(i=1;i<NX;i++)
{
for(j=1;j<NY;j++)
{

a = 2*((dt)/(dx*dx)+(dt)/(dy*dy));
b = -1*(dt/(dx*dx));
c = -1*(dt/(dy*dy));
d = (u[i][j]-u[i-1][j])+(v[i][j]-v[i][j-1]);


if(isnan(d))
    {
        
        d = 0;       //������Сʱ����� nan; 
    }
 else
 {
 	d = d;
 	
   }  
    

p[i][j] = 1e-5*(-1/a)*(b*p[i+1][j]+b*p[i-1][j]+c*p[i][j+1]+c*p[i][j-1]+d); //ѹ�� p ��������ɢ������ͳ���1e-4�������Լ��Գ�����;

}
}

for(i=1;i<NX;i++)
{
for(j=1;j<NY;j++)
{

p1[i][j] = pg[i][j]+0.035*p[i][j]; //�ٳ���0.035���ɳ�����ϵ��; 
}
}

//�߽��ϵ�ѹ��

//���ұ߽�

for(j=1;j<NY;j++)
{
	p1[0][j] = p1[1][j];
	p1[NX][j] = p1[NX-1][j];
 } 

//���±߽�

for(i=0;i<=NX;i++)
{
	p1[i][0] = p1[i][1];
	p1[i][NY] = p1[i][NY-1];
	
 } 


//���� p ֵ

for(i=0;i<=NX;i++)
{
	
	for(j=0;j<=NY;j++)
	{
		
		p[i][j] = p1[i][j];
		
	 } 
	
 } 


//���ø��µ� p �������� 

for(i=1;i<NX;i++)
{
	
	for(j=1;j<NY;j++)
	{ 
		
	 // A =  Pr*((u[i+1][j]-2*u[i][j]+u[i-1][j])/(dx*dx)+(u[i][j+1]-2*u[i][j]+u[i][j-1])/(dy*dy));
	  u1[i][j] = u[i][j]-dt*(p[i+1][j]-p[i][j])/(dx);	
		
	}
	

	}	
	
//Ux�߽�����
//���ұ߽� 
for(j=1;j<NY;j++) 
{
	u1[0][j] = 0;
	u1[NX][j] = 0;
	
}

//���±߽�

for(i=0;i<=NX;i++)
{
	u1[i][0] = 0;
	u1[i][NY] = 0;
	
	
}
 
//����Ux

for(i=0;i<=NX;i++)
{
	for(j=0;j<=NY;j++)
	{
		u[i][j] = u1[i][j];
		
		
	}
	
	
 }  
 
 

//Uy

for(i=1;i<NX;i++)
{
	for(j=1;j<NY;j++)
	{
		
//	B = Pr*((v[i+1][j]-2*v[i][j]+v[i-1][j])/(dx*dx)+(v[i][j+1]-2*v[i][j]+v[i][j-1])/(dy*dy))+Gr*Pr*Pr*T[i][j];	
	v1[i][j] = v[i][j]-dt*(p[i][j+1]-p[i][j])/(dy);	
		
		
	}
	
	
}

//Uy�߽�
//���ұ߽�

for(j=1;j<NY;j++)
{
	
	v1[0][j] = 0;
	v1[NX][j] = 0;
	
 } 

//���±߽�

for(i=0;i<=NX;i++)
{
	v1[i][0] = 0;
	v1[i][NY] = 0;
	
	
 } 
//����Uy 
for(i=0;i<=NX;i++)
{
	for(j=0;j<=NY;j++)
	{
		v[i][j] = v1[i][j];
		
		
	}
	
	
 } 
}
//�¶ȳ�����

for(i=1;i<NX;i++)
{
	for(j=1;j<NY;j++)
	{
	C = ((T[i+1][j]-2*T[i][j]+T[i-1][j])/(dx*dx))+((T[i][j+1]-2*T[i][j]+T[i][j-1])/(dy*dy)); 
   T1[i][j] = T[i][j]+dt*(-1*u[i][j]*((T[i+1][j]-T[i-1][j])/(2*dx))-1*v[i][j]*((T[i][j+1]-T[i][j-1])/(2*dy))+C);	
		
		
		
	}
	

 } 

//�߽�����

//���ұ߽�

for(j=1;j<NY;j++)
{
	
	T1[0][j] = 0.5;
	T1[NX][j] = -0.5;
	
	
 } 


//���±߽�

for(i=0;i<=NX;i++)
{
	T1[i][0] = T1[i][1];
	T1[i][NY] = T1[i][NY-1];    //���±߽����	
	
 } 

//����Tֵ

for(i=0;i<=NX;i++)
{
	for(j=0;j<=NY;j++)
	{
		T[i][j] = T1[i][j];
			
	}
 } 


//����ʸ�����ٶȳ�

for(i=0;i<=NX;i++)
{
	for(j=0;j<=NY;j++)
	{
		
		U1[i][j] = sqrt(u[i][j]*u[i][j]+v[i][j]*v[i][j]); 
		
	}
	
	
	
	
 } 








}

ostringstream name;
	  name<<"cavity_"<<6<<".dat";
	  ofstream out(name.str().c_str());
	  out<< "Title= \"LBM Lid Driven Flow\"\n" << "VARIABLES=\"X\",\"Y\",\"u\",\"v\",\"p\",\"T\",\"U1\"\n" << "ZONE T=\"BOX\",I=" << NX+1 << ",J=" << NY+1 << ",F=POINT" << endl;
	  for(int j=0;j<=NY;j++)
	     for(int i=0;i<=NX;i++)
	     {
	     	
			 
			 out<<double(i*dx)/Lx<<" "<<double(j*dx)/Ly<<" "<<u[i][j]<<" "<<v[i][j]<<" "<<p[i][j]<<" "<<T[i][j]<<" "<<U1[i][j]<<endl;
		       
		
		 }








	return 0;
 } 
