#include<iostream>
#include<math.h>
#include<stdlib.h>
#include<time.h>
//#include<algorithm>
#include<stdio.h>
#include<fstream>
#include <sys/stat.h>
#include "string.h"
#include <vector>
#include "mem_alloc.h"

using namespace std;

#define N 				248 // 
#define MAX_ERR 	1.0e-10
#define EPS 			1.0e-10
#define ERR_UV 		1.0e-9

int PR_ITE=20000, fixx=int(N/2),fixy=10;
double L = 124.0, RHO = 0.2, DELT=0.025; //, UTOP = 1.0, uwall=0.0,vwall=0.0;
int 		ITE 				= 10000; //4000.0 / DELT;
double 	fix_value 	= 0.0 ;
double 	temperature = 1.5 ;
double 	Sigma 			= 1.0 ;
double 	AtomMass 		= 1.0 ;
double	CV					= 1.5 ;
double  timeNow			= 0.0 ;

double packingF,b_enskog,g_enskog;
double **u,**u_new,**u_star;
double **v,**v_new,**v_star;
double **Pr,**sigma,h,**vorticity,**tav1,**tav2;
double *vv,*x,*r,*r0,*p,*t,*s,*Ax,*midU,*midV;
double *strLW,*strTW,*strRW,*strBW;


double delo,deln,alfa,beta,infi_norm,rho,rhon,omega,shearF;
int nsq,UV_norm,Boupoints;
double ***bou_conditions,**LeftBou,**RightBou,**TopBou,**BottomBou;

double **T, Kcorrection, Vcorrection;


inline double VISCO(double cellTemp)
{
	// T=1.5, rho=0.05, Mu=0.222, rho=0.1 Mu=0.233, rho=0.2 Mu=0.274, rho=0.3 Mu=0.348
	
	return 0.274; //viscosity ;
}
inline double conductivity(double cellTemp)
{
	// T=1.5, rho=0.05, k=0.866, rho=0.1 k=0.942, rho=0.2 k=1.162, rho=0.3 k=1.508
	
	return 1.162;//conductK ;
}

inline bool exists_file (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}
void apply_boundary(int m)
{
	switch ( m )
	{
		case 1 :
			// velocities on walls 
		  for(int i=1; i<N+1; i++)
		  {
		     u_star[0][i] = LeftBou[i-1][0];
		     u_star[N][i] = RightBou[i-1][0];
		     v_star[i][0] = BottomBou[i-1][1];
		     v_star[i][N] = TopBou[i-1][1];
		  }
		  // update velocities in ghost cells
		  for(int i = 1; i < N; i++)
		  {
		     v_star[0][i] = LeftBou[i-1][1] * 2.0 - v_star[1][i];
		     v_star[N+1][i] = RightBou[i-1][1] * 2.0 - v_star[N][i];
		     u_star[i][0] = BottomBou[i-1][0] * 2.0 - u_star[i][1];
		     u_star[i][N+1] = TopBou[i-1][0] * 2.0 - u_star[i][N];
		  }
			break;
			
		case 2 :
			for(int i=1; i<N+1; i++)
		  {
		     u[0][i] = u_new[0][i] = LeftBou[i-1][0];
		     u[N][i] = u_new[N][i] = RightBou[i-1][0];
		     v[i][0] = v_new[i][0] = BottomBou[i-1][1];
		     v[i][N] = v_new[i][N] = TopBou[i-1][1];
		  }
		  for(int i = 1; i < N; i++)
		  {
		     u[i][0]=u_new[i][0]; u[i][N+1]=u_new[i][N+1];
		     v[0][i]=v_new[0][i]; v[N+1][i]=v_new[N+1][i];
		     
		     v_new[0][i] = LeftBou[i-1][1] * 2.0 - v_new[1][i];
		     v_new[N+1][i] = RightBou[i-1][1] * 2.0 - v_new[N][i];
		     u_new[i][0] = BottomBou[i-1][0] * 2.0 - u_new[i][1];
		     u_new[i][N+1] = TopBou[i-1][0] * 2.0 - u_new[i][N];
		  }
			break;
			
		case 3 :
			for(int i = 1; i < N+1; i++){
				T[i][0] = -T[i][1] + 2. * BottomBou[i-1][2];
				T[i][N+1] = -T[i][N] + 2. * TopBou[i-1][2];
			}
			
			for(int i = 1; i< N+1; i++){
				T[0][i] = -T[1][i] + 2. * LeftBou[i-1][2];
				T[N+1][i] = -T[N][i] + 2. * RightBou[i-1][2];
			}
			
			T[0][0] = T[0][1];		T[0][N+1] = T[0][N];
			T[N+1][0] = T[N][0];	T[N+1][N+1] = T[N+1][N];
			break;
	}
}
void interpolate_boundary(int side, double& variable, int varPoint, double location)
{
	double delL,delR,delta;
	if( location < bou_conditions[side][0][0] ) { 
			delL = bou_conditions[side][0][0] - location;
			delR = bou_conditions[side][1][0] - location;
			variable = (bou_conditions[side][0][varPoint]/delL - bou_conditions[side][1][varPoint]/delR)/(1.0/delL - 1.0/delR);
		} else if( location > bou_conditions[side][Boupoints-1][0] ) {
			delL = bou_conditions[side][Boupoints-1][0] - location;
			delR = bou_conditions[side][Boupoints-2][0] - location;
			variable = (bou_conditions[side][Boupoints-1][varPoint]/delL - bou_conditions[side][Boupoints-2][varPoint]/delR)/(1.0/delL - 1.0/delR);
		} else {
			int j = 0;
			while( bou_conditions[side][j][0] < location ) { j++; if( j == Boupoints ) break;  }
			delL = (location - bou_conditions[side][j-1][0]);
			delR = (bou_conditions[side][j][0] - location);
			delta=bou_conditions[side][j][0] - bou_conditions[side][j-1][0];
			variable = (delL * bou_conditions[side][j][varPoint] + delR * bou_conditions[side][j-1][varPoint])/delta;
		}
		//if( side == 2 ) cout << "\n location " << location-250 << "  variable " << variable << "  varPoint " << varPoint;
}
double left_wall_exp(double rLoc)
{
	//return 0.2*0.0586*exp(-0.0218*rLoc)*(1.0 - exp(-5.0586*rLoc))*(1.0+3.914*exp(-0.165*rLoc)); // rho=0.2, U=0.2, LW exptl
	return 0.5*0.08*exp(-0.03878*rLoc)*(1.0 - exp(-3.866*rLoc))*(1.0+3.557*exp(-0.265*rLoc)); // rho=0.2, U=0.5, LW exptl
}
double right_wall_exp(double rLoc)
{
	//return -0.2*0.0821*exp(-0.0218*rLoc)*(1.0 - exp(-4.1319*rLoc))*(1.0+2.434*exp(-0.18*rLoc)); // rho=0.2, U=0.2, RW exptl
	return -0.5*0.1388*exp(-0.042*rLoc)*(1.0 - exp(-3.06*rLoc))*(1.0+1.137*exp(-0.276*rLoc)); // rho=0.2, U=0.5, RW exptl
}
double top_plate_exp(double rLoc)
{
	//if( rLoc < L/2.0 ) return 1.0017*0.2*0.98*(1.0-exp(-0.682*pow(rLoc,0.472)))*(1-exp(-3.859*rLoc)); // rho=0.2, U=0.2, exptl
	//else return 0.9982*0.2*0.978*(1.0-exp(-0.58*pow((L-rLoc),0.563)))*(1-exp(-5.288*(L-rLoc))); 
	
	if( rLoc < L/2.0 ) return 1.0006*0.5*0.969*(1.0-exp(-0.672*pow(rLoc,0.488)))*(1-exp(-4.071*rLoc)); // rho=0.2, U=0.5, exptl
	else return 0.9993*0.5*0.965*(1.0-exp(-0.593*pow((L-rLoc),0.6)))*(1-exp(-4.173*(L-rLoc))); 
}
void calculate_boundary()
{
	ofstream Output1, Output2, Output3, Output4;
	Output1.open("recordLeft.dat"); Output2.open("recordBottom.dat");
	Output3.open("recordRight.dat"); Output4.open("recordTop.dat");
	
	double location1, location2;
	for(int i=0; i < N; i++){ // left wall
		location1 = i * h + 0.5 * h; 
		//interpolate_boundary( 0, LeftBou[i][0], 1, location1);
		//interpolate_boundary( 0, LeftBou[i][2], 3, location1);
		LeftBou[i][0] = 0.0;
		LeftBou[i][2] = temperature;
		location2 = i * h +  h;
		//interpolate_boundary( 0, LeftBou[i][1], 2, location2);
		LeftBou[i][1] = left_wall_exp(L-location2);
		
		Output1 << location1 << "\t" << LeftBou[i][0] << "\t" << location2 << "\t" << LeftBou[i][1] << endl;
	}
	Output1.close();
	
	for(int i=0; i < N; i++){ // bottom
		location1 = i * h + 0.5 * h;
		//interpolate_boundary( 1, BottomBou[i][1], 2, location1);
		//interpolate_boundary( 1, BottomBou[i][2], 3, location1);
		BottomBou[i][1] = 0.0;
		BottomBou[i][2] = temperature;
		location2 = i * h +  h;
		//interpolate_boundary( 1, BottomBou[i][0], 1, location2);
		BottomBou[i][0] = 0.0;
		Output2 << location2 << "\t" << BottomBou[i][0] << "\t" << location1 << "\t" << BottomBou[i][1] << endl;
	}
	Output2.close();
	
	for(int i=0; i < N; i++){ // right wall
		location1 = i * h + 0.5 * h;
		//interpolate_boundary( 2, RightBou[i][0], 1, location1);
		//interpolate_boundary( 2, RightBou[i][2], 3, location1);
		RightBou[i][0] = 0.0;
		RightBou[i][2] = temperature;
		location2 = i * h +  h;
		//interpolate_boundary( 2, RightBou[i][1], 2, location2);		
		RightBou[i][1] = right_wall_exp(L-location2);
		
		Output3 << location1 << "\t" << RightBou[i][0] << "\t" << location2 << "\t" << RightBou[i][1] << endl;
	}
	Output3.close();
	
	for(int i=0; i < N; i++){ // moving plate
		location1 = i * h + 0.5 * h;
		//interpolate_boundary( 3, TopBou[i][1], 2, location1);
		//interpolate_boundary( 3, TopBou[i][2], 3, location1);
		TopBou[i][1] = 0.0;
		TopBou[i][2] = temperature;
		location2 = i * h +  h;
		//interpolate_boundary( 3, TopBou[i][0], 1, location2);
		TopBou[i][0] = top_plate_exp(location2);
		
		Output4 << location2 << "\t" << TopBou[i][0] << "\t" << location1 << "\t" << TopBou[i][1] << endl;
	}
	Output4.close();
}
void initialization()
{
		h = L/(double)N;
    nsq = N*N;
    
    allocate(u, N+1, N+2);
		allocate(u_new, N+1, N+2);
		allocate(u_star, N+1, N+2);
		allocate(v, N+2, N+1);
		allocate(v_new, N+2, N+1);
		allocate(v_star, N+2, N+1);
		allocate(Pr, N+2, N+2);
		allocate(sigma, N, N);
		allocate(vorticity, N, N);
		allocate(tav1, N+1, N+1);
		allocate(tav2, N+1, N+1);
		allocate(vv, N*N); allocate(x, N*N); allocate(r, N*N); allocate(r0, N*N); allocate(p, N*N);
		allocate(t, N*N); allocate(s, N*N); allocate(Ax, N*N);
		allocate(midU, N+1); allocate(midV, N+1);
		allocate(strLW, N); allocate(strTW, N); allocate(strRW, N); allocate(strBW, N);
		allocate(T, N+2, N+2);
		allocate(LeftBou,N,3); allocate(RightBou,N,3); allocate(TopBou,N,3); allocate(BottomBou,N,3); 
		
    for(int i=10; i<N-10; i++)
    {
        for(int j=10; j<N-10; j++)
        {
            u[j][i]=u_new[j][i]=u_star[j][i]= 0.01*sin(j*i) ;
            v[i][j]=v_new[i][j]=v_star[i][j]= 0.01*sin(j*i) ;
        }
        for(int j=0; j<N+2; j++)
        {
            Pr[i][j]=0.0;
            T[i][j] = temperature;
        }
    }
    
		/*char *filename = new char[20];
		strcpy(filename,"bou_conditionMod");
		if ( ! exists_file(filename) ) { cout << "Error. boundary conditions file not existent. \n"; exit(0); }
		delete [] filename;
		
		ifstream readBou;
		readBou.open("bou_conditionMod");
		
		int variables;
		readBou >> variables >> Boupoints;
		Allocate(bou_conditions,4,Boupoints,variables); 
		
		// left bottom right top
		for(int n = 0; n < 4; n++)
			for(int i = 0; i < Boupoints ; i++ )
				for(int j = 0; j < variables; j++ )
					readBou >> bou_conditions[n][i][j];
					
		readBou.close(); */
		calculate_boundary();
		
		apply_boundary(1);
		apply_boundary(2);
    for(int i=0; i<N+2; i++)
      for(int j=0; j<N+2; j++)
      {
         Pr[i][j] = 0.0;
         T[i][j] = temperature;
      }
    
    // CG
    for(int i=0; i<nsq; i++)
    {
        vv[i]=p[i]=r[i]=r0[i]=s[i]=t[i]=x[i]=EPS; 
    }
    shearF = 0.;
		
		// packing fraction 
		packingF = M_PI * RHO * Sigma * Sigma * Sigma / 6. ;
		// b = 2 * PI * Sigma^3 / 3
		b_enskog = 2.094395102 * Sigma * Sigma * Sigma;
	
		// g_enskog = ( 1 - PF/2 ) / ( 1 - PF )^3
		g_enskog = ( 1 - packingF/2. ) / pow( 1 - packingF , 3. );
		
		Kcorrection = 0.99 + 0.1597*RHO - 0.7464*RHO*RHO+1.2115*pow(RHO,3.0)-0.5583*pow(RHO,4.0);
		
		Vcorrection = 1.02 + 10.61*(RHO-0.495)*(RHO-0.495)*(RHO-0.495);
		
		/*for(int i = N-10; i < N ; i++ ){
				for(int j = 0; j < variables-1; j++ )
					cout << TopBou[i][j] << "\t";
				cout << endl;
		}*/
}

void provisional_velocity()
{
    for(int i=1; i<N; i++)
    {
        for(int j=1; j<N+1; j++)
        {
            double vconv = (v[i][j]+v[i+1][j]+v[i+1][j-1]+v[i][j-1])/4.0;
            double diffusion, u_x,u_y, avgMuUP,avgMuDN, tempAvg;
            // u * du/dx
            if(u[i][j]>0.) u_x = u[i][j]*(u[i][j]-u[i-1][j])/h;
            else u_x = u[i][j]*(u[i+1][j]-u[i][j])/h;
						
						// v * du/dy
            if(vconv>0.0) u_y = vconv*(u[i][j]-u[i][j-1])/h;
            else u_y = vconv*(u[i][j+1]-u[i][j])/h;
            
            // (2 / RHO) * d( mu * du/dx )/dx
            diffusion = 2 * (VISCO(T[i+1][j])*(u[i+1][j]-u[i][j])-VISCO(T[i][j])*(u[i][j]-u[i-1][j]));
            // (1 / RHO) * d( mu * du/dy )/dy
            tempAvg = ( T[i+1][j+1]+T[i+1][j]+T[i][j]+T[i][j+1] )/4.;
            //avgMuUP = ( VISCO(i+1,j+1)+VISCO(i+1,j)+VISCO(i,j)+VISCO(i,j+1) )/4.;
            avgMuUP = VISCO(tempAvg);
            diffusion += avgMuUP * (u[i][j+1]-u[i][j]) ;
            
            tempAvg = ( T[i+1][j-1]+T[i+1][j]+T[i][j]+T[i][j-1] )/4.;
            //avgMuDN = ( VISCO(i+1,j-1)+VISCO(i+1,j)+VISCO(i,j)+VISCO(i,j-1) )/4.;
            avgMuDN = VISCO(tempAvg);
            diffusion += -1.*avgMuDN * (u[i][j]-u[i][j-1]) ;
            // (1 / RHO) * d( mu * dv/dx )/dy
            diffusion +=  avgMuUP * (v[i+1][j]-v[i][j]) ;
            diffusion +=  -1.*avgMuDN * (v[i+1][j-1]-v[i][j-1]) ;
            // devide by h*h*RHO
            diffusion /= h * h * RHO;
            
            u_star[i][j] = u[i][j] + DELT*( -u_x -u_y + diffusion );
        }
    }
    for(int i=1; i<N+1; i++)
    {
        for(int j=1; j<N; j++)
        {
            double uconv = (u[i][j]+u[i][j+1]+u[i-1][j+1]+u[i-1][j])/4.0;
            double diffusion, v_x,v_y, avgMuUP, avgMuDN, tempAvg;
						
						// u * dv/dx
            if(uconv>0.0)  v_x=uconv*(v[i][j]-v[i-1][j])/h;
            else v_x = uconv*(v[i+1][j]-v[i][j])/h;
						
						// v * dv/dy
            if(v[i][j]>0.0) v_y = v[i][j]*(v[i][j]-v[i][j-1])/h;
            else v_y = v[i][j]*(v[i][j+1]-v[i][j])/h;
						
						// (2 / RHO ) * d( Mu * dv/dy )/dy
						diffusion = 2 * ( VISCO(T[i][j+1])*(v[i][j+1]-v[i][j]) - VISCO(T[i][j])*(v[i][j]-v[i][j-1]) );
						// (1 / RHO) * d( mu * dv/dx )/dx
						tempAvg = ( T[i+1][j+1]+T[i+1][j]+T[i][j]+T[i][j+1] )/4.;
						//avgMuUP = ( VISCO(i+1,j+1)+VISCO(i+1,j)+VISCO(i,j)+VISCO(i,j+1) )/4.;
						avgMuUP = VISCO(tempAvg);
						diffusion += avgMuUP * (v[i+1][j]-v[i][j]) ;
						
						tempAvg = ( T[i-1][j+1]+T[i-1][j]+T[i][j]+T[i][j+1] )/4.;
						//avgMuDN = ( VISCO(i-1,j+1)+VISCO(i-1,j)+VISCO(i,j)+VISCO(i,j+1) )/4.;
						avgMuDN = VISCO(tempAvg);
						diffusion += -1. * avgMuDN * (v[i][j]-v[i-1][j]) ;
						// (1 / RHO) * d( mu * du/dy )/dx
						diffusion += avgMuUP * (u[i][j+1]-u[i][j]);
						diffusion += -1. * avgMuDN * (u[i-1][j+1]-u[i-1][j]);
						// devide by h*h*RHO
            diffusion /= h * h * RHO;
						
            v_star[i][j] = v[i][j] + DELT*( -v_x - v_y + diffusion );
        }
    }
    
    apply_boundary(1);

}

void residual()
{
    if(r[0]<0.0)infi_norm=-1.0*r[0];
    else infi_norm=r[0];
    for(int i=1; i<nsq; i++)
    {
        if(r[i]<0.0){if(-1.0*r[i]>infi_norm) infi_norm=-1.0*r[i];}
        else if(r[i]>infi_norm) infi_norm=r[i];
    }
}
double b(int m)
{
    int yc=int(m/N)+1, xc = m-(yc-1)*N+1;
    double temp = h*RHO*(u_star[xc][yc]-u_star[xc-1][yc]+v_star[xc][yc]-v_star[xc][yc-1])/DELT;
    return temp;
}

void multiA(double* vecOut, double* vecIn)
{
    double termx,termy;
    int gridLoc = fixx-1+(fixy-1)*N;
    
    for(int i=0; i<nsq; i++)
    {
    		if( i == gridLoc ) continue;
        int yc=int((double)(i)/(double)N)+1, xc= i+1-(yc-1)*N;

        if(xc==1) termx = -vecIn[(yc-1)*N+xc-1]+vecIn[(yc-1)*N+xc];
        else if(xc==N) termx = -vecIn[(yc-1)*N+xc-1] + vecIn[(yc-1)*N+xc-2];
        else termx = vecIn[(yc-1)*N+xc-2] - 2.0*vecIn[(yc-1)*N+xc-1] + vecIn[(yc-1)*N+xc];

        if(yc==1) termy = -vecIn[(yc-1)*N+xc-1] + vecIn[(yc)*N+xc-1];
        else if(yc==N) termy = -vecIn[(yc-1)*N+xc-1] + vecIn[(yc-2)*N+xc-1];
        else termy = vecIn[(yc)*N+xc-1] - 2.0*vecIn[(yc-1)*N+xc-1] + vecIn[(yc-2)*N+xc-1];

        vecOut[i] = termx+termy;
    }
    
}
double multi_matrix(double* vec1, double* vec2)
{
  double temp=0.0;
  int gridLoc = fixx-1+(fixy-1)*N;
  
  for(int i=0; i<nsq; i++){ 
  	if( i == gridLoc ) continue;
   	temp += vec1[i] * vec2[i]; 
  }
  return temp;
}
void pressure_iteration()
{
    multiA(Ax,x);
    int gridLoc = fixx-1+(fixy-1)*N;
    
    for(int i=0; i<nsq; i++)
    {
    		//if( i == gridLoc ) continue;
        r0[i] = r[i] = b(i) - Ax[i];
        vv[i] = p[i] = s[i] = t[i] = 0.0;
    } 
    r0[gridLoc] = r[gridLoc] = 0.0;
    p[gridLoc] = s[gridLoc] = x[gridLoc] = 0.0;

    rhon = alfa = omega = 1.0;

    int inthere=0;
    for(int ip=0; ip < PR_ITE; ip++)
    {
        rho = rhon;
        rhon = multi_matrix(r0,r); //cout<<"\n\n rhon  "<<rhon;

				if(ip==0){
					for(int k=0; k<nsq; k++){
						if( k == gridLoc ) continue;
						p[k] = r[k]; }
				}	else{
					beta=(rhon/rho)*(alfa/omega); //cout<<"\n beta  "<<beta;
					for(int k=0; k<nsq; k++) {
						if( k == gridLoc ) continue;
						p[k] = r[k] + beta*(p[k] - omega*vv[k]);
					}
				}

        multiA(vv,p);
        //cout<<"\n new"<<vv[N*N-10];
        alfa = rhon / multi_matrix(r0,vv); //cout<<"\n alfa  "<<alfa;
        for(int k=0; k<nsq; k++){
        	if( k == gridLoc ) continue;
        	s[k] = r[k] - alfa * vv[k]; 
        }
        multiA(t,s);
        omega = multi_matrix(t,s)/multi_matrix(t,t); //cout<<"\n omega  "<<omega;
        for(int k=0; k<nsq; k++)
        {
        		if( k == gridLoc ) continue;
            x[k] += alfa*p[k] + omega * s[k]; //cout<<" xx "<<x[k];
            r[k] = s[k] - omega * t[k];
        }
        //calc_Ax();
        //for(int k=0; k<nsq; k++)
        //{
        	//r[k] = b(k)-Ax[k];
        //}
        //r[fixx-1+(fixy-1)*N] = 0.0;
        //x[fixx-1+(fixy-1)*N] = fix_value;
        //residual();
        if( (multi_matrix(r,r)/(double)(N*N-1)) <MAX_ERR) break;
        inthere=ip; 
    }
    cout<<"\n   BiCGStab iterations   "<<inthere;
    for(int i=0; i<nsq; i++)
    {
        //i++;
        int yc = int((double)i/(double)N)+1,xc = i-(yc-1)*N+1; //cout<<"\n xc "<<xc<<" yc "<<yc<<"  "<<i;
        Pr[xc][yc] = x[i];
    }
    for(int i=1; i<N+1; i++)
    {
        Pr[0][i] =  Pr[1][i];
        Pr[N+1][i] =  Pr[N][i];
        Pr[i][0] =  Pr[i][1];
        Pr[i][N+1] =  Pr[i][N];
    }
    Pr[0][0]=Pr[0][1];Pr[0][N+1]=Pr[0][N];
    Pr[N+1][0]=Pr[N][0];Pr[N+1][N+1]=Pr[N+1][N];
}
void new_velocity()
{
    for(int i=1; i<N; i++)
    {
        for(int j=1; j<=N; j++)
        {
            u[i][j]=u_new[i][j]; v[j][i]=v_new[j][i];
            u_new[i][j]=u_star[i][j]-DELT*( Pr[i+1][j]-Pr[i][j] )/RHO/h;
            v_new[j][i]=v_star[j][i]-DELT*( Pr[j][i+1]-Pr[j][i] )/RHO/h;
        }
    }
    
    apply_boundary( 2 );

}

void temperature_update()
{
	double uconv,vconv,tu,tv,tdiff,source,tempAvg;
	double avgKUP,avgKDN;
	for(int i=1; i<N+1; i++){
		for(int j=1; j<N+1; j++){
			uconv = ( u_new[i][j] + u_new[i-1][j] ) / 2.;
			// u * dT/dx
			if( uconv > 0. ) tu = ( T[i][j] - T[i-1][j] ) / h;
			else tu = ( T[i+1][j] - T[i][j] ) / h;
			tu *= uconv;
			
			vconv = ( v_new[i][j] + v_new[i][j-1] ) / 2.;
			// v * dT/dy
			if( vconv > 0. ) tv = ( T[i][j] - T[i][j-1] ) / h;
			else tv = ( T[i][j+1] - T[i][j] ) / h;
			tv *= vconv;
			
			// ( 1/Cv/RHO ) * d( k * dT/dx )/dx
			tempAvg = (T[i+1][j]+T[i][j])/2.;
			//avgKUP = (conductivity(i+1,j)+conductivity(i,j))/2.;
			avgKUP = conductivity(tempAvg);
			tdiff = avgKUP * (T[i+1][j]-T[i][j]);
			tempAvg = (T[i-1][j]+T[i][j])/2.;
			//avgKDN = (conductivity(i-1,j)+conductivity(i,j))/2.;
			avgKDN = conductivity(tempAvg);
			tdiff += -1. * avgKDN * (T[i][j]-T[i-1][j]);
			
			
			
			// ( 1/Cv/RHO ) * d( k * dT/dy )/dy
			tempAvg = (T[i][j+1]+T[i][j])/2.;
			//avgKUP = (conductivity(i,j+1)+conductivity(i,j))/2.;
			avgKUP = conductivity(tempAvg);
			tdiff += avgKUP * (T[i][j+1]-T[i][j]);
			tempAvg = (T[i][j-1]+T[i][j])/2.;
			//avgKDN = (conductivity(i,j-1)+conductivity(i,j))/2.;
			avgKDN = conductivity(tempAvg);
			tdiff += -1. * avgKDN * (T[i][j]-T[i][j-1]);
			
			
			// divide by (CV*RHO*h*h) as Cv = 3/2
			tdiff /= (CV * RHO * h * h ) ;
			
			// (2 Mu /Cv/RHO)* [ (du/dx)^2 + (dv/dy)^2 + 0.5*(dv/dx + du/dy)^2 ]
			double vbyx = (v_new[i+1][j]+v_new[i+1][j-1]-v_new[i-1][j]-v_new[i-1][j-1])/4.;
			if( i == 1 ) { vbyx = (v_new[i+1][j]+v_new[i+1][j-1]+v_new[i][j]+v_new[i][j-1])/4.0; }
			else if( i == N ) { vbyx = -(v_new[i-1][j]+v_new[i-1][j-1]+v_new[i][j]+v_new[i][j-1])/4.0; }
			
			double ubyy = (u_new[i][j+1]+u_new[i-1][j+1]-u_new[i][j-1]-u_new[i-1][j-1])/4.;
			//if( j == N ) { ubyy = UTOP_Loc(i) - (u_new[i][j]+u_new[i-1][j]+u_new[i][j-1]+u_new[i-1][j-1])/4.0; }
			source = 0.5*(vbyx+ubyy)*(vbyx+ubyy) ;
			double ubyx = (u_new[i][j] - u_new[i-1][j]);
			double vbyy = (v_new[i][j] - v_new[i][j-1]);
			source += ubyx * ubyx + vbyy * vbyy;
			source *= (2.0*VISCO(T[i][j]) / RHO / CV / h / h ); 
			
			
			
			T[i][j] += DELT * (tdiff + source - tv - tu );
		}
	}
	
	apply_boundary( 3 );
	
}
void ave_values()
{
    for(int i=0; i<=N; i++)
    {
			midU[i] = tav1[N/2][i];
			midV[i] = tav2[i][N/2];
    }
    for(int j=0; j<=N; j++)
    {
        for(int i=0; i<=N; i++)
        {
            tav1[i][j] = (u_new[i][j]+u_new[i][j+1])/2.0;
	    			tav2[i][j] = (v_new[i][j]+v_new[i+1][j])/2.0;
        }
    }
}
void write_results()
{
    ofstream print;
    print.open("Value_UV.dat");
    print<<"TITLE = FLOW \nVARIABLES= X, Y, U, V, num, stress, vorti, pres, temp, Mu, k  \nZONE T= Omega I = "<<N<<", J= "<<N;
    for(int i=1; i<=N; i++)
    {
        for(int j=1; j<=N; j++)
        {
            double xx=(double)i*h-h/2.0, yy=(double)j*h-h/2.0;
            double ux=(u_new[i][j]+u_new[i-1][j])/2.0;
	    			double vy=(v_new[i][j]+v_new[i][j-1])/2.0;
            print<<"\n"<<xx<<"\t"<<yy<<"\t"<<0<<"\t"<<ux<<"\t"<<vy<<"\t"<<sigma[i-1][j-1]<<"\t"<<vorticity[i-1][j-1]<<"\t"<<Pr[i-1][j-1]<<"\t"<<T[i-1][j-1]<<"\t"<<VISCO(T[i-1][j-1])<<"\t"<<conductivity(T[i-1][j-1]);
        }
    }

    /*ofstream write;
    write.open("Value_S_P_V.dat");
    write<<"TITLE = FLOW \nVARIABLE= X, Y, stress, Pr, Vorticity \nZONE T= Omega I = "<<N<<", J= "<<N;
    for(int j=1; j<=N; j++)
    {
        for(int i=1; i<=N; i++)
        {
            double xx=(double)i*h-h/2.0, yy=(double)j*h-h/2.0;
            double sig = (sigma[i-1][j-1]+sigma[i-1][j]+sigma[i][j-1]+sigma[i][j])/4.0;
	    double vor = (vorticity[i-1][j-1]+vorticity[i-1][j]+vorticity[i][j-1]+vorticity[i][j])/4.0;
            write<<"\n"<<xx<<"\t"<<yy<<"\t"<<sig<<"\t"<<Pr[i][j]<<"\t"<<vor;
        }
    }*/

		ofstream writewallleft;
		writewallleft.open("Left_wall_stress.dat");
		writewallleft<<"VARIABLES = rByL, PXYleft, PXYtop, T, Mu, k, ut_top";
		for(int i=1; i<N; i++){
			double xx = ((double)i)*h;
			double uut = top_plate_exp( xx );
			//double strVal = VISCO * (u_new[i][N+1]-u_new[i][N]+u_new[i+1][N+1]-u_new[i+1][N])/h/2.;
			writewallleft<<"\n"<< xx <<"\t"<<strLW[N-i]<<"\t"<<strTW[i]<<"\t"<<T[i][N]<<"\t"<<VISCO(T[i][N])<<"\t"<<conductivity(T[i][N])<<"\t"<<uut;
		}
		
		ofstream writewallright;
		writewallright.open("Right_wall_stress.dat");
		writewallright<<"VARIABLES = rByL, right, top";
		for(int i=1; i<N; i++){
			double xx=((double)i)*h;
			//double strVal = VISCO * (u_new[N-i][N+1]-u_new[N-i][N]+u_new[N-i-1][N+1]-u_new[N-i-1][N])/h/2. ;
			writewallright<<"\n"<< xx <<"\t"<<strRW[N-i]<<"\t"<<strTW[N-i];
		}
		
		ofstream writeUw;
		writeUw.open("Wall_uw.dat");
		writeUw<<"VARIABLES = r, leftWall, leftTop, rightTop, rightWall, Cverti, Chori";
		for(int i=1; i<N+1; i++){
			double rr = ((double)i) * h;
			double grad = (u_new[1][N+1-i] - u_new[0][N+1-i])/h;
			grad += (v_new[1][N+1-i]+v_new[1][N+1-i-1] - v_new[0][N+1-i]-v_new[0][N+1-i-1])/2.0/h;
			writeUw << "\n" << rr << "\t" << grad;
			
			grad = (v_new[i][N] - v_new[i][N-1])/h;
			grad += (u_new[i][N+1]-u_new[i-1][N+1] + u_new[i][N]-u_new[i-1][N])/2.0/h;
			writeUw << "\t" << grad;
			
			grad = (v_new[N+1-i][N] - v_new[N+1-i][N-1])/h;
			grad += (u_new[N+1-i][N+1]-u_new[N+1-i-1][N+1] + u_new[N+1-i][N]-u_new[N+1-i-1][N])/2.0/h;
			writeUw << "\t" << grad;
			
			grad = (u_new[N][N+1-i] - u_new[N-1][N+1-i])/h;
			grad += (v_new[N+1][N+1-i]+v_new[N+1][N+1-i-1] - v_new[N][N+1-i]-v_new[N][N+1-i-1])/2.0/h;
			writeUw << "\t" << grad;
			
			grad = (u_new[int(N/2)][N+1-i] - u_new[int(N/2)-1][N+1-i])/h;
			grad += (v_new[int(N/2)][N+1-i] - v_new[int(N/2)][N-1-i-1])/h;
			writeUw << "\t" << grad;
			
			grad = (v_new[i][int(N/2)] - v_new[i][int(N/2)-1])/h;
			grad += (u_new[i][int(N/2)] - u_new[i-1][int(N/2)])/h;
			writeUw << "\t" << grad;
			
			
			
			
			//writeUw << "\n" << rr << "\t" << (v_star[0][N-i]+v_new[1][N-i])/2.0 << "\t" << (u_star[i][N+1]+u_new[i][N])/2.0;
			//writeUw <<"\t" << (u_star[N-i][N+1]+u_new[N-i][N])/2.0 << "\t" << (v_new[N][N-i]+v_star[N+1][N-i])/2.0 ;
		}
		writeUw.close();
		
		ofstream printgen;
		printgen.open("Gen_info.txt");
		printgen<<"\nSystem length : "<<L<<"\t"<<L;
		printgen<<"\nNumber of bins : "<<N;
		printgen<<"\n dx/L ratio : "<<h/L;
		printgen<<"\nShear force on top plate : "<<shearF;
		printgen<<"\nShear Stress on top plate : "<<shearF/L;
}
void calculate_shear()
{
    for(int i=0; i<N; i++)
    {
	for(int j=0; j<N; j++)
	{
		sigma[j][i] = VISCO(T[i][j])*(tav1[j+1][i+1]+tav1[j][i+1]-tav1[j+1][i]-tav1[j][i])/2.0/h;
		sigma[j][i]+= VISCO(T[i][j])*(tav2[j+1][i+1]-tav2[j][i+1]+tav2[j+1][i]-tav2[j][i])/2.0/h;
	}
    }
}
void calculate_vorti()
{
    for(int i=0; i<N; i++)
    {
	for(int j=0; j<N; j++)
	{
		vorticity[j][i] = (tav2[j+1][i+1]-tav2[j][i+1]+tav2[j+1][i]-tav2[j][i])/2.0/h;
		vorticity[j][i]-=(tav1[j+1][i+1]+tav1[j][i+1]-tav1[j+1][i]-tav1[j][i])/2.0/h;
	}
    }
}
void UV_residual()
{
    UV_norm = 0.0;
    for(int i=1; i<N; i++)
    {
	for(int j=1; j<=N; j++)
	{
		double Udiff=u[i][j]-u_new[i][j];
		double Vdiff=v[j][i]-v_new[j][i];
		if(Udiff<0.0) Udiff *= -1.0;
		if(Vdiff<0.0) Vdiff *= -1.0;
		if(Udiff > UV_norm) UV_norm = Udiff;
		if(Vdiff > UV_norm) UV_norm = Vdiff;
	}
    }
}
void calculate_wallStress()
{
	double avgMu;
	for(int i=1; i<N; i++){
		// Mu * dv/dx
		avgMu = (VISCO(T[1][i])+VISCO(T[0][i])+VISCO(T[1][i+1])+VISCO(T[0][i+1]) )/4.;
		strLW[i] = avgMu * (v_new[1][i]-v_new[0][i])/h; //(tav2[1][i])/h;
		avgMu = (VISCO(T[N][i])+VISCO(T[N+1][i])+VISCO(T[N][i+1])+VISCO(T[N+1][i+1]) )/4.;
		strRW[i] = avgMu * (v_new[N+1][i]-v_new[N][i])/h; //(tav2[N-1][i])/h;
		
		// Mu * du/dy
		avgMu = (VISCO(T[i][N])+VISCO(T[i+1][N])+VISCO(T[i][N+1])+VISCO(T[i+1][N+1]) )/4.;
		strTW[i] = avgMu * (u_new[i][N+1]-u_new[i][N])/h;//(UTOP-tav1[i][N-1])/h;		// shear stress on top plate
		avgMu = (VISCO(T[i][1])+VISCO(T[i+1][1])+VISCO(T[i][0])+VISCO(T[i+1][0]) )/4.;
		strBW[i] = (u_new[i][N+1]-u_new[i][N])/h; //(tav1[i][1])/h; // shear rate on top plate
	}
	
		
	// simpson's 3/8 rule
	shearF = strTW[1]  ;
	shearF += strTW[N-1]  ;
	for(int i=1; i<N-2; i++) {
		if(i%3 == 0) shearF += 2. * strTW[i+1]  ;
		else shearF += 3. * strTW[i+1]  ;
	}
	shearF *= ( h * 3. / 8. ) ;
	shearF += strTW[1] * h  * 0.5;
	shearF += strTW[N-1] * h  * 0.5;
	
}
double vel_convergence()
{
		double tmpU = 0.0, tmpV = 0.0;
    for(int i=1; i<N; i++)
    {
        for(int j=1; j<=N; j++)
        {
            tmpU += (u[i][j] - u_new[i][j])*(u[i][j] - u_new[i][j]); 
            tmpV += (v[j][i] - v_new[j][i])*(v[j][i] - v_new[j][i]);
        }
    }
    cout << " residual \n";
    return sqrt( (tmpU + tmpV)/(2.0 * N * N) );

}

int main()
{
		clock_t start,finish; 
    double total; 
    start = clock();
    
    initialization();
    provisional_velocity();
    pressure_iteration();
    new_velocity();
		ofstream write1D;  write1D.open("W_heat.dat");  write1D.close();
	

    for(int ii=1; ii<ITE+2; ii++)
    {
        provisional_velocity();
        pressure_iteration();
        new_velocity();
        //temperature_update();
        cout<<"\n"<<ii;
				timeNow += DELT;
				if ( ii % 1000  == 0 ) {
					if( vel_convergence() < ERR_UV ) break;
				}
    }
    ave_values();
    calculate_shear();
    calculate_vorti();
    calculate_wallStress();
    write_results();
    
    
    deallocate(u, N+1, N+2);
		deallocate(u_new, N+1, N+2);
		deallocate(u_star, N+1, N+2);
		deallocate(v, N+2, N+1);
		deallocate(v_new, N+2, N+1);
		deallocate(v_star, N+2, N+1);
		deallocate(Pr, N+2, N+2);
		deallocate(sigma, N, N);
		deallocate(vorticity, N, N);
		deallocate(tav1, N+1, N+1);
		deallocate(tav2, N+1, N+1);
		deallocate(vv, N*N); deallocate(x, N*N); deallocate(r, N*N); deallocate(r0, N*N); deallocate(p, N*N);
		deallocate(t, N*N); deallocate(s, N*N); deallocate(Ax, N*N);
		deallocate(midU, N+1); deallocate(midV, N+1);
		deallocate(strLW, N); deallocate(strTW, N); deallocate(strRW, N); deallocate(strBW, N);
		deallocate(T, N+2, N+2);
		deallocate(LeftBou,N,3); deallocate(RightBou,N,3); deallocate(TopBou,N,3); deallocate(BottomBou,N,3); 
		
		finish = clock(); total = (double) (finish - start) / (double) CLOCKS_PER_SEC;
    printf("Total = %f\n",total);
}
