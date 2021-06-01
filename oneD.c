																																												
// gcc EarthAndSunTest.c -o temp -lglut -lm -lGLU -lGL
//To stop hit "control c" in the window you launched it from.
#include <GL/glut.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define PI 3.141592654

#define XWindowSize 2000
#define YWindowSize 2000

#define STOP_TIME 30.0
//#define DT  0.000000001
#define   DT  0.00000001
#define N 3

#define GRAVITY 1.0 

#define MASS_OF_SUN 332946.0
#define MASS_OF_EARTH 1.0 
#define EARTH_DIAMETER 0.0000851738
#define EARTH_INT_VELOCITY 576.912 
#define EARTH_TO_SUN_DIS 1.0
#define ANGLE1 0.0
#define ANGLE2 (PI/2.0)	

#define DRAW 10000
#define PRINTS 100000.0
#define WRITE_DATA (STOP_TIME/PRINTS)

// Globals
double px[N], py[N], pz[N], vx[N], vy[N], vz[N], fx[N], fy[N], fz[N], mass[N], G, initail_angular_speed;
double pxold[N], pyold[N], pzold[N]; 
FILE *data_file, *data_file1, *data_file2;
int i;

void set_initail_conditions()
{
    double temp, tempi, tempj, tempk;

	G = GRAVITY;
	mass[0] = MASS_OF_SUN;
	mass[1] = 1.0*MASS_OF_EARTH;
	mass[2] = 1.0*MASS_OF_EARTH;
	
	px[0] = 0.0;
	py[0] = 0.0;
	pz[0] = 0.0;
	
	px[1] = EARTH_TO_SUN_DIS*cos(ANGLE1);
	py[1] = EARTH_TO_SUN_DIS*sin(ANGLE1);
	pz[1] = 0.0;
	
	px[2] = 1.0*EARTH_TO_SUN_DIS*cos(ANGLE2);
	py[2] = 1.0*EARTH_TO_SUN_DIS*sin(ANGLE2);
	pz[2] = 0.0;
	
	vx[0] = 0.0;
	vy[0] = 0.0;
	vz[0] = 0.0;
	
	vx[1] = 1.0*EARTH_INT_VELOCITY*cos(ANGLE1 + PI/2.0);
	vy[1] = 1.0*EARTH_INT_VELOCITY*sin(ANGLE1 + PI/2.0);
	vz[1] = 0.0;
	
	vx[2] = EARTH_INT_VELOCITY*cos(ANGLE2 + PI/2.0);
	vy[2] = EARTH_INT_VELOCITY*sin(ANGLE2 + PI/2.0);
	vz[2] = 0.0;
	
//figure 8	
//initial conditions: (x1,y1) = (-0.97000436, 0.24308753), (x2,y2) = (-x1, -y1), (x3,y3) = (0,0) 
//initial velocites (vx1,vy1) = (vx2, vy2) = -(vx3, vy3)/2; where (vx3,vy3) = (0.93240737, 0.86473146) 
/*
	mass[0] = 1.0;
	mass[1] = 1.0;
	mass[2] = 1.0;
	
	px[0] = -0.97000436;
	py[0] = 0.24308753;
	pz[0] = 0.0;
	
	px[1] = 0.97000436;
	py[1] = -0.24308753;
	pz[1] = 0.0;
	
	px[2] = 0.0;
	py[2] = 0.0;
	pz[2] = 0.0;
	
	vx[0] = -0.93240737/2.0;
	vy[0] = -0.86473146/2.0;
	vz[0] = 0.0;
	
	vx[1] = -0.93240737/2.0;
	vy[1] = -0.86473146/2.0;
	vz[1] = 0.0;
	
	vx[2] = 0.93240737;
	vy[2] = 0.86473146;
	vz[2] = 0.0;
*/
	//for the Varlet
	for(i=0; i<3; i++)
	{
		pxold[i] = px[i];
		pyold[i] = py[i];
		pzold[i] = pz[i];
	}
		 	

    //Calculating the angular velosities of the planets
	tempi =   (py[1]*vz[1] - pz[1]*vy[1]);
	tempj = -((px[1]*vz[1] - pz[1]*vx[1]));
	tempk =   (px[1]*vy[1] - py[1]*vx[1]);

	temp = tempi*tempi + tempj*tempj + tempk*tempk;
	initail_angular_speed = sqrt(temp);
	temp = px[1]*px[1] + py[1]*py[1] + pz[1]*pz[1];
	initail_angular_speed = initail_angular_speed/temp;
}

void print_out_values(double time)
{
	double temp, tempi, tempj, tempk;
	double delta_v1, delta_v2;
	double px1, py1, pz1;
	double px2, py2, pz2;
	double vx1, vy1, vz1;
	double vx2, vy2, vz2;
	double seperation_distance;
	double closing_velocity;
	double angle_between_planets;
	double delta_orbital_distance1, delta_orbital_distance2;
	double angular_speed1, angular_speed2;
	
	//Raw planet velocities calculations
	tempi = vx[1]-vx[2];
	tempj = vy[1]-vy[2];
	tempk = vz[1]-vz[2];

	temp = tempi*tempi + tempj*tempj + tempk*tempk;
	closing_velocity = sqrt(temp);

	delta_v1 = sqrt(vx[1]*vx[1] + vy[1]*vy[1] + vz[1]*vz[1]) - EARTH_INT_VELOCITY; //Change in initial velocity of planet 1

	delta_v2 = sqrt(vx[2]*vx[2] + vy[2]*vy[2] + vz[2]*vz[2]) - EARTH_INT_VELOCITY; //Change in initial velocity of planet 2
	
	//Setting everything relative to the Sun's position and velocity
	px1 = px[1]-px[0];
	py1 = py[1]-py[0];
	pz1 = pz[1]-pz[0];

	px2 = px[2]-px[0];
	py2 = py[2]-py[0];
	pz2 = pz[2]-pz[0];

	vx1 = vx[1]-vx[0];
	vy1 = vy[1]-vy[0];
	vz1 = vz[1]-vz[0];

	vx2 = vx[2]-vx[0];
	vy2 = vy[2]-vy[0];
	vz2 = vz[2]-vz[0];

	//Calculating seperation of planets
 	tempi = px1-px2;
	tempj = py1-py2;
	tempk = pz1-pz2;
	temp = tempi*tempi + tempj*tempj + tempk*tempk;
	seperation_distance = sqrt(temp);

	//Calculating planet's closing velosity
	tempi = vx1-vx2;
	tempj = vy1-vy2;
	tempk = vz1-vz2;
	temp = tempi*tempi + tempj*tempj + tempk*tempk;
	closing_velocity = sqrt(temp);

	//Calculation the change in planets velocities
	delta_v1 = sqrt(vx1*vx1 + vy1*vy1 + vz1*vz1) - EARTH_INT_VELOCITY;
	delta_v2 = sqrt(vx2*vx2 + vy2*vy2 + vz2*vz2) - EARTH_INT_VELOCITY;
	
	//Calculating the angle between the planets
	angle_between_planets = acos(px1*px2 + py1*py2 + pz1*pz2)/(sqrt(px1*px1 + py1*py1 + pz1*pz1)*sqrt(px2*px2 + py2*py2 + pz2*pz2));

	//Calculating the change in the orbital distances
	temp = px1*px1 + py1*py1 + pz1*pz1;
	delta_orbital_distance1 = sqrt(temp) - EARTH_TO_SUN_DIS;

	temp = px2*px2 + py2*py2 + pz2*pz2;
	delta_orbital_distance2 = sqrt(temp) - EARTH_TO_SUN_DIS;

	//Calculating the angular velosities of the planets
	tempi =   (py1*vz1 - pz1*vy1);
	tempj = -((px1*vz1 - pz1*vx1));
	tempk =   (px1*vy1 - py1*vx1);

	temp = tempi*tempi + tempj*tempj + tempk*tempk;
	angular_speed1 = sqrt(temp);
	temp = px1*px1 + py1*py1 + pz1*pz1;
	angular_speed1 = angular_speed1/temp - initail_angular_speed;

	tempi =   (py2*vz2 - pz2*vy2);
	tempj = -((px2*vz2 - pz2*vx2));
	tempk =   (px2*vy2 - py2*vx2);

	temp = tempi*tempi + tempj*tempj + tempk*tempk;
	angular_speed2 = sqrt(temp);
	temp = px2*px2 + py2*py2 + pz2*pz2;
	angular_speed2 = angular_speed2/temp - initail_angular_speed;


	printf("\n years = %lf", time*91.84280);
	//printf("\n closing velocity = %lf km/hr seperation = %lf km" closing_velocity*185.817, seperation_distance*149600000.0);
	//printf("\n orbit of planet1 = %lf orbit of planet2 = %lf", delta_orbital_distance1*149600000.0, delta_orbital_distance2*149600000.0);
	printf("\n angle between planets = %lf", angle_between_planets*180.0/PI);
	//printf("\n angle velocity planets1 = %lf  angle velocity planets2 = %lf", angular_speed1, angular_speed2);


	//printf("\n Time = %lf", time);
	//printf("\n Angle between planets in degrees = %lf", angle_between_planets*180.0/PI);
	//printf("\n Seperation between planets= %lf", seperation_distance);
	//printf("\n closing velocity = %lf\n", closing_velocity);
	//printf("\n Delta orbit of planet1 = %lf Delta orbit of planet2 = %lf", delta_orbital_distance1, delta_orbital_distance2);
	//printf("\n Delta velocity planets1 = %lf  Delta velocity planets2 = %lf", delta_v1, delta_v2);
	//printf("\n angle velocity planets1 = %lf  angle velocity planets2 = %lf", angular_speed1, angular_speed2);

	//printf("\n         Orbital  difference between plant1 to planet2 = %lf", delta_orbital_distance1 - delta_orbital_distance2);
	//printf("\n         Velosity difference between plant1 to planet2 = %lf", delta_v1 - delta_v2);
	//printf("\n Angular Velosity difference between plant1 to planet2 = %lf", angular_speed1 - angular_speed2);

    fprintf(data_file1, "%lf %lf %lf %lf %lf\n", time, angle_between_planets*180.0/PI, delta_orbital_distance1, delta_v1, angular_speed1);
    fprintf(data_file2, "%lf %lf %lf %lf %lf\n", time, angle_between_planets*180.0/PI, delta_orbital_distance2, delta_v2, angular_speed2);
    fprintf(data_file,  "%lf %lf %lf %lf %lf\n", time, angle_between_planets*180.0/PI, delta_orbital_distance1 - delta_orbital_distance2, delta_v1 - delta_v2, angular_speed1 - angular_speed2);
//printf("%lf %lf %lf %lf %lf\n", time, angle_between_planets*180.0/PI, delta_orbital_distance1 - delta_orbital_distance2, delta_v1 - delta_v2, angular_speed1 - angular_speed2 );
}

void draw_picture()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glClear(GL_DEPTH_BUFFER_BIT);
	
	//gluLookAt(0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
	
	glColor3d(1.0,1.0,0.5);
	glPushMatrix();
	glTranslatef(px[0], py[0], pz[0]);
	glutSolidSphere(0.1,20,20);
	glPopMatrix();
	
	glColor3d(1.0,0.0,0.0);
	glPushMatrix();
	glTranslatef(px[1], py[1], pz[1]);
	glutSolidSphere(0.03,20,20);
	glPopMatrix();
	
	glColor3d(0.0,0.0,1.0);
	glPushMatrix();
	glTranslatef(px[2], py[2], pz[2]);
	glutSolidSphere(0.03,20,20);
	glPopMatrix();
	
	glutSwapBuffers();
}

euler(double dt)
{
	double fx[3], fy[3], fz[3], f; 
	double dx,dy,dz,d;
	double dvx,dvy,dvz;
	int i,j;
	
	for(i=0; i<3; i++)
	{
		fx[i] = 0.0;
		fy[i] = 0.0;
		fz[i] = 0.0;
	}

	for(i=0; i<3; i++)
	{
		for(j=i+1; j<3; j++)
		{
			dx = px[j]-px[i];
			dy = py[j]-py[i];
			dz = pz[j]-pz[i];
			d  = sqrt(dx*dx + dy*dy + dz*dz);

			if(d < EARTH_DIAMETER)
			{
				return(1);
			}
			else
			{
				f  = (GRAVITY*mass[i]*mass[j])/(d*d);
			}

			fx[i] += f*dx/d;
			fx[j] -= f*dx/d;
			fy[i] += f*dy/d;
			fy[j] -= f*dy/d;
			fz[i] += f*dz/d;
			fz[j] -= f*dz/d;
		}
	}
}
	
improvedEuler(double dt)
{
	double fx[3], fy[3], fz[3],fxEuler[3], fyEuler[3], fzEuler[3],f; 
	double pxEuler[3], pyEuler[3], pzEuler[3], vxEuler[3], vyEuler[3], vzEuler[3];
	
	double dx,dy,dz,d;
	double dvx,dvy,dvz;
	int i,j;
	
	for(i=0; i<3; i++)
	{
		fx[i] = 0.0;
		fy[i] = 0.0;
		fz[i] = 0.0;
	}

	for(i=0; i<3; i++)
	{
		for(j=i+1; j<3; j++)
		{
			dx = px[j]-px[i];
			dy = py[j]-py[i];
			dz = pz[j]-pz[i];
			d  = sqrt(dx*dx + dy*dy + dz*dz);

			if(d < EARTH_DIAMETER)
			{
				return(1);
			}
			else
			{
				f  = (GRAVITY*mass[i]*mass[j])/(d*d);
			}

			fx[i] += f*dx/d;
			fx[j] -= f*dx/d;
			fy[i] += f*dy/d;
			fy[j] -= f*dy/d;
			fz[i] += f*dz/d;
			fz[j] -= f*dz/d;
		}
	}

	for(i=0; i<3; i++)
	{
		pxEuler[i] = px[i] + vx[i]*dt;
		pyEuler[i] = py[i] + vy[i]*dt;
		pzEuler[i] = pz[i] + vz[i]*dt;
		
		vxEuler[i] = vx[i] + (fx[i]/mass[i])*dt;
		vyEuler[i] = vy[i] + (fy[i]/mass[i])*dt;
		vzEuler[i] = vz[i] + (fz[i]/mass[i])*dt;
	}
	
	for(i=0; i<3; i++)
	{
		fxEuler[i] = 0.0;
		fyEuler[i] = 0.0;
		fzEuler[i] = 0.0;
	}

	for(i=0; i<3; i++)
	{
		for(j=i+1; j<3; j++)
		{
			dx = pxEuler[j]-pxEuler[i];
			dy = pyEuler[j]-pyEuler[i];
			dz = pzEuler[j]-pzEuler[i];
			d  = sqrt(dx*dx + dy*dy + dz*dz);

			if(d < EARTH_DIAMETER)
			{
				return(1);
			}
			else
			{
				f  = (GRAVITY*mass[i]*mass[j])/(d*d);
			}

			fxEuler[i] += f*dx/d;
			fxEuler[j] -= f*dx/d;
			fyEuler[i] += f*dy/d;
			fyEuler[j] -= f*dy/d;
			fzEuler[i] += f*dz/d;
			fzEuler[j] -= f*dz/d;
		}
	}
	
	for(i=0; i<3; i++)
	{
		px[i] += (vx[i] + vxEuler[i])*(dt/2.0);
		py[i] += (vy[i] + vyEuler[i])*(dt/2.0);
		pz[i] += (vz[i] + vzEuler[i])*(dt/2.0);
		
		vx[i] += (fx[i]/mass[i] + fxEuler[i]/mass[i])*(dt/2.0);
		vy[i] += (fy[i]/mass[i] + fyEuler[i]/mass[i])*(dt/2.0);
		vz[i] += (fz[i]/mass[i] + fzEuler[i]/mass[i])*(dt/2.0);
	}
}

RK4(double dt)
{
	double fx[3], fy[3], fz[3],f;
	double m0x[3],m0y[3],m0z[3]; 
	double m1x[3],m1y[3],m1z[3];
	double m2x[3],m2y[3],m2z[3]; 
	double m3x[3],m3y[3],m3z[3]; 
	double tempx[3], tempy[3], tempz[3];
	
	double dx,dy,dz,d;
	double dvx,dvy,dvz;
	int i,j;
	
	for(i=0; i<3; i++)
	{
		fx[i] = 0.0;
		fy[i] = 0.0;
		fz[i] = 0.0;
		tempx[i] = px[i];
		tempy[i] = py[i];
		tempz[i] = pz[i];
	}

	for(i=0; i<3; i++)
	{
		for(j=i+1; j<3; j++)
		{
			dx = tempx[j]-tempx[i];
			dy = tempy[j]-tempy[i];
			dz = tempz[j]-tempz[i];
			d  = sqrt(dx*dx + dy*dy + dz*dz);

			if(d < EARTH_DIAMETER)
			{
				return(1);
			}
			else
			{
				f  = (GRAVITY*mass[i]*mass[j])/(d*d);
			}

			fx[i] += f*dx/d;
			fx[j] -= f*dx/d;
			fy[i] += f*dy/d;
			fy[j] -= f*dy/d;
			fz[i] += f*dz/d;
			fz[j] -= f*dz/d;
		}
	}

	for(i=0; i<3; i++)
	{
		m0x[i] = (fx[i]/mass[i])*dt;
		m0y[i] = (fy[i]/mass[i])*dt;
		m0z[i] = (fz[i]/mass[i])*dt;
	}
	
	for(i=0; i<3; i++)
	{
		fx[i] = 0.0;
		fy[i] = 0.0;
		fz[i] = 0.0;
		tempx[i] = px[i] + (dt/2.0)*vx[i];
		tempy[i] = py[i] + (dt/2.0)*vy[i];
		tempz[i] = pz[i] + (dt/2.0)*vz[i];
	}

	for(i=0; i<3; i++)
	{
		for(j=i+1; j<3; j++)
		{
			dx = tempx[j]-tempx[i];
			dy = tempy[j]-tempy[i];
			dz = tempz[j]-tempz[i];
			d  = sqrt(dx*dx + dy*dy + dz*dz);

			if(d < EARTH_DIAMETER)
			{
				return(1);
			}
			else
			{
				f  = (GRAVITY*mass[i]*mass[j])/(d*d);
			}

			fx[i] += f*dx/d;
			fx[j] -= f*dx/d;
			fy[i] += f*dy/d;
			fy[j] -= f*dy/d;
			fz[i] += f*dz/d;
			fz[j] -= f*dz/d;
		}
	}
	
	for(i=0; i<3; i++)
	{	
		m1x[i] = (fx[i]/mass[i])*dt;
		m1y[i] = (fy[i]/mass[i])*dt;
		m1z[i] = (fz[i]/mass[i])*dt;
	}
	
	for(i=0; i<3; i++)
	{
		fx[i] = 0.0;
		fy[i] = 0.0;
		fz[i] = 0.0;
		tempx[i] = px[i] + (dt/2.0)*vx[i] + (dt/4.0)*m0x[i];
		tempy[i] = py[i] + (dt/2.0)*vy[i] + (dt/4.0)*m0y[i];
		tempz[i] = pz[i] + (dt/2.0)*vz[i] + (dt/4.0)*m0z[i];
	}

	for(i=0; i<3; i++)
	{
		for(j=i+1; j<3; j++)
		{
			dx = tempx[j]-tempx[i];
			dy = tempy[j]-tempy[i];
			dz = tempz[j]-tempz[i];
			d  = sqrt(dx*dx + dy*dy + dz*dz);

			if(d < EARTH_DIAMETER)
			{
				return(1);
			}
			else
			{
				f  = (GRAVITY*mass[i]*mass[j])/(d*d);
			}

			fx[i] += f*dx/d;
			fx[j] -= f*dx/d;
			fy[i] += f*dy/d;
			fy[j] -= f*dy/d;
			fz[i] += f*dz/d;
			fz[j] -= f*dz/d;
		}
	}
	
	for(i=0; i<3; i++)
	{	
		m2x[i] = (fx[i]/mass[i])*dt;
		m2y[i] = (fy[i]/mass[i])*dt;
		m2z[i] = (fz[i]/mass[i])*dt;
	}
	
	for(i=0; i<3; i++)
	{
		fx[i] = 0.0;
		fy[i] = 0.0;
		fz[i] = 0.0;
		tempx[i] = px[i] + dt*vx[i] + (dt/2.0)*m1x[i];
		tempy[i] = py[i] + dt*vy[i] + (dt/2.0)*m1y[i];
		tempz[i] = pz[i] + dt*vz[i] + (dt/2.0)*m1z[i];
	}

	for(i=0; i<3; i++)
	{
		for(j=i+1; j<3; j++)
		{
			dx = tempx[j]-tempx[i];
			dy = tempy[j]-tempy[i];
			dz = tempz[j]-tempz[i];
			d  = sqrt(dx*dx + dy*dy + dz*dz);

			if(d < EARTH_DIAMETER)
			{
				return(1);
			}
			else
			{
				f  = (GRAVITY*mass[i]*mass[j])/(d*d);
			}

			fx[i] += f*dx/d;
			fx[j] -= f*dx/d;
			fy[i] += f*dy/d;
			fy[j] -= f*dy/d;
			fz[i] += f*dz/d;
			fz[j] -= f*dz/d;
		}
	}
	
	for(i=0; i<3; i++)
	{	
		m3x[i] = (fx[i]/mass[i])*dt;
		m3y[i] = (fy[i]/mass[i])*dt;
		m3z[i] = (fz[i]/mass[i])*dt;
	}
	

	for(i=0; i<3; i++)
	{
		px[i] = px[i] + dt*vx[i] + (dt/6.0)*(m0x[i] + m1x[i] + m2x[i]);
		py[i] = py[i] + dt*vy[i] + (dt/6.0)*(m0y[i] + m1y[i] + m2y[i]);
		pz[i] = pz[i] + dt*vz[i] + (dt/6.0)*(m0z[i] + m1z[i] + m2z[i]);
		
		vx[i] = vx[i] + (1.0/6.0)*(m0x[i] + 2.0*m1x[i] + 2.0*m2x[i] + m3x[i]);
		vy[i] = vy[i] + (1.0/6.0)*(m0y[i] + 2.0*m1y[i] + 2.0*m2y[i] + m3y[i]);
		vz[i] = vz[i] + (1.0/6.0)*(m0z[i] + 2.0*m1z[i] + 2.0*m2z[i] + m3z[i]);
	}
}

Verlet(double dt, double time)
{
	double fx[3], fy[3], fz[3], f;
	double tempx[3], tempy[3], tempz[3]; 
	double dx,dy,dz,d;
	double dvx,dvy,dvz;
	int i,j;
	
	for(i=0; i<3; i++)
	{
		fx[i] = 0.0;
		fy[i] = 0.0;
		fz[i] = 0.0;
	}

	for(i=0; i<3; i++)
	{
		for(j=i+1; j<3; j++)
		{
			dx = px[j]-px[i];
			dy = py[j]-py[i];
			dz = pz[j]-pz[i];
			d  = sqrt(dx*dx + dy*dy + dz*dz);

			if(d < EARTH_DIAMETER)
			{
				return(1);
			}
			else
			{
				f  = (GRAVITY*mass[i]*mass[j])/(d*d);
			}

			fx[i] += f*dx/d;
			fx[j] -= f*dx/d;
			fy[i] += f*dy/d;
			fy[j] -= f*dy/d;
			fz[i] += f*dz/d;
			fz[j] -= f*dz/d;
		}
	}

	for(i=0; i<3; i++)
	{	 
		if(time == 0.0)
		{
			px[i] = px[i] + vx[i]*dt + (fx[i]/mass[i])*0.5*dt*dt;
			py[i] = py[i] + vy[i]*dt + (fy[i]/mass[i])*0.5*dt*dt;
			pz[i] = pz[i] + vz[i]*dt + (fz[i]/mass[i])*0.5*dt*dt;
		}
		else
		{
			tempx[i] = px[i];
			tempy[i] = py[i];
			tempz[i] = pz[i];
			px[i] = 2.0*px[i] - pxold[i] + (fx[i]/mass[i])*dt*dt;
			py[i] = 2.0*py[i] - pyold[i] + (fy[i]/mass[i])*dt*dt;
			pz[i] = 2.0*pz[i] - pzold[i] + (fz[i]/mass[i])*dt*dt;
			pxold[i] = tempx[i];
			pyold[i] = tempy[i];
			pzold[i] = tempz[i];
		}
	}
}

leapFrog(double dt, double time)
{
	double fx[3], fy[3], fz[3], f; 
	double dx,dy,dz,d;
	double dvx,dvy,dvz;
	int i,j;
	
	for(i=0; i<3; i++)
	{
		fx[i] = 0.0;
		fy[i] = 0.0;
		fz[i] = 0.0;
	}

	for(i=0; i<3; i++)
	{
		for(j=i+1; j<3; j++)
		{
			dx = px[j]-px[i];
			dy = py[j]-py[i];
			dz = pz[j]-pz[i];
			d  = sqrt(dx*dx + dy*dy + dz*dz);

			if(d < EARTH_DIAMETER)
			{
				return(1);
			}
			else
			{
				f  = (GRAVITY*mass[i]*mass[j])/(d*d);
			}

			fx[i] += f*dx/d;
			fx[j] -= f*dx/d;
			fy[i] += f*dy/d;
			fy[j] -= f*dy/d;
			fz[i] += f*dz/d;
			fz[j] -= f*dz/d;
		}
	}

	for(i=0; i<3; i++)
	{
		if(time == 0.0)
		{
			vx[i] += (fx[i]/mass[i])*0.5*dt;
			vy[i] += (fy[i]/mass[i])*0.5*dt;
			vz[i] += (fz[i]/mass[i])*0.5*dt;
		}
		else
		{
			vx[i] += (fx[i]/mass[i])*dt;
			vy[i] += (fy[i]/mass[i])*dt;
			vz[i] += (fz[i]/mass[i])*dt;
		}

		px[i] += vx[i]*dt;
		py[i] += vy[i]*dt;
		pz[i] += vz[i]*dt;
	}
}

conservativeModel(double dt, double time)
{
	double fx[3], fy[3], fz[3], f; 
	double dx,dy,dz,d;
	double dvx,dvy,dvz;
	double pxNew[3], pyNew[3], pzNew[3];
	double vxNew[3], vyNew[3], vzNew[3];
	double pxOld[3], pyOld[3], pzOld[3];
	double vxOld[3], vyOld[3], vzOld[3];
	double w =1.0;
	int i,j;
	
	
	for(i=0; i<3; i++)
	{
		pxOld[i] = px[i];
		pyOld[i] = py[i];
		pzOld[i] = pz[i];
		
		vxOld[i] = vx[i];
		vyOld[i] = vy[i];
		vzOld[i] = vz[i];
	}
	
	while(tolerance < test)
	{
		for(i=0; i<3; i++)
		{
			pxNew[i] = pxOld[i];
			pyNew[i] = pyOld[i];
			pzNew[i] = pzOld[i];
		
			vxNew[i] = vxOld[i];
			vyNew[i] = vyOld[i];
			vzNew[i] = vzOld[i];
		}
		
		for(i=0; i<3; i++)
		{
			pxNew[i] = pxNew[i] - w*(pxNew[i] - px[i] -dt*(vxNew[i] + vx[i])/2.0);
			pyNew[i] = pyNew[i] - w*(pyNew[i] - py[i] -dt*(vyNew[i] + vy[i])/2.0);
			pzNew[i] = pzNew[i] - w*(pzNew[i] - pz[i] -dt*(vzNew[i] + vz[i])/2.0);
		
			vxNew[i] = vxNew[i] - w*();
			vyNew[i] = vyNew[i] - w*();
			vzNew[i] = vzNew[i] - w*();
		}
	}
	
	
	
	
	for(i=0; i<3; i++)
	{
		fx[i] = 0.0;
		fy[i] = 0.0;
		fz[i] = 0.0;
	}

	for(i=0; i<3; i++)
	{
		for(j=i+1; j<3; j++)
		{
			dx = px[j]-px[i];
			dy = py[j]-py[i];
			dz = pz[j]-pz[i];
			d  = sqrt(dx*dx + dy*dy + dz*dz);

			if(d < EARTH_DIAMETER)
			{
				return(1);
			}
			else
			{
				f  = (GRAVITY*mass[i]*mass[j])/(d*d);
			}

			fx[i] += f*dx/d;
			fx[j] -= f*dx/d;
			fy[i] += f*dy/d;
			fy[j] -= f*dy/d;
			fz[i] += f*dz/d;
			fz[j] -= f*dz/d;
		}
	}

	for(i=0; i<3; i++)
	{
		if(time == 0.0)
		{
			vx[i] += (fx[i]/mass[i])*0.5*dt;
			vy[i] += (fy[i]/mass[i])*0.5*dt;
			vz[i] += (fz[i]/mass[i])*0.5*dt;
		}
		else
		{
			vx[i] += (fx[i]/mass[i])*dt;
			vy[i] += (fy[i]/mass[i])*dt;
			vz[i] += (fz[i]/mass[i])*dt;
		}

		px[i] += vx[i]*dt;
		py[i] += vy[i]*dt;
		pz[i] += vz[i]*dt;
	}
}


int n_body()
{
	//double fx[3], fy[3], fz[3], f; 
	//double dx,dy,dz,d, dt;
	//double dvx,dvy,dvz,close_seperate;
	int    tdraw = 0; 
    double tprint = 0.0;
	double  time = 0.0;
	//int i,j;

    data_file = fopen("data", "w");
    data_file1 = fopen("data1", "w");
    data_file2 = fopen("data2", "w");
    if (data_file == NULL) 
    {
         printf("I couldn't open data_file.\n");
         exit(0);
    }
    if (data_file1 == NULL) 
    {
         printf("I couldn't open data_file1.\n");
         exit(0);
    }
    if (data_file2 == NULL) 
    {
         printf("I couldn't open data_file2.\n");
         exit(0);
    }
    fprintf(data_file1, "%lf\n", STOP_TIME);

	while(time < STOP_TIME)
	{
		//euler(DT);
		//improvedEuler(DT);
		//leapFrog(DT, time);
		//RK4(DT);
		Verlet(DT, time);

		if(tdraw == DRAW) 
		{
			draw_picture();
			tdraw = 0;
		}

		if(tprint >= WRITE_DATA) 
		{
			print_out_values(time);
			tprint = 0.0;
		}

		time += DT;
		tdraw++;
		tprint += DT;
	}
    fclose(data_file);
    fclose(data_file1);
    fclose(data_file2);
	return(2);
}

void control()
{	
	int    tdraw = 0;
	double  time = 0.0;

	glClear(GL_COLOR_BUFFER_BIT);
	glClear(GL_DEPTH_BUFFER_BIT);

	set_initail_conditions();
	
	draw_picture();
	
	if(n_body() == 1) printf("\n Done Planets hit \n");
	else printf("\n Done timed out \n");
	
	printf("\n DONE \n");
	while(1);
}

void Display(void)
{
	gluLookAt(0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT);
	glClear(GL_DEPTH_BUFFER_BIT);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glutSwapBuffers();
	glFlush();
	control();
}

void reshape(int w, int h)
{
	glViewport(0, 0, (GLsizei) w, (GLsizei) h);

	glMatrixMode(GL_PROJECTION);

	glLoadIdentity();

	glFrustum(-0.2, 0.2, -0.2, 0.2, 0.2, 50.0);

	glMatrixMode(GL_MODELVIEW);
}

int main(int argc, char** argv)
{
	glutInit(&argc,argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGB);
	glutInitWindowSize(XWindowSize,YWindowSize);
	glutInitWindowPosition(0,0);
	glutCreateWindow("2 Body 3D");
	GLfloat light_position[] = {1.0, 1.0, 1.0, 0.0};
	GLfloat light_ambient[]  = {0.0, 0.0, 0.0, 1.0};
	GLfloat light_diffuse[]  = {1.0, 1.0, 1.0, 1.0};
	GLfloat light_specular[] = {1.0, 1.0, 1.0, 1.0};
	GLfloat lmodel_ambient[] = {0.2, 0.2, 0.2, 1.0};
	GLfloat mat_specular[]   = {1.0, 1.0, 1.0, 1.0};
	GLfloat mat_shininess[]  = {10.0};
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glShadeModel(GL_SMOOTH);
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_DEPTH_TEST);
	glutDisplayFunc(Display);
	glutReshapeFunc(reshape);
	glutMainLoop();
	return 0;
}






