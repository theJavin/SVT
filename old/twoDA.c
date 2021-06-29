// gcc oneD.c -o temp -lglut -lm -lGLU -lGL
//To stop hit "control c" in the window you launched it from. stuff

// Length will be in millimeters
// Time will be in milliseconds
// Mass will be in ???

// Fiber length 100 micrometers or 0.1 millimeters
// Sodium wave speed .5 meters/sec or 0.5 millimeters/millisec
#include <GL/glut.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define PI 3.141592654

#define XWindowSize 2000
#define YWindowSize 2000 

#define STOP_TIME 60000.0
#define DT  0.004
#define N 20
	
#define DRAW 1000

// Globals
float Px[N], Vx[N], Fx[N], Mass[N];
float Py[N], Vy[N], Fy[N];

float FiberLength = 0.1;
float FiberStrength = 10.0;
float FiberCompresionMultiplier = 10.0;
float FiberTentionMultiplier = 10.0;
float FiberCompresionStopFraction = 0.6;

float TendonLength = 0.1;
float TendonStrength = 1.0;
float TendonCompresionMultiplier = 10.0;
float TendonTentionMultiplier = 10.0;
float TendonCompresionStopFraction = 0.6;

float APWaveSpeed[N]; //0.5 is a good value.
float APWaveFrontLx, APWaveFrontLy, APWaveFrontRx, APWaveFrontRy;
int APWaveAmunity;

float ContractionStrength[N-1]; // 5.0 is a good value
int ContractionOn[N-1];
float ContractionTime[N-1];
float ContractionDuration[N-1]; // 100.0 is a good value
float RelaxationDuration[N-1]; // 200.0 is a good value

float BeatPeriod = 400.0;

float Viscosity = 10.0;

float Radius = 2.0;
float CentralPush = 5.0;

void set_initial_conditions()
{
	// Node Values
	for(int i = 0; i < N; i++)
	{
		Mass[i] = 1.0;
		//Px[i] = (float)(i+1)*FiberLength - centerX;
		//Px[i] = (float)(i+1)*FiberLength - centerX;
		Px[i] = Radius*cos((float)i*2.0*PI/(float)N+PI/2.0);
		Py[i] = Radius*sin((float)i*2.0*PI/(float)N+PI/2.0);
		Vx[i] = 0.0;
		Vy[i] = 0.0;	
	}
	
	// Muscle (connector) Values
	for(int i = 0; i < N-1; i++)
	{	
		APWaveSpeed[i] = 0.02;
		RelaxationDuration[i] = 200.0;
		ContractionStrength[i] = 5.0;
		ContractionDuration[i] = 100.0;
	}
	
	//RelaxationDuration[3] = 400.0;
	//RelaxationDuration[4] = 400.0;
	//RelaxationDuration[5] = 400.0;
	//RelaxationDuration[6] = 400.0;
}

void draw_picture()
{
	float d, dx, dy;
	glClear(GL_COLOR_BUFFER_BIT);
	glClear(GL_DEPTH_BUFFER_BIT);
	
	// Drawing nodes
	for(int i = 0; i < N; i++)
	{
		glColor3d(1.0,1.0,1.0);
		glPushMatrix();
		glTranslatef(Px[i], Py[i], 0.0);
		glutSolidSphere(0.005,20,20);
		glPopMatrix();	
	}
	
	// Drawing muscles
	glColor3d(1.0,0.0,0.0);
	for(int i = 0; i < N-1; i++)
	{
		dx = Px[i+1]-Px[i];
		dy = Py[i+1]-Py[i];
		d = sqrt(dx*dx + dy*dy)
		glLineWidth(1.0/d);
		glBegin(GL_LINES);
			glVertex3f(Px[i], 0.0, 0.0);
			glVertex3f(Px[i+1], 0.0, 0.0);
		glEnd();
	}

	// Drawing sodium wave fronts
	glColor3d(0.0,0.0,1.0);
	glPushMatrix();
	glTranslatef(APWaveFrontLx, APWaveFrontLy, 0.0);
	glutSolidSphere(0.005,20,20);
	glPopMatrix();
	
	glColor3d(0.0,0.0,1.0);
	glPushMatrix();
	glTranslatef(APWaveFrontRx, APWaveFrontRy, 0.0);
	glutSolidSphere(0.005,20,20);
	glPopMatrix();
	
	glutSwapBuffers();
}

void generalMuscleForces()
{
	float f; 
	float dx,dy,d;
	
	// Getting forces for the muscle fiber without contraction	
	for(int i = 0; i < N-1; i++)
	{
		dx = Px[i+1]-Px[i];
		dy = Py[i+1]-Py[i];
		d = sqrt(dx*dx + dy*dy)
		if(d < FiberCompresionStopFraction*FiberLength)
		{
			f  = FiberStrength*FiberCompresionMultiplier*(d - FiberLength);
		}
		else if(d < FiberLength)
		{
			f  = FiberStrength*(d - FiberLength);
		}
		else
		{
			f  = FiberStrength*FiberTentionMultiplier*(d - FiberLength);
		}
		
		Fx[i]   += f*dx/d;
		Fx[i+1] -= f*dx/d;
		Fy[i]   += f*dy/d;
		Fy[i+1] -= f*dy/d;
	}
}

int contractionForces(float dt, float time)
{
	float f; 
	float dx,dy,d;
	int flag;
	float ratio;
	int aPWaveFrontInMuscle;
	
	// Checking the sodium wave location.
	for(int i = 0; i < N-1; i++)
	{
		if(Px[i] <= APWaveFrontLx && APWaveFrontLx < Px[i+1])
		{
			if(Py[i] <= APWaveFrontLy && APWaveFrontLy < Py[i+1])
			{
				if((APWaveFrontLy-Py[i])*(Px[i+1]-Px[i]) - (APWaveFrontLx-Px[i])*(Py[i+1]-Py[i]))
				{
					ratio = (APWaveFront - Px[i])/(Px[i+1]-Px[i]);
					aPWaveFrontInMuscle = i;
					if(ContractionOn[i] == 0)
					{
						ContractionOn[i] = 1;
						APWaveAmunity = i;
					}
				}
			}
		}
	}
	
	// Getting forces for the muscle fiber contraction
	for(int i = 0; i < N-1; i++)
	{
		if(ContractionOn[i] == 1)
		{
			if(ContractionTime[i] < ContractionDuration[i])
			{
				Fx[i]   += ContractionStrength[i];
				Fx[i+1] -= ContractionStrength[i];
				ContractionTime[i] += dt;
			}
			else if(ContractionTime[i] < ContractionDuration[i] + RelaxationDuration[i])
			{
				Fx[i]   += 0.0;
				Fx[i+1] -= 0.0;
				ContractionTime[i] += dt;
			}
			else
			{
				ContractionOn[i] = 0;
				ContractionTime[i] = 0.0;
			}
		}
	}

	// Moving the system forward in time with leap-frog.
	for(int i = 0; i < N; i++)
	{
		if(time == 0.0)
		{
			Vx[i] += (Fx[i]/Mass[i])*0.5*dt;
		}
		else
		{
			Vx[i] += (Fx[i]/Mass[i])*dt;
		}

		Px[i] += Vx[i]*dt;
	}
	
	// Moving sodium wave front
	if(flag == -1)
	{
		APWaveFront = AttachmentLeft;
	}
	else
	{
		APWaveFront = Px[flag] + (Px[flag+1]-Px[flag])*ratio + APWaveSpeed[flag]*DT;
	}
}

int n_body()
{
	int   tdraw = 0; 
	double time = 0.0;
	float beatTimer = 0.0;
	
	for(int i = 0; i < N-1; i++)
	{
		ContractionOn[i] = 0;
		ContractionTime[i] = 0.0;
	}
	APWaveFrontLx = Px[0];
	APWaveFrontLy = Py[0];
	APWaveAmunity = -1;
	
	while(time < STOP_TIME)
	{
		// Zeroing out the nodal forces.
		for(int i = 0; i < N; i++)
		{
			Fx[i] = 0.0;
			Fy[i] = 0.0;
		}
	
		// Adding some damping to the system.
		for(int i = 0; i < N; i++)
		{	
			Fx[i]   += -Viscosity*Vx[i];
			Fy[i]   += -Viscosity*Vy[i];
		}
		
		generalMuscleForces();
		
		contractionForces(DT, time);

		if(tdraw == DRAW) 
		{
			draw_picture();
			tdraw = 0;
			printf("\n Time = %f", time);
		}
		else
		{
			tdraw++;
		}
		
		if(BeatPeriod < beatTimer)
		{
			APWaveFrontLx = Px[0];
			APWaveFrontLy = Py[0];
			APWaveFrontRx = Px[0];
			APWaveFrontRy = Py[0];
			beatTimer = 0.0;
		}
		else
		{
			beatTimer += DT;
		}
		
		time += DT;
	}
	return(1);
}

void control()
{	
	int    tdraw = 0;
	float  time = 0.0;

	glClear(GL_COLOR_BUFFER_BIT);
	glClear(GL_DEPTH_BUFFER_BIT);

	set_initial_conditions();
	
	draw_picture();
	
	if(n_body() == 1) printf("\n N-body success \n");
	
	printf("\n DONE \n");
	while(1);
}

void Display(void)
{
	gluLookAt(0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
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

	glFrustum(-0.2, 0.2, -0.2, 0.2, 0.2, 80.0);

	glMatrixMode(GL_MODELVIEW);
}

int main(int argc, char** argv)
{
	glutInit(&argc,argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGB);
	glutInitWindowSize(XWindowSize,YWindowSize);
	glutInitWindowPosition(0,0);
	glutCreateWindow("1D Myocardium");
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

