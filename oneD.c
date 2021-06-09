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
#define DT  0.0001
#define N 10
	
#define DRAW 1000

// Globals
float Px[N], Vx[N], Fx[N], Mass[N];

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

float SodiumWaveSpeed = 0.5; //0.01; //0.5;
float SodiumWaveFront;

float ContractionStrength = 5.0;
int ContractionOn[N-1];
float ContractionTime[N-1];
float ContractionStopTime = 100.0;
float ContractionRelaxationTime = 200.0;

float BeatPeriod = 200.0;

float Viscosity = 10.0;
float AttachmentLeft, AttachmentRight;

void set_initial_conditions()
{
	float centerX = FiberLength*(N+1)/2.0;
	AttachmentLeft = 0.0 - centerX;
	for(int i = 0; i < N; i++)
	{
		Mass[i] = 1.0;
		Px[i] = (float)(i+1)*FiberLength - centerX;
		Vx[i] = 0.0;	
	}
	AttachmentRight = (float)(N+1)*FiberLength - centerX;
}

void draw_picture()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glClear(GL_DEPTH_BUFFER_BIT);
	
	// Drawing nodes
	for(int i = 0; i < N; i++)
	{
		glColor3d(1.0,1.0,1.0);
		glPushMatrix();
		glTranslatef(Px[i], 0.0, 0.0);
		glutSolidSphere(0.005,20,20);
		glPopMatrix();	
	}
	
	// Drawing muscles
	glColor3d(1.0,0.0,0.0);
	for(int i = 0; i < N-1; i++)
	{
		glLineWidth(1.0/(Px[i+1]-Px[i]));
		glBegin(GL_LINES);
			glVertex3f(Px[i], 0.0, 0.0);
			glVertex3f(Px[i+1], 0.0, 0.0);
		glEnd();
	}
	
	// Drawing Tendons
	glColor3d(1.0,1.0,1.0);
	glLineWidth(1.0/(Px[0]-AttachmentLeft));
	glBegin(GL_LINES);
		glVertex3f(AttachmentLeft, 0.0, 0.0);
		glVertex3f(Px[0], 0.0, 0.0);
	glEnd();
	
	glLineWidth(1.0/(AttachmentRight - Px[N-1]));
	glBegin(GL_LINES);
		glVertex3f(Px[N-1], 0.0, 0.0);
		glVertex3f(AttachmentRight, 0.0, 0.0);
	glEnd();

	// Drawing sodium wave front
	glColor3d(1.0,1.0,0.0);
	glLineWidth(2.0);
	glBegin(GL_LINES);
		glVertex3f(SodiumWaveFront, -0.5, 0.0);
		glVertex3f(SodiumWaveFront, 0.5, 0.0);
	glEnd();
	
	glutSwapBuffers();
}

void generalMuscleForces()
{
	float f; 
	float dx,d;
	
	// Getting forces for the muscle fiber without contraction	
	for(int i = 0; i < N-1; i++)
	{
		dx = Px[i+1]-Px[i];
		d  = sqrt(dx*dx);
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
	}
	
	dx = AttachmentLeft-Px[0];
	d  = sqrt(dx*dx);
	if(d < TendonCompresionStopFraction*TendonLength)
	{
		f  = TendonStrength*TendonCompresionMultiplier*(d - TendonLength);
	}
	else if(d < FiberLength)
	{
		f  = TendonStrength*(d - FiberLength);
	}
	else
	{
		f  = TendonStrength*TendonTentionMultiplier*(d - TendonLength);
	}
	Fx[0]   += f*dx/d;
	
	dx = AttachmentRight-Px[N-1];
	d  = sqrt(dx*dx);
	if(d < TendonCompresionStopFraction*FiberLength)
	{
		f  = TendonStrength*TendonCompresionMultiplier*(d - TendonLength);
	}
	else if(d < FiberLength)
	{
		f  = TendonStrength*(d - TendonLength);
	}
	else
	{
		f  = TendonStrength*TendonTentionMultiplier*(d - TendonLength);
	}
	Fx[N-1]   += f*dx/d;
}

int contractionForces(float dt, float time)
{
	float f; 
	float dx,d;
	int flag;
	float ratio;
	
	// Checking the sodium wave location.
	for(int i = 0; i < N-1; i++)
	{
		if(Px[i] <= SodiumWaveFront && SodiumWaveFront < Px[i+1])
		{
			if(ContractionOn[i] == 0)
			{
				ContractionOn[i] = 1;
			}
		}
	}
	
	// Getting forces for the muscle fiber contraction
	for(int i = 0; i < N-1; i++)
	{
		if(ContractionOn[i] == 1)
		{
			if(ContractionTime[i] < ContractionStopTime)
			{
				Fx[i]   += ContractionStrength;
				Fx[i+1] -= ContractionStrength;
				ContractionTime[i] += dt;
			}
			else if(ContractionTime[i] < ContractionStopTime + ContractionRelaxationTime)
			{
				Fx[i]   += 0.0;
				Fx[i+1] -= 0.0;
				ContractionTime[i] += dt;
			}
			else
			{
				ContractionOn[i] == 0;
				ContractionTime[i] = 0.0;
			}
		}
	}
	
	// Fiding which cell the sodium wave front is in. Then finding it's ratio between them.
	flag = -1;
	for(int i = 0; i < N-1; i++)
	{
		if(Px[i] <= SodiumWaveFront && SodiumWaveFront < Px[i+1])
		{
			ratio = (SodiumWaveFront - Px[i])/(Px[i+1]-Px[i]);
			flag = i;
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
		SodiumWaveFront = AttachmentLeft;
	}
	else
	{
		SodiumWaveFront = Px[flag] + (Px[flag+1]-Px[flag])*ratio + SodiumWaveSpeed*DT;
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
	SodiumWaveFront = Px[0];
	
	while(time < STOP_TIME)
	{
		// Zeroing out the nodal forces.
		for(int i = 0; i < N; i++)
		{
			Fx[i] = 0.0;
		}
	
		// Adding some damping to the system.
		for(int i = 0; i < N; i++)
		{	
			Fx[i]   += -Viscosity*Vx[i];
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
			SodiumWaveFront = Px[0];
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

	glFrustum(-0.2, 0.2, -0.2, 0.2, 0.2, 50.0);

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

