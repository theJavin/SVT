// gcc oneD1.1.c -o temp -lglut -lm -lGLU -lGL
//To stop hit "control c" in the window you launched it from. stuff
// Length will be in millimeters
// Time will be in seconds
// Mass will be in ???
// Fiber length 100 micrometers or 0.1 millimeters
// Sodium wave speed .5 meters/sec or 500 millimeters/sec
#include <GL/glut.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define PI 3.141592654

#define XWindowSize 2000
#define YWindowSize 2000 

#define STOP_TIME 60.0
#define DT  0.00001
#define N 10
	
#define DRAW 100

// Globals
float Px[N], Vx[N], Fx[N], Mass[N];
float FiberLength = 0.1;
float FiberStrength = 50.0;
float SodiumWaveSpeed = 500.0;
float SodiumWaveFront;
int ContractionOn[N-1];
float ContractionTime[N-1];
float ContractionStrength = 20.0;
float ContractionStopTime = 0.5;
float Viscosity = 10.0;
float BeatPeriod = 1.0;
float AttachmentStart, AttachmentStop;

void set_initial_conditions()
{
	float centerX = FiberLength*(N+1)/2.0;
	AttachmentStart = 0.0 - centerX;
	for(int i = 0; i < N; i++)
	{
		Mass[i] = 1.0;
		Px[i] = (float)(i+1)*FiberLength - centerX;
		Vx[i] = 0.0;	
	}
	AttachmentStop = (float)(N+1)*FiberLength - centerX;
}

void draw_picture()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glClear(GL_DEPTH_BUFFER_BIT);
	
	glColor3d(1.0,1.0,0.5);
	glPushMatrix();
	glTranslatef(Px[0], 0.0, 0.0);
	glutSolidSphere(0.005,20,20);
	glPopMatrix();
	
	for(int i = 1; i < N; i++)
	{
		glColor3d(1.0,1.0,1.0);
		glPushMatrix();
		glTranslatef(Px[i], 0.0, 0.0);
		glutSolidSphere(0.005,20,20);
		glPopMatrix();	
	}
	
	glColor3d(1.0,0.0,0.0);
	for(int i = 0; i < N-1; i++)
	{
		glLineWidth(1.0/(Px[i+1]-Px[i]));
		glBegin(GL_LINES);
			glVertex3f(Px[i], 0.0, 0.0);
			glVertex3f(Px[i+1], 0.0, 0.0);
		glEnd();
	}
	
	glColor3d(1.0,1.0,1.0);
	//glLineWidth(6.0);
	glLineWidth(1.0/(Px[0]-AttachmentStart));
	glBegin(GL_LINES);
		glVertex3f(AttachmentStart, 0.0, 0.0);
		glVertex3f(Px[0], 0.0, 0.0);
	glEnd();
	
	glLineWidth(1.0/(AttachmentStop - Px[N-1]));
	glBegin(GL_LINES);
		glVertex3f(Px[N-1], 0.0, 0.0);
		glVertex3f(AttachmentStop, 0.0, 0.0);
	glEnd();

	glColor3d(1.0,1.0,0.0);
	glLineWidth(2.0);
	glBegin(GL_LINES);
		glVertex3f(SodiumWaveFront, -0.5, 0.0);
		glVertex3f(SodiumWaveFront, 0.5, 0.0);
	glEnd();
	
	glutSwapBuffers();
}

int leapFrog(float dt, float time, float SodiumWaveFront)
{
	float f; 
	float dx,d;
	
	// Zeroing out the nodal forces.
	for(int i = 0; i < N; i++)
	{
		Fx[i] = 0.0;
	}
	
	// Getting forces for the muscle fiber without contraction	
	for(int i = 0; i < N-1; i++)
	{
		dx = Px[i+1]-Px[i];
		d  = sqrt(dx*dx);
		if(d < 0.5*FiberLength)
		{
			f  = FiberStrength*10.0*(d - FiberLength);
		}
		else if(d < FiberLength)
		{
			f  = FiberStrength*(d - FiberLength);
		}
		else
		{
			f  = FiberStrength*10.0*(d - FiberLength);
		}
		
		Fx[i]   += f*dx/d;
		Fx[i+1] -= f*dx/d;
	}
	
	dx = AttachmentStart-Px[0];
	d  = sqrt(dx*dx);
	if(d < 0.5*FiberLength)
	{
		f  = FiberStrength*10.0*(d - FiberLength);
	}
	else if(d < FiberLength)
	{
		f  = 0.0*(d - FiberLength);
	}
	else
	{
		f  = FiberStrength*1.0*(d - FiberLength);
	}
	Fx[0]   += f*dx/d;
	
	dx = AttachmentStop-Px[N-1];
	d  = sqrt(dx*dx);
	if(d < 0.5*FiberLength)
	{
		f  = FiberStrength*10.0*(d - FiberLength);
	}
	else if(d < FiberLength)
	{
		f  = 0.0*(d - FiberLength);
	}
	else
	{
		f  = FiberStrength*1.0*(d - FiberLength);
	}
	Fx[N-1]   += f*dx/d;
	
	// 
	for(int i = 0; i < N-1; i++)
	{
		if(Px[i] < SodiumWaveFront)
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
		if(ContractionOn[i] == 1 && ContractionTime[i] < ContractionStopTime)
		{
			Fx[i]   += ContractionStrength;
			Fx[i+1] -= ContractionStrength;
			ContractionTime[i] += dt;
		}
	}
	
	// Adding some damping to the system.
	for(int i = 0; i < N; i++)
	{	
		Fx[i]   += -Viscosity*Vx[i];
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
}

int n_body()
{
	int    tdraw = 0; 
	float  time = 0.0;
	float beatTimer = 0.0;
	
	for(int i = 0; i < N-1; i++)
	{
		ContractionOn[i] = 0;
		ContractionTime[i] = 0.0;
	}
	SodiumWaveFront = Px[0];
	
	while(time < STOP_TIME)
	{
		leapFrog(DT, time, SodiumWaveFront);

		if(tdraw == DRAW) 
		{
			draw_picture();
			tdraw = 0;
		}
		
		SodiumWaveFront += SodiumWaveSpeed*DT;
		time += DT;
		tdraw++;
		
		if(BeatPeriod < beatTimer)
		{
			for(int i = 0; i < N-1; i++)
			{
				ContractionOn[i] = 0;
				ContractionTime[i] = 0.0;
			}
			SodiumWaveFront = Px[0];
			beatTimer = 0.0;
		}
		beatTimer += DT;
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

