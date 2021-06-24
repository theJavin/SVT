// nvcc SVT.cu -o svt -lglut -lm -lGLU -lGL
//To stop hit "control c" in the window you launched it from.

// Length will be in millimeters
// Time will be in milliseconds
// Mass will be in... mass units?

// Fiber length 100 micrometers or 0.1 millimeters
// Sodium wave speed .5 meters/sec or 0.5 millimeters/millisec

#include <GL/glut.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define PI 3.141592654

#define STOP_TIME 60000.0
#define DT  0.001

#define NUMBER_OF_NODES 266 //266 //62
#define NUMBER_OF_MUSCLES 552 //552 //132
#define LINKS_PER_NODE 24 //24

// Globals
int DrawRate;

float4 NodePosition[NUMBER_OF_NODES], NodeVelocity[NUMBER_OF_NODES], NodeForce[NUMBER_OF_NODES];
float NodeMass[NUMBER_OF_NODES];
int NodeLinks[NUMBER_OF_NODES][LINKS_PER_NODE]; // The nodes that this node is connected to
int NodeMuscles[NUMBER_OF_NODES][LINKS_PER_NODE]; // The muscle that connects this node to ther other nodes
int NodeAblatedYesNo[NUMBER_OF_NODES];

// How the muscle will act without contraction.
int MuscleConectionA[NUMBER_OF_MUSCLES];
int MuscleConectionB[NUMBER_OF_MUSCLES];
float MuscleMass[NUMBER_OF_MUSCLES];
float MuscleLength[NUMBER_OF_MUSCLES];
float MuscleRelaxedStrength[NUMBER_OF_MUSCLES];
float MuscleCompresionMultiplier;
float MuscleTentionMultiplier;
float MuscleCompresionStopFraction[NUMBER_OF_MUSCLES];  // 0.6 is the standard value
float Viscosity;
float3 MuscleColor[NUMBER_OF_MUSCLES];

// Muscle contraction parameters
int ContractionOnOff[NUMBER_OF_MUSCLES];
float ContractionTimer[NUMBER_OF_MUSCLES];
float ActionPotentialSpeed[NUMBER_OF_MUSCLES]; //0.5 is a good value.
float ActionPotentialDuration[NUMBER_OF_MUSCLES];
float ContractionDuration[NUMBER_OF_MUSCLES]; // 100.0 is a good value
float RelaxationDuration[NUMBER_OF_MUSCLES]; // 200.0 is a good value
float ContractionStrength[NUMBER_OF_MUSCLES]; // 5.0 is a good value

float BloodPresure;
float BeatPeriod;

// Prototyping functions
void initializeNodesAndLinksSphere62();
void linkNodesToMuscles();
void setMuscleAtributesAndNodeMasses();
void ablatedNodes();
void draw_picture();
void generalMuscleForces();
void outwardPresure();
void turnOnNodeMuscles(int);
int contractionForces(float, float);
void dampingForce();
void moveNodes(float, float);
void ectopicEvents(float, float);
int n_body(float);
void control();
void mymouse(int, int, int, int);
void Display(void);
void reshape(int, int);

#include "./setNodesAndLinks.h"

void linkNodesToMuscles()
{
	//Setting the ends of the muscles
	int index = 0;
	for(int i = 0; i < NUMBER_OF_NODES; i++)
	{
		for(int j = 0; j < LINKS_PER_NODE; j++)
		{
			if(NodeLinks[i][j] != -1)
			{
				if(i < NodeLinks[i][j])
				{
					if(NUMBER_OF_MUSCLES <= index)
					{
						printf("\nTSU Error: number of muscles is out of bounds\n");
						exit(0);
					} 
					MuscleConectionA[index] = i;
					MuscleConectionB[index] = NodeLinks[i][j];
					index++;
				}
			}
		}
	}
	
	// Setting the node muscles. Each node will have a list of nodes they are attached to (NodeLinks[][]) and the muscle that attaches it to that node (NodeMuscles[][]).
	// Setting them all to -1 first.
	for(int i = 0; i < NUMBER_OF_NODES; i++)
	{
		for(int j = 0; j < LINKS_PER_NODE; j++)
		{
			NodeMuscles[i][j] = -1;
		}	
	}
	for(int i = 0; i < NUMBER_OF_NODES; i++)
	{
		for(int j = 0; j < LINKS_PER_NODE; j++)
		{
			if(NodeLinks[i][j] != -1)
			{
				for(int k = 0; k < NUMBER_OF_MUSCLES; k++)
				{
					if((MuscleConectionA[k] == i && MuscleConectionB[k] == NodeLinks[i][j]) || (MuscleConectionA[k] == NodeLinks[i][j] && MuscleConectionB[k] == i))
					{
						NodeMuscles[i][j] = k;
					}
				}
			}
		}
	}
}

void setMuscleAtributesAndNodeMasses()
{	
	float dx, dy, dz;
	float sum;
	
	MuscleCompresionMultiplier = 50.0;
	MuscleTentionMultiplier = 50.0;
	Viscosity = 5.0;
	BloodPresure = 0.02;
	
	// Setting other parameters
	for(int i = 0; i < NUMBER_OF_MUSCLES; i++)
	{	
		MuscleMass[i] = 1.0;
		dx = NodePosition[MuscleConectionA[i]].x - NodePosition[MuscleConectionB[i]].x;
		dy = NodePosition[MuscleConectionA[i]].y - NodePosition[MuscleConectionB[i]].y;
		dz = NodePosition[MuscleConectionA[i]].z - NodePosition[MuscleConectionB[i]].z;
		MuscleLength[i] = sqrt(dx*dx + dy*dy + dz*dz);;
		MuscleRelaxedStrength[i] = 0.1;
		MuscleCompresionStopFraction[i] = 0.6;
		ContractionOnOff[i] = 0;
		ContractionTimer[i] = 0.0;
		ActionPotentialSpeed[i] = 0.02; // 0.2
		ActionPotentialDuration[i] = MuscleLength[i]/ActionPotentialSpeed[i];
		ContractionDuration[i] = 20.0;  // 100.0
		RelaxationDuration[i] = 60.0;  // 200.0
		ContractionStrength[i] = 0.2; //0.1;
		MuscleColor[i].x = 1.0;
		MuscleColor[i].y = 0.0;
		MuscleColor[i].z = 0.0;
	}
	
	// Setting the node masses
	for(int i = 0; i < NUMBER_OF_NODES; i++)
	{
		sum = 0.0;
		for(int j = 0; j < LINKS_PER_NODE; j++)
		{
			if(NodeMuscles[i][j] != -1)
			{
				sum += MuscleMass[NodeMuscles[i][j]];
			}
		}
		NodeMass[i] = sum/2.0;
	}
}

void draw_picture()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glClear(GL_DEPTH_BUFFER_BIT);
	
	glColor3d(1.0,1.0,0.0);
	glPushMatrix();
	glTranslatef(NodePosition[0].x, NodePosition[0].y, NodePosition[0].z);
	glutSolidSphere(0.03,20,20);
	glPopMatrix();
	
	// Drawing nodes
	for(int i = 1; i < NUMBER_OF_NODES; i++)
	{
		if(NodeAblatedYesNo[i] == 0)
		{
			glColor3d(0.0,1.0,0.0);
		}
		else
		{
			glColor3d(1.0,1.0,1.0);
		}
		glPushMatrix();
		glTranslatef(NodePosition[i].x, NodePosition[i].y, NodePosition[i].z);
		glutSolidSphere(0.01,20,20);
		glPopMatrix();	
	}
	
	/*
	glColor3d(1.0,1.0,0.0);
	glPushMatrix();
	glTranslatef(NodePosition[NUMBER_OF_NODES-1].x, NodePosition[NUMBER_OF_NODES-1].y, NodePosition[NUMBER_OF_NODES-1].z);
	glutSolidSphere(0.03,20,20);
	glPopMatrix();
	*/
	
	// Drawing muscles
	glColor3d(1.0,0.0,0.0);
	for(int i = 0; i < NUMBER_OF_NODES; i++)
	{
		for(int j = 0; j < LINKS_PER_NODE; j++)
		{
			if(NodeLinks[i][j] != -1)
			{
				glColor3d(MuscleColor[NodeMuscles[i][j]].x, MuscleColor[NodeMuscles[i][j]].y, MuscleColor[NodeMuscles[i][j]].z);
				//glLineWidth(1.0/(Px[i]-Px[NodeLinks[i][j]]));
				glBegin(GL_LINES);
					glVertex3f(NodePosition[i].x, NodePosition[i].y, NodePosition[i].z);
					glVertex3f(NodePosition[NodeLinks[i][j]].x, NodePosition[NodeLinks[i][j]].y, NodePosition[NodeLinks[i][j]].z);
				glEnd();
			}
			
		}	
	}
	
	glutSwapBuffers();
}

void generalMuscleForces()
{
	float f; 
	float dx, dy, dz, d;
	int muscleNumber, nodeNumber;
	
	// Getting forces on the nodes from the muscle fiber without contraction	
	for(int i = 0; i < NUMBER_OF_NODES; i++)
	{
		for(int j = 0; j < LINKS_PER_NODE; j++)
		{
			muscleNumber = NodeMuscles[i][j];
			nodeNumber = NodeLinks[i][j];
			if(nodeNumber != -1)
			{
				dx = NodePosition[nodeNumber].x - NodePosition[i].x;
				dy = NodePosition[nodeNumber].y - NodePosition[i].y;
				dz = NodePosition[nodeNumber].z - NodePosition[i].z;
				d  = sqrt(dx*dx + dy*dy + dz*dz);
				if(d < 0.00001) 
				{
					printf("\n TSU Error: In generalMuscleForces d is very small between nodeNumbers = %d %d seperation = %f. Take a look at this\n", i, nodeNumber, d);
					glColor3d(0.0,0.0,1.0);
					glPushMatrix();
					glTranslatef(NodePosition[i].x, NodePosition[i].y, NodePosition[i].z);
					glutSolidSphere(0.03,20,20);
					glPopMatrix();
					glPushMatrix();
					glTranslatef(NodePosition[nodeNumber].x, NodePosition[nodeNumber].y, NodePosition[nodeNumber].z);
					glutSolidSphere(0.03,20,20);
					glPopMatrix();
					glutSwapBuffers();
					while(1);
				}
				if(d < MuscleCompresionStopFraction[muscleNumber]*MuscleLength[muscleNumber])
				{
					f  = MuscleRelaxedStrength[muscleNumber]*MuscleCompresionMultiplier*(d - MuscleLength[muscleNumber]);
				}
				else if(d < MuscleLength[muscleNumber])
				{
					f  = MuscleRelaxedStrength[muscleNumber]*(d - MuscleLength[muscleNumber]);
				}
				else
				{
					f  = MuscleRelaxedStrength[muscleNumber]*MuscleTentionMultiplier*(d - MuscleLength[muscleNumber]);
				}
				NodeForce[i].x  += f*dx/d;
				NodeForce[i].y  += f*dy/d;
				NodeForce[i].z  += f*dz/d;
			}
		}
	}
}

void outwardPresure()
{
	float f; 
	float dx, dy, dz, d;
	float4 centerOfMass;
	
	for(int i = 0; i < NUMBER_OF_NODES; i++)
	{
		 centerOfMass.x += NodePosition[i].x*NodeMass[i];
		 centerOfMass.y += NodePosition[i].y*NodeMass[i];
		 centerOfMass.z += NodePosition[i].z*NodeMass[i];
		 centerOfMass.w += NodeMass[i];
	}
	
	centerOfMass.x /= centerOfMass.w;
	centerOfMass.y /= centerOfMass.w;
	centerOfMass.z /= centerOfMass.w;
		 
	for(int i = 0; i < NUMBER_OF_NODES; i++)
	{
		 NodePosition[i].x -= centerOfMass.x;
		 NodePosition[i].y -= centerOfMass.y;
		 NodePosition[i].z -= centerOfMass.z;
	}
	
	// Getting forces on the nodes from the presure of the blood pushing out	
	for(int i = 0; i < NUMBER_OF_NODES; i++)
	{
		dx = 0.0 - NodePosition[i].x;
		dy = 0.0 - NodePosition[i].y;
		dz = 0.0 - NodePosition[i].z;
		d  = sqrt(dx*dx + dy*dy + dz*dz);
		if(d < 0.0001) 
		{
			printf("\nTSU Error: In outwardPresure d is very small. Take a look at this\n");
			exit(0);
		}
		
		f  = -BloodPresure;
		
		NodeForce[i].x  += f*dx/d;
		NodeForce[i].y  += f*dy/d;
		NodeForce[i].z  += f*dz/d;
	}
}

void turnOnNodeMuscles(int index)
{
	for(int j = 0; j < LINKS_PER_NODE; j++)
	{
		if((NodeLinks[index][j] != -1) && (ContractionOnOff[NodeMuscles[index][j]] == 0))
		{
			ContractionOnOff[NodeMuscles[index][j]] = 1;
			ContractionTimer[NodeMuscles[index][j]] = 0.0;
		}
	}
}

int contractionForces(float dt, float time)
{
	float dx, dy, dz, d;
	int muscleNumber, nodeNumber;
	
	// Getting forces for the muscle fiber contraction
	for(int i = 0; i < NUMBER_OF_NODES; i++) 
	{
		for(int j = 0; j < LINKS_PER_NODE; j++)
		{
			muscleNumber = NodeMuscles[i][j];
			nodeNumber = NodeLinks[i][j];
			if(nodeNumber != -1)
			{	
				if(ContractionOnOff[muscleNumber] == 1)
				{
					if((ActionPotentialDuration[muscleNumber] - dt < ContractionTimer[muscleNumber]) && (ContractionTimer[muscleNumber] < ActionPotentialDuration[muscleNumber] + dt))
					{
						if(NodeAblatedYesNo[nodeNumber] == 0)
						{
							turnOnNodeMuscles(nodeNumber);
						}
					}
					
					if(ContractionTimer[muscleNumber] < ActionPotentialDuration[muscleNumber])
					{
						MuscleColor[muscleNumber].x = 1.0;
						MuscleColor[muscleNumber].y = 1.0;
						MuscleColor[muscleNumber].z = 1.0;
					}
					else
					{
						MuscleColor[muscleNumber].x = 1.0;
						MuscleColor[muscleNumber].y = 0.0;
						MuscleColor[muscleNumber].z = 0.0;
					}
					
					if(ContractionTimer[muscleNumber] < ContractionDuration[muscleNumber])
					{
						dx = NodePosition[nodeNumber].x - NodePosition[i].x;
						dy = NodePosition[nodeNumber].y - NodePosition[i].y;
						dz = NodePosition[nodeNumber].z - NodePosition[i].z;
						d  = sqrt(dx*dx + dy*dy + dz*dz);
						if(d < 0.00001) 
						{
							printf("\n TSU Error: In contractionForces d is very small between nodeNumbers = %d %d seperation = %f. Take a look at this\n", i, nodeNumber, d);
							glColor3d(0.0,0.0,1.0);
							glPushMatrix();
							glTranslatef(NodePosition[i].x, NodePosition[i].y, NodePosition[i].z);
							glutSolidSphere(0.03,20,20);
							glPopMatrix();
							glPushMatrix();
							glTranslatef(NodePosition[nodeNumber].x, NodePosition[nodeNumber].y, NodePosition[nodeNumber].z);
							glutSolidSphere(0.03,20,20);
							glPopMatrix();
							glutSwapBuffers();
							while(1);
						}
						
						NodeForce[i].x   += ContractionStrength[muscleNumber]*dx/d;
						NodeForce[i].y   += ContractionStrength[muscleNumber]*dy/d;
						NodeForce[i].z   += ContractionStrength[muscleNumber]*dz/d;
					
						ContractionTimer[muscleNumber] += dt;
					}
					else if(ContractionTimer[muscleNumber] < ContractionDuration[muscleNumber] + RelaxationDuration[muscleNumber])
					{
						NodeForce[i].x   += 0.0;
						NodeForce[i].y   += 0.0;
						NodeForce[i].z   += 0.0;
						
						ContractionTimer[muscleNumber] += dt;
					}
					else
					{
						ContractionOnOff[muscleNumber] = 0;
						ContractionTimer[muscleNumber] = 0.0;
					}
				}
			}
		}
	}
	return(1);
}

void dampingForce()
{
	for(int i = 0; i < NUMBER_OF_NODES; i++)
	{	
		NodeForce[i].x   += -Viscosity*NodeVelocity[i].x;
		NodeForce[i].y   += -Viscosity*NodeVelocity[i].y;
		NodeForce[i].z   += -Viscosity*NodeVelocity[i].z;
	}
}

void moveNodes(float dt, float time)  // LeapFrog
{
	// Moving the system forward in time with leap-frog.
	for(int i = 0; i < NUMBER_OF_NODES; i++)
	{
		if(time == 0.0)
		{
			NodeVelocity[i].x += (NodeForce[i].x/MuscleMass[i])*0.5*dt;
			NodeVelocity[i].y += (NodeForce[i].y/MuscleMass[i])*0.5*dt;
			NodeVelocity[i].z += (NodeForce[i].z/MuscleMass[i])*0.5*dt;
		}
		else
		{
			NodeVelocity[i].x += (NodeForce[i].x/MuscleMass[i])*dt;
			NodeVelocity[i].y += (NodeForce[i].y/MuscleMass[i])*dt;
			NodeVelocity[i].z += (NodeForce[i].z/MuscleMass[i])*dt;
		}

		NodePosition[i].x += NodeVelocity[i].x*dt;
		NodePosition[i].y += NodeVelocity[i].y*dt;
		NodePosition[i].z += NodeVelocity[i].z*dt;
	}
}

int n_body(float dt)
{
	int   tdraw = 0; 
	double time = 0.0;
	float beatTimer = 0.0;
	
	while(time < STOP_TIME)
	{
		if(BeatPeriod <= beatTimer)
		{
			turnOnNodeMuscles(0);
			beatTimer = 0.0;
		}
		else beatTimer += dt;
		
		ectopicEvents(time, dt);
		
		// Zeroing out the nodal forces.
		for(int i = 0; i < NUMBER_OF_NODES; i++)
		{
			NodeForce[i].x   = 0.0;
			NodeForce[i].y   = 0.0;
			NodeForce[i].z   = 0.0;
		}
		
		generalMuscleForces();
		
		contractionForces(dt, time);
		
		outwardPresure();
		
		dampingForce();
		
		moveNodes(dt, time);

		if(tdraw == DrawRate) 
		{
			draw_picture();
			tdraw = 0;
			printf("\n Time = %f", time);
		}
		else tdraw++;
		
		time += dt;
	}
	return(1);
}

void ablatedNodes()
{
	// Setting all nodes as not ablated
	for(int i = 0; i < NUMBER_OF_NODES; i++)
	{
		NodeAblatedYesNo[i] = 0;
	}
	
	//Nodes to ablate
	for(int i = 1; i < 23; i++)
	{	
		//NodeAblatedYesNo[i] = 1;
	}
	
	
	NodeAblatedYesNo[1] = 1;
	NodeAblatedYesNo[13] = 1;
	NodeAblatedYesNo[25] = 1;
	NodeAblatedYesNo[37] = 1;
	NodeAblatedYesNo[49] = 1;
	
	NodeAblatedYesNo[2] = 1;
	NodeAblatedYesNo[3] = 1;
	NodeAblatedYesNo[4] = 1;
	NodeAblatedYesNo[5] = 1;
	NodeAblatedYesNo[6] = 1;
	NodeAblatedYesNo[7] = 1;
	NodeAblatedYesNo[8] = 1;
	
	NodeAblatedYesNo[9] = 1;
	NodeAblatedYesNo[10] = 1;
	NodeAblatedYesNo[11] = 1;
	NodeAblatedYesNo[14] = 1;
	NodeAblatedYesNo[17] = 1;
	
}

void ectopicEvents(float time, float dt)
{
	float er = dt/2.0;
	
	if((50.0 - er <= time) && (time < 220.0 + er))
	{
		turnOnNodeMuscles(31);
	}
	
	if((51.0 - er <= time) && (time < 230.0 + er))
	{
		//turnOnNodeMuscles(41);
	}
	
	if((240.0 - er <= time) && (time < 230.0 + er))
	{
		//turnOnNodeMuscles(59);
	}
}

void control()
{	
	
	//initializeNodesAndLinksSphere62();
	//initializeNodesAndLinksSphere266();
	initializeNodesAndLinksSphere(24);
	
	linkNodesToMuscles();
	setMuscleAtributesAndNodeMasses();
	ablatedNodes();
	
	draw_picture();
	
	DrawRate = 1000;
	BeatPeriod = 100;
	
	if(n_body(DT) == 1) printf("\n N-body success \n");
	
	printf("\n DONE \n");
	while(1);
}

// Window globals
int XWindowSize = 1000;
int YWindowSize = 1000; 

// Clip plains
double Near = 0.2;
double Far = 80.0;

//Direction here your eye is located location
double EyeX = 0.0;
double EyeY = 2.0;
double EyeZ = 2.0;

//Where you are looking
double CenterX = 0.0;
double CenterY = 0.0;
double CenterZ = 0.0;

//Up vector for viewing
double UpX = 0.0;
double UpY = 1.0;
double UpZ = 0.0;

void Display(void)
{
	gluLookAt(EyeX, EyeY, EyeZ, CenterX, CenterY, CenterZ, UpX, UpY, UpZ);
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

	glFrustum(-0.2, 0.2, -0.2, 0.2, Near, Far);

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

