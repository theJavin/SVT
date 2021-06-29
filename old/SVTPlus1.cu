// nvcc SVTPlus.cu -o svtplus -lglut -lm -lGLU -lGL -lmenu -lncurses
//To stop hit "control c" in the window you launched it from.

// Length will be in millimeters
// Time will be in milliseconds
// Mass will be in ???

// Fiber length 100 micrometers or 0.1 millimeters
// Sodium wave speed .5 meters/sec or 0.5 millimeters/millisec

/*
#include <GL/glut.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <curses.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
using namespace std;
*/

#include <GL/glut.h>
#include <GL/glu.h>
#include <GL/gl.h>
#include <math.h>
#include <stdio.h>
#include "stdio.h"
#include <stdlib.h>
#include <cuda.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <curand.h>
#include <curand_kernel.h>
#include <signal.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

#define PI 3.141592654

#define STOP_TIME 60000.0
#define DT  0.001

#define NUMBER_OF_NODES 62
#define NUMBER_OF_MUSCLES 132

// Globals
float4 NodePosition[NUMBER_OF_NODES], NodeVelocity[NUMBER_OF_NODES], NodeForce[NUMBER_OF_NODES];
float NodeMass[NUMBER_OF_NODES];
int NodeLinks[NUMBER_OF_NODES][12]; // The nodes that this node is connected to
int NodeMuscles[NUMBER_OF_NODES][12]; // The muscle that connects this node to ther other nodes
int NodeAblatedYesNo[NUMBER_OF_NODES];

// How the muscle will act without contraction.
int MuscleConectionA[NUMBER_OF_MUSCLES];
int MuscleConectionB[NUMBER_OF_MUSCLES];
float MuscleMass[NUMBER_OF_MUSCLES];
float MuscleLength[NUMBER_OF_MUSCLES];
float MuscleRelaxedStrength[NUMBER_OF_MUSCLES];
float MuscleCompresionMultiplier = 10.0;
float MuscleTentionMultiplier = 10.0;
float MuscleCompresionStopFraction[NUMBER_OF_MUSCLES];  // 0.6 is the standard value
float Viscosity = 10.0;
float3 MuscleColor[NUMBER_OF_MUSCLES];

// Muscle contraction parameters
int ContractionOnOff[NUMBER_OF_MUSCLES];
float ContractionTimer[NUMBER_OF_MUSCLES];
float ActionPotentialSpeed[NUMBER_OF_MUSCLES]; //0.5 is a good value.
float ActionPotentialDuration[NUMBER_OF_MUSCLES];
float ContractionDuration[NUMBER_OF_MUSCLES]; // 100.0 is a good value
float RelaxationDuration[NUMBER_OF_MUSCLES]; // 200.0 is a good value
float ContractionStrength[NUMBER_OF_MUSCLES]; // 5.0 is a good value

float BloodPresure = 0.05;

float BeatPeriod = 400.0;

int DrawRate = 1000;

//Globals for setting up the viewing window 
int XWindowSize = 1000;
int YWindowSize = 1000; 
double Near = 0.2;
double Far = 80.0;

/*
double ViewBoxSize = 300.0;

GLdouble Left = -ViewBoxSize;
GLdouble Right = ViewBoxSize;
GLdouble Bottom = -ViewBoxSize;
GLdouble Top = ViewBoxSize;
GLdouble Front = ViewBoxSize;
GLdouble Back = -ViewBoxSize;
*/

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

// Prototyping functions
void draw_picture();
void generalMuscleForces();
void outwardPresure();
void turnOnNodeMuscles(int);
int contractionForces(float, float);
void dampingForce();
void moveNodes(float, float);
int n_body();
void control();
void mymouse(int, int, int, int);
void Display(void);
void reshape(int, int);
static void signalHandler(int);

int set_initial_conditions()
{	
	int index;
	float dx, dy, dz;
	float sum;
	
	// Node position values for a sphere with 62 nodes//0.5 is a good value.
	NodePosition[0].x = 0.0;
	NodePosition[0].y = 1.0;
	NodePosition[0].z = 0.0;
	NodePosition[NUMBER_OF_NODES-1].x = 0.0;
	NodePosition[NUMBER_OF_NODES-1].y = -1.0;
	NodePosition[NUMBER_OF_NODES-1].z = 0.0;
	
	index = 1;
	for(int i = 1; i < 6; i++)
	{
		for(int j = 0; j < 12; j++)
		{
			if((NUMBER_OF_NODES-1) <= index)
			{
				printf("\nTSU Error: number of nodes is out of bounds\n");
				return(0);
			} 
			NodePosition[index].y = sin(PI/2.0 -i*PI/6.0);
			NodePosition[index].x = cos(PI/2.0 -i*PI/6.0)*cos(j*PI/6.0);
			NodePosition[index].z = cos(PI/2.0 -i*PI/6.0)*sin(j*PI/6.0);
			
			index++;
		}	
	}
	
	// Zeroing out velocity and acceleration
	for(int i = 0; i < NUMBER_OF_NODES; i++)
	{
		NodeVelocity[index].y = 0.0;
		NodeVelocity[index].x = 0.0;
		NodeVelocity[index].z = 0.0;
		
		NodeForce[index].y = 0.0;
		NodeForce[index].x = 0.0;
		NodeForce[index].z = 0.0;
	}
	

	// Below are the edges for the links connecting 62 node sphere.
	// 0: 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12
	
	// 1:  0  12 2  13		13: 1  24 14 25		25: 13 36 26 37		37: 25 48 38 49		49: 37 60 50 61
	// 2:  0  1  3  14		14: 2  13 15 26		26: 14 25 27 38		38: 26 37 39 50		50: 38 49 51 61
	// 3:  0  2  4  15		15: 3  14 16 27		27: 15 26 27 39		39: 27 38 40 51		51: 39 50 52 61
	// 4:  0  3  5  16		16: 4  15 17 28		28: 16 27 27 40		40: 28 39 41 52		52: 40 51 53 61
	// 5:  0  4  6  17		17: 5  16 18 29		29: 17 28 27 41		41: 29 40 42 53		53: 41 52 54 61
	// 6:  0  5  7  18		18: 6  17 19 30		30: 18 29 27 42		42: 30 41 43 54		54: 42 53 55 61
	// 7:  0  6  8  19		19: 7  18 20 31		31: 19 30 27 43		43: 31 41 44 55		55: 43 54 56 61
	// 8:  0  7  9  20		20: 8  19 21 32		32: 20 31 27 44		44: 32 43 45 56		56: 44 55 57 61
	// 9:  0  8  10 21		21: 9  20 22 33		33: 21 32 27 45		45: 33 44 46 57		57: 45 56 58 61
	// 10: 0  9  11 22		22: 10 21 23 34		34: 22 33 27 46		46: 34 45 47 58		58: 46 57 59 61
	// 11: 0  10 12 23		23: 11 22 24 35		35: 23 34 27 47		47: 35 46 48 59		59: 47 58 60 61
	// 12: 0  11 1  24		24: 12 23 13 36		36: 24 35 25 48		48: 36 47 37 60		60: 48 59 49 61
	
	// 61: 49 50 51 52 53 54 55 56 57 58 59 60
	
	// Setting the nodes to -1 so you can tell the nodes that where not used.
	// The first and the last nodes had 12 links so I made them all have 12.
	// The rest only had 4 so you may want to revisit this.
	for(int i = 0; i < NUMBER_OF_NODES; i++)
	{
		for(int j = 0; j < 12; j++)
		{
			NodeLinks[i][j] =  -1;
			NodeMuscles[i][j] = -1;
		}	
	}
	
	// Setting all nodes as not ablated
	for(int i = 0; i < NUMBER_OF_NODES; i++)
	{
		NodeAblatedYesNo[i] = 0;
	}
	
	// Setting edges for the 0th node.
	for(int i = 0; i < 12; i++)
	{
		NodeLinks[0][i] =  i + 1;
	}
	
	// Setting the edges that are connected to the 0th node
	for(int j = 0; j < 4; j++)
	{
		for(int i = 1; i < 13; i++)
		{
			// Connect to node above
			if(j == 0)
			{
				NodeLinks[i][j] =  0;
			}
			
			// Connect to the node to the left
			if(j == 1)
			{
				NodeLinks[i][j] =  (i+10)%12 + 1;
			}
			
			// Connect to the node to the right
			if(j == 2)
			{
				NodeLinks[i][j] =  (i+12)%12 + 1;
			}
			
			// Connect to the node below
			if(j == 3)
			{
				NodeLinks[i][j] =  i + 12;
			}
		}
		
		// Setting the middle 3 sections
		for(int k = 0; k <= 3*12; k += 12)
		{
			for(int i = 13 + k; i < 25 + k; i++)
			{
				// Connect to node above
				if(j == 0)
				{
					NodeLinks[i][j] =  i - 12;
				}
				
				// Connect to the node to the left
				if(j == 1)
				{
					NodeLinks[i][j] =  (i+10)%12 + 13 + k;
				}
				
				// Connect to the node to the right
				if(j == 2)
				{
					NodeLinks[i][j] =  (i+12)%12 + 13 + k;
				}
				
				// Connect to the node below
				if(j == 3)
				{
					NodeLinks[i][j] =  i + 12;
				}
			}
		}
		
		// Setting the edges that are linked to the 61st node
		for(int i = 49; i < 61; i++)
		{
			// Connect to node above
			if(j == 0)
			{
				NodeLinks[i][j] =  i - 12;
			}
			
			// Connect to the node to the left
			if(j == 1)
			{
				NodeLinks[i][j] =  (i+10)%12 + 49;
			}
			
			// Connect to the node to the right
			if(j == 2)
			{
				NodeLinks[i][j] =  (i+12)%12 + 49;
			}
			
			// Connect to the node below
			if(j == 3)
			{
				NodeLinks[i][j] =  NUMBER_OF_NODES - 1;
			}
		}
		
		// Setting the 61st node.
		for(int i = 0; i < 12; i++)
		{
			NodeLinks[61][i] =  i + 49;
		}
	}
	
	//Setting the ends of the muscles
	index = 0;
	for(int i = 0; i < NUMBER_OF_NODES; i++)
	{
		for(int j = 0; j < 12; j++)
		{
			if(NodeLinks[i][j] != -1)
			{
				if(i < NodeLinks[i][j])
				{
					if(NUMBER_OF_MUSCLES <= index)
					{
						printf("\nTSU Error: number of edges is out of bounds\n");
						return(0);
					} 
					MuscleConectionA[index] = i;
					MuscleConectionB[index] = NodeLinks[i][j];
					index++;
				}
			}
		}
	}
	
	// Setting the node muscles. Each node will have a list of nodes they are attached to (NodeLinks[][]) and the muscle that attaches it to that node (NodeMuscles[][]).
	for(int i = 0; i < NUMBER_OF_NODES; i++)
	{
		for(int j = 0; j < 12; j++)
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
		ContractionDuration[i] = 100.0;  // 100.0
		RelaxationDuration[i] = 200.0;  // 200.0
		ContractionStrength[i] = 0.1;
		MuscleColor[i].x = 1.0;
		MuscleColor[i].y = 0.0;
		MuscleColor[i].z = 0.0;
	}
	
	//Nodes to ablate
	for(int i = 35; i < 50; i++)
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
	
	// Setting the node masses
	for(int i = 0; i < NUMBER_OF_NODES; i++)
	{
		sum = 0.0;
		for(int j = 0; j < 12; j++)
		{
			if(NodeMuscles[i][j] != -1)
			{
				sum += MuscleMass[NodeMuscles[i][j]];
			}
		}
		NodeMass[i] = sum/2.0;
		printf("\nNodeMass[%d] = %f", i, NodeMass[i]);
	}

	return(1);
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
	
	// Drawing muscles
	glColor3d(1.0,0.0,0.0);
	for(int i = 0; i < NUMBER_OF_NODES; i++)
	{
		for(int j = 0; j < 12; j++)
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

void moveView()
{
	cout << "\nEnter the desired draw rate: ";
	cin >> DrawRate;
}

void generalMuscleForces()
{
	float f; 
	float dx, dy, dz, d;
	int muscleNumber, nodeNumber;
	
	// Getting forces on the nodes from the muscle fiber without contraction	
	for(int i = 0; i < NUMBER_OF_NODES; i++)
	{
		for(int j = 0; j < 12; j++)
		{
			muscleNumber = NodeMuscles[i][j];
			nodeNumber = NodeLinks[i][j];
			if(nodeNumber != -1)
			{
				dx = NodePosition[nodeNumber].x - NodePosition[i].x;
				dy = NodePosition[nodeNumber].y - NodePosition[i].y;
				dz = NodePosition[nodeNumber].z - NodePosition[i].z;
				d  = sqrt(dx*dx + dy*dy + dz*dz);
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
		
		f  = -BloodPresure;
		
		NodeForce[i].x  += f*dx/d;
		NodeForce[i].y  += f*dy/d;
		NodeForce[i].z  += f*dz/d;
	}
}

void turnOnNodeMuscles(int index)
{
	for(int j = 0; j < 12; j++)
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
		for(int j = 0; j < 12; j++)
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

int n_body()
{
	int   tdraw = 0; 
	double time = 0.0;
	float beatTimer = 0.0;
	
	//mousemask(ALL_MOUSE_EVENTS | REPORT_MOUSE_POSITION, NULL);
	//MEVENT event;
	
	//int mouse;
	
	while(time < STOP_TIME)
	{
		if(BeatPeriod <= beatTimer)
		{
			turnOnNodeMuscles(0);
			beatTimer = 0.0;
		}
		else beatTimer += DT;
		
		// Zeroing out the nodal forces.
		for(int i = 0; i < NUMBER_OF_NODES; i++)
		{
			NodeForce[i].x   = 0.0;
			NodeForce[i].y   = 0.0;
			NodeForce[i].z   = 0.0;
		}
		
		generalMuscleForces();
		
		contractionForces(DT, time);
		
		outwardPresure();
		
		dampingForce();
		
		moveNodes(DT, time);

		if(tdraw == DrawRate) 
		{
			draw_picture();
			tdraw = 0;
			printf("\n Time = %f", time);
		}
		else tdraw++;
		
		time += DT;
	}
	return(1);
}

void control()
{	
	//int    tdraw = 0;
	//float  time = 0.0;
	
	struct sigaction sa;
	
	sa.sa_handler = signalHandler;
	sigemptyset(&sa.sa_mask);
	sa.sa_flags = SA_RESTART; // Restart functions if interrupted by handler
	if (sigaction(SIGINT, &sa, NULL) == -1)
	{
		printf("\nTSU Error: sigaction error\n");
	}

	glClear(GL_COLOR_BUFFER_BIT);
	glClear(GL_DEPTH_BUFFER_BIT);

	set_initial_conditions();
	
	draw_picture();
	
	if(n_body() == 1) printf("\n N-body success \n");
	
	printf("\n DONE \n");
	while(1);
}

void mymouse(int button, int state, int x, int y)
{	
	if(state == GLUT_DOWN)
	{
		if(button == GLUT_LEFT_BUTTON)
		{
			printf("\n  x = %d", x);
		}
		else
		{
			printf("\n  y = %d", y);
		}
	}
}

void Display(void)
{
	gluLookAt(EyeX, EyeY, EyeZ, CenterX, CenterY, CenterZ, UpX, UpY, UpZ);
	glClear(GL_COLOR_BUFFER_BIT);
	glClear(GL_DEPTH_BUFFER_BIT);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glutSwapBuffers();
	glFlush();
	glutMouseFunc(mymouse);
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

static void signalHandler(int signum)
{
	int command;
    
	cout << "\n\n******************************************************" << endl;
	cout << "Enter:666 to kill the run." << endl;
	cout << "Enter:1 to change the draw rate." << endl;
	cout << "Enter:2 Move view" << endl;
	cout << "Enter:3 to continue the run." << endl;
	cout << "******************************************************\n\nCommand: ";
    
	cin >> command;
    
	if(command == 666)
	{
		exit(0);
	}
	else if(command == 1)
	{
		cout << "\nEnter the desired draw rate: ";
		cin >> DrawRate;
		cout << "\nDrawRate: " << DrawRate << endl;
	}
	else if (command == 2)
	{
		moveView();
	}
	else if (command == 3)
	{
		cout << "\nRun continued." << endl;
	}
	else
	{
		cout <<"\n\n Invalid Command\n" << endl;
	}
}

int main(int argc, char** argv)
{
	glutInit(&argc,argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGB);
	glutInitWindowSize(XWindowSize,YWindowSize);
	glutInitWindowPosition(0,0);
	glutCreateWindow("SVT Plus");
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

