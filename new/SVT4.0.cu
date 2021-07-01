// nvcc SVT4.0.cu -o svt4.0 -lglut -lm -lGLU -lGL
//To stop hit "control c" in the window you launched it from. stuff

// Length will be in millimeters
// Time will be in milliseconds
// Mass will be in ???

// Fiber length 100 micrometers or 0.1 millimeters
// Sodium wave speed .5 meters/sec or 0.5 millimeters/millisec

#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <GL/glut.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
using namespace std;

#define PI 3.141592654

//#define NumberOfNodes 266 //266 //62
//#define NumberOfMuscles 552 //552 //132
//#define LinksPerNode 24 //24

// Globals
float Dt;
int DrawRate;
int Pause;

int TypeOfShape;
int Divisions;

float Viscosity;
float BloodPresure;

float BeatPeriod;

float MassOfAtria;
float RadiusOfAtria = 1.0; // Should be 31.83098862

int NumberOfNodes;
int NumberOfMuscles;
int LinksPerNode;

float4 *NodePosition, *NodeVelocity, *NodeForce;
float *NodeMass;
int *NodeLinks; // The nodes that this node is connected to
int *NodeMuscles; // The muscle that connects this node to ther other nodes
float *NodeArea; // The surface area a node covers
int *NodeAblatedYesNo;

// How the muscle will act without contraction.
int *MuscleConectionA, *MuscleConectionB;
float *MuscleMass, *MuscleLength;
float *MuscleRelaxedStrength, BaseMuscleRelaxedStrength;
float *MuscleCompresionStopFraction, BaseMuscleCompresionStopFraction;
float3 *MuscleColor;
float MuscleCompresionMultiplier;
float MuscleTentionMultiplier;

// Muscle contraction parameters
int *ContractionOnOff;
float *ContractionTimer;
float *ActionPotentialSpeed, BaseActionPotentialSpeed; //0.5 is a good value.
float *ActionPotentialDuration;
float *ContractionDuration, BaseContractionDuration; // 100.0 is a good value
float *RelaxationDuration, BaseRelaxationDuration; // 200.0 is a good value
float *ContractionStrength, BaseContractionStrength; // 5.0 is a good value	

int   DrawTimer; 
float RunTime;
float BeatTimer;
float4 CenterOfSimulation;

// Window globals
static int Window;
int XWindowSize;
int YWindowSize; 
double Near;
double Far;
double EyeX;
double EyeY;
double EyeZ;
double CenterX;
double CenterY;
double CenterZ;
double UpX;
double UpY;
double UpZ;

// Prototyping functions
void allocateMemory(int, int);
void setNodesAndMusclesCircle(int); 
void setNodesAndMusclesSphere(int);
void linkNodesToMuscles();
void setMuscleAttributesAndNodeMasses(int, int);
void drawPicture();
void generalMuscleForces();
void outwardPresure();
void turnOnNodeMuscles(int);
int contractionForces(float, float);
void dampingForce();
void moveNodes(float, float);
void hardCodedAblatedNodes();
void hardCodedEctopicEvents(float, float);
void n_body(float);
void setup();
void KeyPressed(unsigned char, int, int);
void mymouse(int, int, int, int);
void Display(void);
void reshape(int, int);
void allocateMemory(int, int);
void simulationScript();
void readSimulationParameters();

#include "./setNodesAndMuscles.h"
#include "./callBackFunctions.h"

void readSimulationParameters()
{
/*
	MassOfAtria = 1.0;// Need to look this up. ????????????
	MuscleCompresionMultiplier = 50.0;  // How hard a muscle resists being compressed past its compression limit.
	MuscleTentionMultiplier = 50.0;  // How hard a musle pulls back when it is stretched past its natural length.
	Viscosity = 10.0; // Jsut something to give resistance to movement.This will be divided by the number of nodes for scalling. ????????
	BloodPresure = 1.0; This will be scaled by the number of noddes. ?????????????????
	
	BaseMuscleRelaxedStrength = 0.1; // This will be scaled by multiplying by muscle length. This is the standard but can be adjusted for each muscle.
	BaseMuscleCompresionStopFraction = 0.7; // The percentage a muscles length can contract. This is the standard but can be adjusted for each muscle.
	BaseActionPotentialSpeed = 0.2; // This is the speed of the action potential. This is the standard but can be adjusted for each muscle.	
		
	BaseContractionDuration = 20.0;  // 100.0
	BaseRelaxationDuration = 40.0;  // 200.0
	BaseContractionStrength = 0.2; // This will be scaled by multipling by muscle length.
*/
	
	ifstream data;
	string name;
	
	data.open("./simulationSetup");
	
	if(data.is_open() == 1)
	{
		getline(data,name,'=');
		data >> TypeOfShape;
		
		getline(data,name,'=');
		data >> Divisions;
		
		getline(data,name,'=');
		data >> Viscosity;
		
		getline(data,name,'=');
		data >> BloodPresure;
		
		getline(data,name,'=');
		data >> MassOfAtria;
		
		getline(data,name,'=');
		data >> MuscleCompresionMultiplier;
		
		getline(data,name,'=');
		data >> MuscleTentionMultiplier;
		
		getline(data,name,'=');
		data >> BaseMuscleRelaxedStrength;
		
		getline(data,name,'=');
		data >> BaseContractionStrength;
		
		getline(data,name,'=');
		data >> BaseMuscleCompresionStopFraction;
		
		getline(data,name,'=');
		data >> BaseContractionDuration;
		
		getline(data,name,'=');
		data >> BaseRelaxationDuration;
		
		getline(data,name,'=');
		data >> BaseActionPotentialSpeed;
		
		getline(data,name,'=');
		data >> BeatPeriod;
		
		getline(data,name,'=');
		data >> DrawRate;
		
		getline(data,name,'=');
		data >> Dt;
	}
	else
	{
		printf("\nTSU Error could not open simulationSetup file\n");
		exit(0);
	}
	data.close();
	
	if(TypeOfShape == 1)
	{	
		if(Divisions == 0)
		{
			printf("\n So you want to run a simulation with nothing in it.");
			printf("\n That's easy just look at a blank screen. \n");
			printf("\n Good Bye. \n");
			exit(0);
		}
		if(Divisions == 1)
		{
			printf("\n Seriously a circle of 1!");
			printf("\n This is sad. You need to get out make some friends. \n");
			printf("\n Good Bye. \n");
			exit(0);
		}
		printf("\n You will be simulating a circle with %d divisions\n", Divisions);
	}
	else if(TypeOfShape == 2)
	{
		if(Divisions%2 != 0)
		{
			printf("\n I said the number had to be even!");
			printf("\n Beem me up Scotty. There is no intelligent life down here. \n");
			printf("\n Good Bye. \n");
			exit(0);
		}
		else if(Divisions < 5)
		{
			printf("\n Yo Einstien! I said the number had to be even and greater than or equal to 4.\n");
			printf("\n Good Bye. \n");
			exit(0);
		}
	}
	else
	{
		printf("\n Type of simulation is incorrect. \n");
		printf("\n Good Bye. \n");
		exit(0);
	}
}

void allocateMemory(int type, int divisions)
{
	if(type == 1) // Circle
	{
		NumberOfNodes = divisions;
		NumberOfMuscles = divisions;
		LinksPerNode = 2;
	}
	else if(type == 2) //Sphere
	{
		NumberOfNodes = divisions*(divisions/2 - 1) + 2;
		NumberOfMuscles = divisions + (divisions*2)*(divisions/2 - 1);
		LinksPerNode = divisions;
	}
	else if(type == 3) //Sphere with thickness
	{
		printf("\n Thick spheres are not in yet.\n");
		exit(0);
	}
	else
	{
		printf("\n Bad object type.\n");
		exit(0);
	}
	
	NodePosition = (float4*)malloc(NumberOfNodes*sizeof(float4));
	NodeVelocity = (float4*)malloc(NumberOfNodes*sizeof(float4));
	NodeForce    = (float4*)malloc(NumberOfNodes*sizeof(float4));
	
	NodeMass = (float*)malloc(NumberOfNodes*sizeof(float));
	NodeArea = (float*)malloc(NumberOfNodes*sizeof(float));
	NodeLinks = (int*)malloc(NumberOfNodes*LinksPerNode*sizeof(int));
	NodeMuscles = (int*)malloc(NumberOfNodes*LinksPerNode*sizeof(int));
	NodeAblatedYesNo = (int*)malloc(NumberOfNodes*sizeof(int));

	// How the muscle will act without contraction.
	MuscleConectionA = (int*)malloc(NumberOfMuscles*sizeof(int));
	MuscleConectionB = (int*)malloc(NumberOfMuscles*sizeof(int));
	MuscleMass = (float*)malloc(NumberOfMuscles*sizeof(float));
	MuscleLength = (float*)malloc(NumberOfMuscles*sizeof(float));
	MuscleRelaxedStrength = (float*)malloc(NumberOfMuscles*sizeof(float));
	MuscleCompresionStopFraction = (float*)malloc(NumberOfMuscles*sizeof(float));
	MuscleColor = (float3*)malloc(NumberOfMuscles*sizeof(float3));


	// Muscle contraction parameters
	ContractionOnOff = (int*)malloc(NumberOfMuscles*sizeof(int));
	ContractionTimer = (float*)malloc(NumberOfMuscles*sizeof(float));
	ActionPotentialSpeed = (float*)malloc(NumberOfMuscles*sizeof(float));
	ActionPotentialDuration = (float*)malloc(NumberOfMuscles*sizeof(float));
	ContractionDuration = (float*)malloc(NumberOfMuscles*sizeof(float));
	RelaxationDuration = (float*)malloc(NumberOfMuscles*sizeof(float));
	ContractionStrength = (float*)malloc(NumberOfMuscles*sizeof(float));
}

void linkNodesToMuscles()
{
	//Setting the ends of the muscles
	int index = 0;
	for(int i = 0; i < NumberOfNodes; i++)
	{
		for(int j = 0; j < LinksPerNode; j++)
		{
			if(NodeLinks[i*LinksPerNode + j] != -1)
			{
				if(i < NodeLinks[i*LinksPerNode + j])
				{
					if(NumberOfMuscles <= index)
					{
						printf("\nTSU Error: number of muscles is out of bounds\n");
						exit(0);
					} 
					MuscleConectionA[index] = i;
					MuscleConectionB[index] = NodeLinks[i*LinksPerNode + j];
					index++;
				}
			}
		}
	}
	
	// Setting the node muscles. Each node will have a list of nodes they are attached to (NodeLinks[][]) and the muscle that attaches it to that node (NodeMuscles[][]).
	// Setting them all to -1 first.
	for(int i = 0; i < NumberOfNodes; i++)
	{
		for(int j = 0; j < LinksPerNode; j++)
		{
			NodeMuscles[i*LinksPerNode + j] = -1;
		}	
	}
	for(int i = 0; i < NumberOfNodes; i++)
	{
		for(int j = 0; j < LinksPerNode; j++)
		{
			if(NodeLinks[i*LinksPerNode + j] != -1)
			{
				for(int k = 0; k < NumberOfMuscles; k++)
				{
					if((MuscleConectionA[k] == i && MuscleConectionB[k] == NodeLinks[i*LinksPerNode + j]) || (MuscleConectionA[k] == NodeLinks[i*LinksPerNode + j] && MuscleConectionB[k] == i))
					{
						if(NumberOfNodes*LinksPerNode <= (i*LinksPerNode + j))
						{
							printf("\nTSU Error: number of muscles is out of bounds\n");
							exit(0);
						} 
						NodeMuscles[i*LinksPerNode + j] = k;
					}
				}
			}
		}
	}
}

void setMuscleAttributesAndNodeMasses(int type, int divisions)
{	
	float dx, dy, dz, d, d1, d2;
	float sum, totalLengthOfAllMuscles;
	float bloodPresureScaling;
	float surfaceArea;
	
	Viscosity /= NumberOfNodes;
	
	/*
	for(int i = 0; i < NumberOfNodes; i++)
	{
		dx = NodePosition[NodeLinks[i*LinksPerNode + 0]].x - NodePosition[i*LinksPerNode + 3].x;
		dy = NodePosition[NodeLinks[i*LinksPerNode + 0]].y - NodePosition[i*LinksPerNode + 3].y;
		dz = NodePosition[NodeLinks[i*LinksPerNode + 0]].z - NodePosition[i*LinksPerNode + 3].z;
		d1 = sqrt(dx*dx + dy*dy + dz*dz)/2.0;
		dx = NodePosition[NodeLinks[i*LinksPerNode + 1]].x - NodePosition[i*LinksPerNode + 2].x;
		dy = NodePosition[NodeLinks[i*LinksPerNode + 1]].y - NodePosition[i*LinksPerNode + 2].y;
		dz = NodePosition[NodeLinks[i*LinksPerNode + 1]].z - NodePosition[i*LinksPerNode + 2].z;
		d2 = sqrt(dx*dx + dy*dy + dz*dz)/2.0;
		NodeArea[i] = d1*d2;
		printf("\n node area[%d] = %f", i, NodeArea[i]);
	}
	*/
	
	surfaceArea = 4.0*PI*RadiusOfAtria*RadiusOfAtria;
	bloodPresureScaling = surfaceArea/NumberOfNodes; // Need to scale by density too ??????????????
	BloodPresure *= bloodPresureScaling;
	
	CenterOfSimulation.x = 0.0;
	CenterOfSimulation.y = 0.0;
	CenterOfSimulation.z = 0.0;
	
	//Finding the length of each muscle and the total length of all muscles.
	totalLengthOfAllMuscles = 0.0;
	for(int i = 0; i < NumberOfMuscles; i++)
	{	
		dx = NodePosition[MuscleConectionA[i]].x - NodePosition[MuscleConectionB[i]].x;
		dy = NodePosition[MuscleConectionA[i]].y - NodePosition[MuscleConectionB[i]].y;
		dz = NodePosition[MuscleConectionA[i]].z - NodePosition[MuscleConectionB[i]].z;
		d = sqrt(dx*dx + dy*dy + dz*dz);
		MuscleLength[i] = d;
		totalLengthOfAllMuscles += d;
	}
	
	// Setting the mass of all muscles.
	if(type == 1)
	{
		MassOfAtria /= divisions; // If you are on a circle. There are division circle that make up the sphere so the circle is 1/divsiions the total mass.
	}
	for(int i = 0; i < NumberOfMuscles; i++)
	{	
		MuscleMass[i] = (MuscleLength[i]/totalLengthOfAllMuscles)*MassOfAtria;
	}
	
	// Setting other parameters
	for(int i = 0; i < NumberOfMuscles; i++)
	{	
		MuscleRelaxedStrength[i] = BaseMuscleRelaxedStrength*MuscleLength[i];
		MuscleCompresionStopFraction[i] = BaseMuscleCompresionStopFraction;
		ContractionOnOff[i] = 0;
		ContractionTimer[i] = 0.0;
		ActionPotentialSpeed[i] = BaseActionPotentialSpeed; //0.2; // 0.2
		ActionPotentialDuration[i] = MuscleLength[i]/ActionPotentialSpeed[i];
		ContractionDuration[i] = BaseContractionDuration; //20.0;  // 100.0
		RelaxationDuration[i] = BaseRelaxationDuration; //40.0;  // 200.0
		ContractionStrength[i] = BaseContractionStrength*MuscleLength[i]; //0.1;
		
		MuscleColor[i].x = 1.0;
		MuscleColor[i].y = 0.0;
		MuscleColor[i].z = 0.0;
	}
	
	for(int i = 0; i < NumberOfNodes; i++)
	{
		NodeAblatedYesNo[i] = 0; // Setting all nodes to not ablated.
	}
	
	// Setting the node masses
	for(int i = 0; i < NumberOfNodes; i++)
	{
		sum = 0.0;
		for(int j = 0; j < LinksPerNode; j++)
		{
			if(NodeMuscles[i*LinksPerNode + j] != -1)
			{
				sum += MuscleMass[NodeMuscles[i*LinksPerNode + j]];
			}
		}
		NodeMass[i] = sum/2.0;
	}
}

void hardCodedAblatedNodes()
{
	//Nodes to ablate
	for(int i = 1; i < 23; i++)
	{	
		//NodeAblatedYesNo[i] = 1;
	}
	
	//NodeAblatedYesNo[49] = 1;
}

void hardCodedEctopicEvents(float time, float dt)
{
	float er = dt/2.0;
	
	if((50.0 - er <= time) && (time < 220.0 + er))
	{
		//turnOnNodeMuscles(31);
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

void drawPicture()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glClear(GL_DEPTH_BUFFER_BIT);
	
	glColor3d(1.0,1.0,0.0);
	glPushMatrix();
	glTranslatef(NodePosition[0].x, NodePosition[0].y, NodePosition[0].z);
	glutSolidSphere(0.03,20,20);
	glPopMatrix();
	
	// Drawing nodes
	for(int i = 1; i < NumberOfNodes; i++)
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
	glTranslatef(NodePosition[NumberOfNodes-1].x, NodePosition[NumberOfNodes-1].y, NodePosition[NumberOfNodes-1].z);
	glutSolidSphere(0.03,20,20);
	glPopMatrix();
	*/
	
	// Drawing muscles
	glColor3d(1.0,0.0,0.0);
	for(int i = 0; i < NumberOfNodes; i++)
	{
		for(int j = 0; j < LinksPerNode; j++)
		{
			if(NodeLinks[i*LinksPerNode + j] != -1)
			{
				//glColor3d(MuscleColor[NodeMuscles[i*LinksPerNode + j]].x, MuscleColor[NodeMuscles[i*LinksPerNode + j]].y, MuscleColor[NodeMuscles[i*LinksPerNode + j]].z);
				//glLineWidth(1.0/(Px[i]-Px[NodeLinks[i*LinksPerNode + j]]));
				glBegin(GL_LINES);
					glVertex3f(NodePosition[i].x, NodePosition[i].y, NodePosition[i].z);
					glVertex3f(NodePosition[NodeLinks[i*LinksPerNode + j]].x, NodePosition[NodeLinks[i*LinksPerNode + j]].y, NodePosition[NodeLinks[i*LinksPerNode + j]].z);
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
	for(int i = 0; i < NumberOfNodes; i++)
	{
		for(int j = 0; j < LinksPerNode; j++)
		{
			muscleNumber = NodeMuscles[i*LinksPerNode + j];
			nodeNumber = NodeLinks[i*LinksPerNode + j];
			if(nodeNumber != -1)
			{
				dx = NodePosition[nodeNumber].x - NodePosition[i].x;
				dy = NodePosition[nodeNumber].y - NodePosition[i].y;
				dz = NodePosition[nodeNumber].z - NodePosition[i].z;
				d  = sqrt(dx*dx + dy*dy + dz*dz);
				
				// Grabbing numeric overflow before it happens.
				if(d < 0.00001) 
				{
					printf("\n TSU Error: In generalMuscleForces d is very small between nodeNumbers %d and %d the seperation is %f. Take a look at this!\n", i, nodeNumber, d);
					glColor3d(0.0,0.0,1.0);
					
					// Displaying where the problem occured.
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
	double centerOfMassX, centerOfMassY, centerOfMassZ, mass;
	
	centerOfMassX = 0.0;
	centerOfMassY = 0.0;
	centerOfMassZ = 0.0;
	mass = 0.0;
	for(int i = 0; i < NumberOfNodes; i++)
	{
		 centerOfMassX += NodePosition[i].x*NodeMass[i];
		 centerOfMassY += NodePosition[i].y*NodeMass[i];
		 centerOfMassZ += NodePosition[i].z*NodeMass[i];
		 mass += NodeMass[i];
	}
	centerOfMassX /= mass;
	centerOfMassY /= mass;
	centerOfMassZ /= mass;
	
	// Getting forces on the nodes from the presure of the blood pushing out	
	for(int i = 0; i < NumberOfNodes; i++)
	{
		dx = centerOfMassX - NodePosition[i].x;
		dy = centerOfMassY - NodePosition[i].y;
		dz = centerOfMassZ - NodePosition[i].z;
		d  = sqrt(dx*dx + dy*dy + dz*dz);
		
		// Grabbing numeric overflow before it happens.
		if(d < 0.0001) 
		{
			printf("\nTSU Error: In outwardPresure d is very small. Take a look at this\n");
			exit(0);
		}
		
		f  = -BloodPresure;   //*NodeArea[i];
		
		NodeForce[i].x  += f*dx/d;
		NodeForce[i].y  += f*dy/d;
		NodeForce[i].z  += f*dz/d;
	}
}

void turnOnNodeMuscles(int index)
{
	for(int j = 0; j < LinksPerNode; j++)
	{
		if((NodeLinks[index*LinksPerNode + j] != -1) && (ContractionOnOff[NodeMuscles[index*LinksPerNode + j]] == 0))
		{
			ContractionOnOff[NodeMuscles[index*LinksPerNode + j]] = 1;
			ContractionTimer[NodeMuscles[index*LinksPerNode + j]] = 0.0;
		}
	}
}

int contractionForces(float dt, float time)
{
	float dx, dy, dz, d;
	int muscleNumber, nodeNumber;
	
	// Getting forces for the muscle fiber contraction
	for(int i = 0; i < NumberOfNodes; i++)
	{
		for(int j = 0; j < LinksPerNode; j++)
		{
			muscleNumber = NodeMuscles[i*LinksPerNode + j];
			nodeNumber = NodeLinks[i*LinksPerNode + j];
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
						
						// Grabbing numeric overflow before it happens.
						if(d < 0.00001) 
						{
							printf("\n TSU Error: In contractionForces d is very small between nodeNumbers = %d %d seperation = %f. Take a look at this\n", i, nodeNumber, d);
							glColor3d(0.0,0.0,1.0);
							
							// Displaying where the problem occured.
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
	for(int i = 0; i < NumberOfNodes; i++)
	{	
		NodeForce[i].x   += -Viscosity*NodeVelocity[i].x;
		NodeForce[i].y   += -Viscosity*NodeVelocity[i].y;
		NodeForce[i].z   += -Viscosity*NodeVelocity[i].z;
	}
}

void moveNodes(float dt, float time)  // LeapFrog
{
	// Moving the system forward in time with leap-frog.
	for(int i = 0; i < NumberOfNodes; i++)
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

void n_body(float dt)
{	
	if(Pause != 1)
	{
		if(BeatPeriod <= BeatTimer)
		{
			turnOnNodeMuscles(0);
			BeatTimer = 0.0;
		}
		else BeatTimer += dt;
		
		hardCodedEctopicEvents(RunTime, dt);
		
		// Zeroing out the nodal forces.
		for(int i = 0; i < NumberOfNodes; i++)
		{
			NodeForce[i].x   = 0.0;
			NodeForce[i].y   = 0.0;
			NodeForce[i].z   = 0.0;
		}
		
		generalMuscleForces();
		
		contractionForces(dt, RunTime);
		
		outwardPresure();
		
		dampingForce();
		
		moveNodes(dt, RunTime);

		if(DrawTimer == DrawRate) 
		{
			drawPicture();
			DrawTimer = 0;
			printf("\n Time = %f", RunTime);
		}
		else DrawTimer++;
		
		RunTime += dt;
	}
}

void simulationScript()
{
	printf("\n\n\n The Particle Modeling Group hopes you injoy your interactive right atriam simulation.\n\n");
	printf("\n The simulation is paused.");
	printf("\n Move to the mouse over the simulation window and type the following commands.\n");
	printf("\n To run the simulation type r.");
	printf("\n To pause the simulation type p.");
	printf("\n The positive x-axis is to the right.");
	printf("\n The positive y-axis is up.");
	printf("\n The positive z-axis is towards you.");
	printf("\n For an orthoganal view type 0.");
	printf("\n For a fulstrum view type f");
	printf("\n To do a positive spin on the x-axis type x (negative X).");
	printf("\n To do a positive spin on the y-axis type y (negative Y).");
	printf("\n To do a positive spin on the Z-axis type z (negative Z).");
	printf("\n To zoom in type w (out W). Note zoom is meaningless in orthoganal mode.");
	printf("\n To center type c");
	printf("\n To center and out the sinus node up type n");
	printf("\n To ablate or unablate right click the mouse on the node you are interested in");
	printf("\n For best ablation results, pause the simulation and put it in orthaganal mode.");
	printf("\n To quit the simulation type q or hit the kill button on the window.");
	printf("\n\n Happy ablating!\n");
}

void setup()
{	
	readSimulationParameters();
	
	allocateMemory(TypeOfShape, Divisions);
	
	if(TypeOfShape == 1) setNodesAndMusclesCircle(Divisions);
	else if(TypeOfShape == 2) setNodesAndMusclesSphere(Divisions);
	
	linkNodesToMuscles();
	
	setMuscleAttributesAndNodeMasses(TypeOfShape, Divisions);
	
	hardCodedAblatedNodes();
	
	DrawTimer = 0; 
	RunTime = 0.0;
	BeatTimer = 0.0;
	Pause = 1;
	
	simulationScript();
}

int main(int argc, char** argv)
{
	XWindowSize = 1000;
	YWindowSize = 1000; 

	// Clip plains
	Near = 0.2;
	Far = 80.0;

	//Direction here your eye is located location
	EyeX = 0.0;
	EyeY = 0.0;
	EyeZ = 2.0;

	//Where you are looking
	CenterX = 0.0;
	CenterY = 0.0;
	CenterZ = 0.0;

	//Up vector for viewing
	UpX = 0.0;
	UpY = 1.0;
	UpZ = 0.0;
	
	//setup();
	glutInit(&argc,argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGB);
	glutInitWindowSize(XWindowSize,YWindowSize);
	glutInitWindowPosition(5,5);
	Window = glutCreateWindow("SVT");
	gluLookAt(EyeX, EyeY, EyeZ, CenterX, CenterY, CenterZ, UpX, UpY, UpZ);
	
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glFrustum(-0.2, 0.2, -0.2, 0.2, Near, Far);
	glMatrixMode(GL_MODELVIEW);
	glClearColor(0.0, 0.0, 0.0, 0.0);
	
	GLfloat light_position[] = {1.0, 1.0, 1.0, 0.0};
	GLfloat light_ambient[]  = {0.0, 0.0, 0.0, 1.0};
	GLfloat light_diffuse[]  = {1.0, 1.0, 1.0, 1.0};
	GLfloat light_specular[] = {1.0, 1.0, 1.0, 1.0};
	GLfloat lmodel_ambient[] = {0.2, 0.2, 0.2, 1.0};
	GLfloat mat_specular[]   = {1.0, 1.0, 1.0, 1.0};
	GLfloat mat_shininess[]  = {10.0};
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
	glutMouseFunc(mymouse);
	glutKeyboardFunc(KeyPressed);
	glutIdleFunc(idle);
	setup();
	glutMainLoop();
	return 0;
}

