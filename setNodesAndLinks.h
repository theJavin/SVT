
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

void initializeNodesAndLinksSphere(int divitions) 
{	
	int index;
	
	int NumberOfNodes = divitions*(divitions/2 - 1) + 2;
	int NumberOfEdges = divitions + (divitions*2)*(divitions/2 - 1);
	int LinksPerNode = divitions;
	
	printf("\n number of nodes = %d", NumberOfNodes);
	printf("\n number of edges = %d", NumberOfEdges);
	printf("\n number of LinksPerNode = %d", LinksPerNode);
	printf("\n");
	//float NumberOfLinks = ???;
	
	// Node position values for a sphere with 62 nodes//0.5 is a good value.
	NodePosition[0].x = 0.0;
	NodePosition[0].y = 1.0;
	NodePosition[0].z = 0.0;
	NodePosition[NumberOfNodes-1].x = 0.0;
	NodePosition[NumberOfNodes-1].y = -1.0;
	NodePosition[NumberOfNodes-1].z = 0.0;
	
	index = 1;
	for(int i = 1; i < divitions/2; i++)
	{
		for(int j = 0; j < divitions; j++)
		{
			if((NumberOfNodes-1) <= index)
			{
				printf("\nTSU Error: number of nodes is out of bounds\n");
				exit(0);
			} 
			NodePosition[index].y = sin(PI/2.0 -i*PI/((float)divitions/2.0));
			NodePosition[index].x = cos(PI/2.0 -i*PI/((float)divitions/2.0))*cos(j*PI/((float)divitions/2.0));
			NodePosition[index].z = cos(PI/2.0 -i*PI/((float)divitions/2.0))*sin(j*PI/((float)divitions/2.0));
			
			index++;
		}	
	}
	
	// Zeroing out velocity and acceleration
	for(int i = 0; i < NumberOfNodes; i++)
	{
		NodeVelocity[index].y = 0.0;
		NodeVelocity[index].x = 0.0;
		NodeVelocity[index].z = 0.0;
		
		NodeForce[index].y = 0.0;
		NodeForce[index].x = 0.0;
		NodeForce[index].z = 0.0;
	}
	
	// Setting the nodes to -1 so you can tell the nodes that where not used.
	// The first and the last nodes had 12 links so I made them all have 12.
	// The rest only had 4 so you may want to revisit this.
	for(int i = 0; i < NumberOfNodes; i++)
	{
		for(int j = 0; j < LinksPerNode; j++)
		{
			NodeLinks[i][j] =  -1;
		}	
	}
	
	// Setting edges for the 0th node.
	for(int j = 0; j < divitions; j++)
	{
		NodeLinks[0][j] =  j + 1;
	}
	
	// Setting the edges that are connected to the 0th node
	for(int j = 0; j < 4; j++)
	{
		for(int i = 1; i < divitions + 1; i++)
		{
			// Connect to node above
			if(j == 0)
			{
				NodeLinks[i][j] =  0;
			}
			
			// Connect to the node to the left
			if(j == 1)
			{NUMBER_OF_MUSCLES
				NodeLinks[i][j] = (divitions + i - 2)%divitions + 1;
			}
			
			// Connect to the node to the right
			if(j == 2)
			{
				NodeLinks[i][j] =  (divitions + i)%divitions + 1;
			}
			
			// Connect to the node below
			if(j == 3)
			{
				NodeLinks[i][j] =  i + divitions;
			}
		}
		
		// Setting the middle sections
		for(int k = 0; k < (divitions/2 - 3)*divitions; k += divitions)
		{
			for(int i = divitions + 1 + k; i <= 2*divitions + k; i++)
			{
				// Connect to node above
				if(j == 0)
				{
					NodeLinks[i][j] =  i - divitions;
				}
				
				// Connect to the node to the left
				if(j == 1)
				{
					NodeLinks[i][j] =  (i + divitions - 2)%divitions + divitions + 1 + k;
				}
				
				// Connect to the node to the right
				if(j == 2)
				{
					NodeLinks[i][j] =  (i + divitions)%divitions + divitions + 1 + k;
				}
				
				// Connect to the node below
				if(j == 3)
				{
					NodeLinks[i][j] =  i + divitions;
				}
			}
		}
		
		// Setting the edges that are linked to the last node
		for(int i = NumberOfNodes -1 - divitions; i < NumberOfNodes - 1; i++)
		{
			// Connect to node above
			if(j == 0)
			{
				NodeLinks[i][j] =  i - divitions;
			}
			
			// Connect to the node to the left
			if(j == 1)
			{
				NodeLinks[i][j] = (i + divitions - 2)%divitions + NumberOfNodes -1 - divitions;
			}
			
			// Connect to the node to the right
			if(j == 2)
			{
				NodeLinks[i][j] =  (i+divitions)%divitions +NumberOfNodes -1 - divitions;
			}
			
			// Connect to the node below
			if(j == 3)
			{
				NodeLinks[i][j] =  NumberOfNodes - 1;
			}
		}
	}
	
	// Setting the last node.
	for(int j = 0; j < divitions; j++)
	{
		NodeLinks[NumberOfNodes - 1][j] =  NumberOfNodes - 1 - divitions + j;
	}
	
	/*
	for(int i = 0; i < NumberOfNodes; i++)
	{
		for(int j = 0; j < LinksPerNode; j++)
		{
			printf("\n NodeLinks[%d][%d] = %d", i, j, NodeLinks[i][j]);
		}
	}
	*/
	printf("\n Nodes created");
	//while(1);
}