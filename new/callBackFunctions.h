void KeyPressed(unsigned char key, int x, int y)
{
	float dAngle = 0.01;
	float zoom = 0.01;
	float temp;
	float4 lookVector;
	float d;
	float4 centerOfMass;
	
	lookVector.x = CenterX - EyeX;
	lookVector.y = CenterY - EyeY;
	lookVector.z = CenterZ - EyeZ;
	d = sqrt(lookVector.x*lookVector.x + lookVector.y*lookVector.y + lookVector.z*lookVector.z);
	lookVector.x /= d;
	lookVector.y /= d;
	lookVector.z /= d;
	
	centerOfMass.x = 0.0;
	centerOfMass.y = 0.0;
	centerOfMass.z = 0.0;
	for(int i = 0; i < NumberOfNodes; i++)
	{
		 centerOfMass.x += NodePosition[i].x*NodeMass[i];
		 centerOfMass.y += NodePosition[i].y*NodeMass[i];
		 centerOfMass.z += NodePosition[i].z*NodeMass[i];
		 centerOfMass.w += NodeMass[i];
	}
	centerOfMass.x /= centerOfMass.w;
	centerOfMass.y /= centerOfMass.w;
	centerOfMass.z /= centerOfMass.w;
	
	
	if(key == 'q')
	{
		glutDestroyWindow(Window);
		printf("\nw Good Bye\n");
		exit(0);
	}
	if(key == 'x')
	{
		for(int i = 0; i < NumberOfNodes; i++)
		{
			NodePosition[i].y -= centerOfMass.y;
			NodePosition[i].z -= centerOfMass.z;
			temp 			   = cos(dAngle)*NodePosition[i].y - sin(dAngle)*NodePosition[i].z;
			NodePosition[i].z  = sin(dAngle)*NodePosition[i].y + cos(dAngle)*NodePosition[i].z;
			NodePosition[i].y  = temp;
			NodePosition[i].y += centerOfMass.y;
			NodePosition[i].z += centerOfMass.z;
		}
		drawPicture();
	}
	if(key == 'X')
	{
		for(int i = 0; i < NumberOfNodes; i++)
		{
			NodePosition[i].y -= centerOfMass.y;
			NodePosition[i].z -= centerOfMass.z;
			temp			   = cos(-dAngle)*NodePosition[i].y - sin(-dAngle)*NodePosition[i].z;
			NodePosition[i].z  = sin(-dAngle)*NodePosition[i].y + cos(-dAngle)*NodePosition[i].z;
			NodePosition[i].y  = temp; 
			NodePosition[i].y += centerOfMass.y;
			NodePosition[i].z += centerOfMass.z;
		}
		drawPicture();
	}
	if(key == 'y')
	{
		for(int i = 0; i < NumberOfNodes; i++)
		{
			NodePosition[i].x -= centerOfMass.x;
			NodePosition[i].z -= centerOfMass.z;
			temp			   = cos(dAngle)*NodePosition[i].x + sin(dAngle)*NodePosition[i].z;
			NodePosition[i].z  = -sin(dAngle)*NodePosition[i].x + cos(dAngle)*NodePosition[i].z;
			NodePosition[i].x  = temp;
			NodePosition[i].x += centerOfMass.x;
			NodePosition[i].z += centerOfMass.z;
		}
		drawPicture();
	}
	if(key == 'Y')
	{
		for(int i = 0; i < NumberOfNodes; i++)
		{
			NodePosition[i].x -= centerOfMass.x;
			NodePosition[i].z -= centerOfMass.z;
			temp			   = cos(-dAngle)*NodePosition[i].x + sin(-dAngle)*NodePosition[i].z;
			NodePosition[i].z  = -sin(-dAngle)*NodePosition[i].x + cos(-dAngle)*NodePosition[i].z;
			NodePosition[i].x  = temp;
			NodePosition[i].x += centerOfMass.x;
			NodePosition[i].z += centerOfMass.z;
		}
		drawPicture();
	}
	if(key == 'z')
	{
		for(int i = 0; i < NumberOfNodes; i++)
		{
			NodePosition[i].x -= centerOfMass.x;
			NodePosition[i].y -= centerOfMass.y;
			temp			   = cos(dAngle)*NodePosition[i].x - sin(dAngle)*NodePosition[i].y;
			NodePosition[i].y  = sin(dAngle)*NodePosition[i].x + cos(dAngle)*NodePosition[i].y;
			NodePosition[i].x  = temp;
			NodePosition[i].x += centerOfMass.x;
			NodePosition[i].y += centerOfMass.y;
		}
		drawPicture();
	}
	if(key == 'Z')
	{
		for(int i = 0; i < NumberOfNodes; i++)
		{
			NodePosition[i].x -= centerOfMass.x;
			NodePosition[i].y -= centerOfMass.y;
			temp			   = cos(-dAngle)*NodePosition[i].x - sin(-dAngle)*NodePosition[i].y;
			NodePosition[i].y  = sin(-dAngle)*NodePosition[i].x + cos(-dAngle)*NodePosition[i].y;
			NodePosition[i].x  = temp;
			NodePosition[i].x += centerOfMass.x;
			NodePosition[i].y += centerOfMass.y;
		}
		drawPicture();
	}
	if(key == 'w')
	{
		printf("\n look x y z %f  %f  %f", lookVector.x, lookVector.y, lookVector.z);
		for(int i = 0; i < NumberOfNodes; i++)
		{
			NodePosition[i].x -= zoom*lookVector.x;
			NodePosition[i].y -= zoom*lookVector.y;
			NodePosition[i].z -= zoom*lookVector.z;
			
			CenterOfSimulation.x -= zoom*lookVector.x;
			CenterOfSimulation.y -= zoom*lookVector.y;
			CenterOfSimulation.z -= zoom*lookVector.z;
		}
		drawPicture();
	}
	if(key == 's')
	{
		for(int i = 0; i < NumberOfNodes; i++)
		{
			NodePosition[i].x += zoom*lookVector.x;
			NodePosition[i].y += zoom*lookVector.y;
			NodePosition[i].z += zoom*lookVector.z;
			
			CenterOfSimulation.x -= zoom*lookVector.x;
			CenterOfSimulation.y -= zoom*lookVector.y;
			CenterOfSimulation.z -= zoom*lookVector.z;
		}
		drawPicture();
	}
	
	if(key == 'o')
	{
		centerOfMass.x = 0.0;
		centerOfMass.y = 0.0;
		centerOfMass.z = 0.0;
		for(int i = 0; i < NumberOfNodes; i++)
		{
			 centerOfMass.x += NodePosition[i].x*NodeMass[i];
			 centerOfMass.y += NodePosition[i].y*NodeMass[i];
			 centerOfMass.z += NodePosition[i].z*NodeMass[i];
			 centerOfMass.w += NodeMass[i];
		}
		centerOfMass.x /= centerOfMass.w;
		centerOfMass.y /= centerOfMass.w;
		centerOfMass.z /= centerOfMass.w;
		
		for(int i = 0; i < NumberOfNodes; i++)
		{
			NodePosition[i].x += centerOfMass.x;
			NodePosition[i].y -= centerOfMass.y;
			NodePosition[i].z -= centerOfMass.z;
		}
		
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glOrtho(-1.0, 1.0, -1.0, 1.0, 0.2, 80.0);
		glMatrixMode(GL_MODELVIEW);
		drawPicture();
	}
	if(key == 'f')
	{
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glFrustum(-0.2, 0.2, -0.2, 0.2, 0.2, 80.0);
		glMatrixMode(GL_MODELVIEW);
		drawPicture();
	}
	
	if(key == 'p')
	{
		Pause = 1;
	}
	if(key == 'r')
	{
		Pause = 0;
	}
	
	if(key == 'c')
	{
		for(int i = 0; i < NumberOfNodes; i++)
		{
			NodePosition[i].x -= centerOfMass.x;
			NodePosition[i].y -= centerOfMass.y;
			NodePosition[i].z -= centerOfMass.z;
		}
		drawPicture();
	}
	
	if(key == 'n')
	{
		// Moving to zero
		for(int i = 0; i < NumberOfNodes; i++)
		{
			NodePosition[i].x -= centerOfMass.x;
			NodePosition[i].y -= centerOfMass.y;
			NodePosition[i].z -= centerOfMass.z;
		}
		
		// Rotating until Sinus Node is on x y plane
		dAngle = atan(NodePosition[0].z/NodePosition[0].x);
		for(int i = 0; i < NumberOfNodes; i++)
		{
			temp 			  = cos(dAngle)*(double)NodePosition[i].x + sin(dAngle)*(double)NodePosition[i].z;
			NodePosition[i].z = -sin(dAngle)*(double)NodePosition[i].x + cos(dAngle)*(double)NodePosition[i].z;
			NodePosition[i].x = temp;
		}
		
		
		// Rotating until Sinus Node up on the positive y axis.
		if(0.0 <= NodePosition[0].x)
		{
			dAngle = PI/2.0 - atan(NodePosition[0].y/NodePosition[0].x);
			for(int i = 0; i < NumberOfNodes; i++)
			{
				temp			   = cos(dAngle)*NodePosition[i].x - sin(dAngle)*NodePosition[i].y;
				NodePosition[i].y  = sin(dAngle)*NodePosition[i].x + cos(dAngle)*NodePosition[i].y;
				NodePosition[i].x  = temp;
			}
		}
		else
		{
			dAngle = -(PI/2.0 + atan(NodePosition[0].y/NodePosition[0].x));
			for(int i = 0; i < NumberOfNodes; i++)
			{
				temp			   = cos(dAngle)*NodePosition[i].x - sin(dAngle)*NodePosition[i].y;
				NodePosition[i].y  = sin(dAngle)*NodePosition[i].x + cos(dAngle)*NodePosition[i].y;
				NodePosition[i].x  = temp;
			}
		}
		
		
		drawPicture();
		printf("\n look x y z of 0 %f  %f  %f", NodePosition[0].x, NodePosition[0].y, NodePosition[0].z);
	}
	
}

void mymouse(int button, int state, int x, int y)
{	
	float myX, myY;
	int index = -1;
	
	if(state == GLUT_DOWN)
	{
		if(button == GLUT_LEFT_BUTTON)
		{
			//printf("\n Left mouse button down");
			//printf("\n mouse x = %d mouse y = %d\n", x, y);
			
			myX = 2.0*x/XWindowSize - 1.0;
			myY = -2.0*y/YWindowSize + 1.0;
			
			//printf("\n myX = %f myY = %f\n", myX, myY);
			//printf("\n SNX = %f SNY = %f SNZ = %f\n", NodePosition[0].x, NodePosition[0].y, NodePosition[0].z);
			
			glColor3d(0.0,0.0,1.0);
			glPushMatrix();
				glTranslatef(myX, myY, 0.0);
				glutSolidSphere(0.03,20,20);
			glPopMatrix();
			glutSwapBuffers();
			
			for(int i = 0; i < NumberOfNodes; i++)
			{
				if(myX - 0.03 < NodePosition[i].x && NodePosition[i].x < myX + 0.03 && myY - 0.03 < NodePosition[i].y && NodePosition[i].y < myY + 0.03)
				{
					if(index == -1)
					{
						index = i;
					}
					else
					{
						if(NodePosition[index].z < NodePosition[i].z)
						{
							index = i;
						}
					}
					//drawPicture();
				}
			}
			if(index != -1)
			{
				if(NodeAblatedYesNo[index] == 0)
				{
					NodeAblatedYesNo[index] = 1;
				}
				else
				{
					NodeAblatedYesNo[index] = 0;
				}
			}
		}
		else
		{
			//printf("\nRight mouse button down");
			//printf("\nmouse x = %d mouse y = %d\n", x, y);
			drawPicture();
		}
		//printf("\nSNx = %f SNy = %f SNz = %f\n", NodePosition[0].x, NodePosition[0].y, NodePosition[0].z);
	}
	//drawPicture();
}

void Display(void)
{
	drawPicture();
}

void idle()
{
	n_body(Dt);
}

void reshape(int w, int h)
{
	glViewport(0, 0, (GLsizei) w, (GLsizei) h);
}
