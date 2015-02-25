#include "sgct.h"
#include <sstream>
#include <iostream>

sgct::Engine * gEngine;

//-----------------------
// function declarations 
//-----------------------
void myInitOGLFun();
void myPreSyncFun();
void myDrawFun();
//totte gör en funktion


// Add comment

//for syncing variables across a cluster
void myEncodeFun();
void myDecodeFun();

void drawAxes(float size);
void drawWireCube(float size);

//Tillagt
void drawSphere(bool b);
bool isTouched(glm::vec3 vec);
void drawCubesAroundWand(glm::mat4 T);

//-----------------------
// variable declarations 
//-----------------------

//store each device's transform 4x4 matrix in a shared vector
sgct::SharedVector<glm::mat4> sharedTransforms;
sgct::SharedString sharedText;
sgct::SharedObject<size_t> sharedHeadSensorIndex(0);

//Declare spheres
sgct_utils::SGCTSphere * sun = NULL;
sgct_utils::SGCTSphere * earth = NULL;
sgct_utils::SGCTSphere * moon = NULL;
sgct_utils::SGCTSphere * wand = NULL;

//pointer to a device
sgct::SGCTTrackingDevice * devicePtr = NULL;
//pointer to a tracker
sgct::SGCTTracker * trackerPtr = NULL;

sgct::SharedDouble curr_time(0.0);
int time1 = 0;

const float SPHERE_SIZE = 0.04f;

int main( int argc, char* argv[] )
{
	// Allocate
	gEngine = new sgct::Engine( argc, argv );

	// Bind your functions
	gEngine->setInitOGLFunction( myInitOGLFun );
	gEngine->setPreSyncFunction( myPreSyncFun );
	gEngine->setDrawFunction( myDrawFun );

	// Init the engine
	if( !gEngine->init() )
	{
		delete gEngine;
		return EXIT_FAILURE;
	}

	sgct::SharedData::instance()->setEncodeFunction( myEncodeFun );
	sgct::SharedData::instance()->setDecodeFunction( myDecodeFun );

	// Main loop
	gEngine->render();

	// Clean up
	delete gEngine;

	// Exit program
	exit( EXIT_SUCCESS );
}

void myInitOGLFun()
{
	glEnable(GL_DEPTH_TEST);
	
    sun = new sgct_utils::SGCTSphere(0.4f, 512);
	earth = new sgct_utils::SGCTSphere(0.2f, 512);
	moon = new sgct_utils::SGCTSphere(0.1f, 512);
	wand = new sgct_utils::SGCTSphere(0.08f, 512);
  	
	
	//only store the tracking data on the master node
	if( gEngine->isMaster() )
	{
		size_t index = 0;
		
		//allocate shared data
		for(size_t i = 0; i < sgct::Engine::getTrackingManager()->getNumberOfTrackers(); i++)
		{
			trackerPtr = sgct::Engine::getTrackingManager()->getTrackerPtr(i);
			
			//init the shared vector with identity matrixes
			for(size_t j=0; j<trackerPtr->getNumberOfDevices(); j++)
			{
				devicePtr = trackerPtr->getDevicePtr(j);
			
				if( devicePtr->hasSensor() )
				{
					sharedTransforms.addVal( glm::mat4(1.0f) );
					
					//find the head sensor
					if( sgct::Engine::getTrackingManager()->getHeadDevicePtr() == devicePtr ){
						sharedHeadSensorIndex.setVal(index);
					}

					sgct_utils::SGCTSphere * sun = NULL;
					index++;
				}
			}
		}
	}
}


/*
	This callback is called once per render loop iteration.
*/
void myPreSyncFun()
{
  
	//std::cout << "test";
	/*
	Store all transforms in the array by looping through all trackers and all devices.

	Storing values from the tracker in the pre-sync callback will guarantee
	that the values are equal for all draw calls within the same frame.
	This prevents the application from getting different tracked data for
	left and right draws using a stereoscopic display. Remember to get
	all sensor, button and analog data that will affect the rendering in this stage.
	*/
	
	
	//only store the tracking data on the master node
	if( gEngine->isMaster() )
	{	double speed = 25.0;
		size_t index = 0;
		std::stringstream ss;
		
		/*
			Loop trough all trackers (like intersense IS-900, Microsoft Kinect, PhaseSpace etc.)
		*/
		for(size_t i = 0; i < sgct::Engine::getTrackingManager()->getNumberOfTrackers(); i++)
		{
			trackerPtr = sgct::Engine::getTrackingManager()->getTrackerPtr(i);
		
			
			/*
				Loop trough all tracking devices (like headtracker, wand, stylus etc.)
			*/
			for(size_t j = 0; j < trackerPtr->getNumberOfDevices(); j++)
			{
				devicePtr = trackerPtr->getDevicePtr(j);
				
				ss << "Device " << i <<  "-" << j << ": " << devicePtr->getName() << "\n";
				
				if( devicePtr->hasSensor() )
				{
					sharedTransforms.setValAt( index, devicePtr->getWorldTransform() );
					index++;

					double trackerTime = devicePtr->getTrackerDeltaTime();
					ss << "     Sensor id: " << devicePtr->getSensorId()
						<< " freq: " << (trackerTime <= 0.0 ? 0.0 : 1.0/trackerTime) << " Hz\n";

					ss << "\n     Pos\n"
						<< "          x=" << devicePtr->getPosition().x << "\n"
						<< "          y=" << devicePtr->getPosition().y << "\n"
						<< "          z=" << devicePtr->getPosition().z << "\n";

					ss << "\n     Rot\n"
						<< "          rx=" << devicePtr->getEulerAngles().x << "\n"
						<< "          ry=" << devicePtr->getEulerAngles().y << "\n"
						<< "          rz=" << devicePtr->getEulerAngles().z << "\n";
				}

				if( devicePtr->hasButtons() )
				{
					ss << "\n     Buttons\n";
					
					for(size_t k=0; k < devicePtr->getNumberOfButtons(); k++)
					{
						ss << "          Button " << k << ": " << (devicePtr->getButton(k) ? "pressed" : "released") << "\n";
					}
				}

				if( devicePtr->hasAnalogs() )
				{
					ss << "\n     Analog axes\n";
					
					for(size_t k=0; k < devicePtr->getNumberOfAxes(); k++)
					{
						ss << "          Axis " << k << ": " << devicePtr->getAnalog(k) << "\n";
					}
				}

				ss << "\n";
			}
		}

		//store the string stream into the shared string
		sharedText.setVal( ss.str() );
	}
}

/*
	This callback can be called several times per render loop iteration.
	Using a single viewport in stereo (3D) usually results in refresh rate of 120 Hz.
*/


void myDrawFun()
{
 
	double speed = 25.0;
	bool touched;
	std::stringstream ss;
	
	//create scene transform (animation)
	glm::mat4 scene_mat = glm::translate( glm::mat4(0.5f), glm::vec3( 0.0f, 0.0f, -1.0f) );
	scene_mat = glm::rotate( scene_mat, static_cast<float>( curr_time.getVal() * speed ), glm::vec3(0.0f, -1.0f, 0.0f));
	scene_mat = glm::rotate( scene_mat, static_cast<float>( curr_time.getVal() * (speed/2.0) ), glm::vec3(1.0f, 0.0f, 0.0f));
	
	glm::mat4 MVP = gEngine->getActiveModelViewProjectionMatrix() * scene_mat;
	

	glPushMatrix();
	glRotated(sgct::Engine::getTime() * 25.0, 0.0, 1.0, 0.0);
	glTranslatef(0.8f, 0.0f, 0.0f);
	glColor3f(0.0f,1.0f,0.0f);
	earth->draw();
	glPushMatrix();
	glRotated(sgct::Engine::getTime() * 80.0, 0.0, 1.0, 0.0);
	glTranslatef(0.4f, 0.0f, 0.0f);
	glColor3f(0.0f,1.0f,1.0f);
	moon->draw();
	glPopMatrix();
	glPopMatrix();
	
	glPushMatrix();
	glTranslatef(0.8f, 0.0f, 0.0f);
	glColor3f(1.0f,0.5f,0.0f);
	wand->draw();
	glPopMatrix();


	/*
	glPushMatrix();
	glColor3f(1.0f,1.0f,0.0f);
	sun->draw();
	glPopMatrix();
*/

	//Loop through each wand
	for(size_t i = 0; i < sharedTransforms.getSize(); i++)
	{
		if(i != sharedHeadSensorIndex.getVal()) 
		{

			//Get transformation matrix of the wand
			glm::mat4 transMatrix = sharedTransforms.getValAt(i);

			//Draw some cubes around wand
			drawCubesAroundWand(transMatrix);

			//Pos of the wand
			glm::vec4 wandPos = transMatrix*glm::vec4(0.f, 0.f, 0.f, 1.f);

			glm::vec3 wandPosXYZ(wandPos);
			
			//get direction with (0,0,-1,0) and get the first three elements
			glm::vec4 wandDir =  transMatrix*glm::vec4(0.f, 0.f, -1.f, 0.f);
			
			//Vector that goes from sphere to the line
			glm::vec3 sphereDir = glm::vec3(wandPos[0] + wandDir[0], wandPos[1] + wandDir[1], wandPos[2] + wandDir[2]);

			//Take the dot product
			float scalar = glm::dot(sphereDir, wandPosXYZ);

			//Calculate final vec to compare against the sphere radius
			glm::vec3 finalVec = sphereDir*scalar;

			//Check if we have collision
			touched = isTouched(finalVec);

		}

	}


	//Draw sphere
	drawSphere(touched);

	//draw text
	float textVerticalPos = static_cast<float>(gEngine->getActiveWindowPtr()->getYResolution()) - 100.0f;
	int fontSize = 12;
	
	glColor3f(1.0f, 1.0f, 1.0f);
	sgct_text::print(sgct_text::FontManager::instance()->getFont( "SGCTFont", fontSize ),
		120.0f, textVerticalPos,
		sharedText.getVal().c_str() );
}


void drawSphere(bool b){

	if(b){

		glPushMatrix();
		glColor3f(1.0f,1.0f,1.0f);
		sun->draw();
		glPopMatrix();

	}

	else{

		glPushMatrix();
		glColor3f(1.0f,0.0f,0.0f);
		sun->draw();
		glPopMatrix();
	}

}

bool isTouched(glm::vec3 vec){

	if(glm::length(vec) < 0.1f)
		return true;

	else
		return false;

}


void drawCubesAroundWand(glm::mat4 T){

	glLineWidth(2.0);

	glPushMatrix();

	glMultMatrixf( glm::value_ptr( T ) );

	glColor3f(0.5f,0.5f,0.5f);
	drawWireCube(0.1f);

	drawAxes(0.1f);

	//draw pointer line
	glBegin(GL_LINES);
	glColor3f(1.0f,1.0f,0.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, -5.0f);
	glEnd();

	glPopMatrix();
	
}



void myEncodeFun()
{
	sgct::SharedData::instance()->writeVector( &sharedTransforms );
	sgct::SharedData::instance()->writeString( &sharedText );
	sgct::SharedData::instance()->writeObj( &sharedHeadSensorIndex );
}

void myDecodeFun()
{
	sgct::SharedData::instance()->readVector( &sharedTransforms );
	sgct::SharedData::instance()->readString( &sharedText );
	sgct::SharedData::instance()->readObj( &sharedHeadSensorIndex );
}


void drawAxes(float size)
{
	glLineWidth(2.0);
	glBegin(GL_LINES);

	//x-axis
	glColor3f(1.0f,0.0f,0.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(size, 0.0f, 0.0f);

	//y-axis
	glColor3f(0.0f,1.0f,0.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, size, 0.0f);

	//z-axis
	glColor3f(0.0f,0.0f,1.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, size);

	glEnd();
}


void drawWireCube(float size)
{
	//bottom
	glBegin(GL_LINE_STRIP);
	glVertex3f( -size, -size, -size);
	glVertex3f( size, -size, -size);
	glVertex3f( size, -size, size);
	glVertex3f( -size, -size, size);
	glVertex3f( -size, -size, -size);
	glEnd();

	//top
	glBegin(GL_LINE_STRIP);
	glVertex3f( -size, size, -size);
	glVertex3f( size, size, -size);
	glVertex3f( size, size, size);
	glVertex3f( -size, size, size);
	glVertex3f( -size, size, -size);
	glEnd();

	//sides
	glBegin(GL_LINES);
	glVertex3f( -size, -size, -size);
	glVertex3f( -size, size, -size);

	glVertex3f( size, -size, -size);
	glVertex3f( size, size, -size);

	glVertex3f( size, -size, size);
	glVertex3f( size, size, size);

	glVertex3f( -size, -size, size);
	glVertex3f( -size, size, size);
	glEnd();
}