﻿//  
// Submit by ofek shachar and yonathan gov
// הוגש על ידי אופק שחר ויונתן גוב
//
// Created by Ran Dror on April 2018
// Copyright (c) Ran Dror. All right reserved.
// Code can not be copied and/or distributed without the express permission of author
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> 
#include <stdarg.h>
#include "GL/glut.h" 
#include "glm.h"
#include "MatrixAlgebraLib.h"

//////////////////////////////////////////////////////////////////////
// Grphics Pipeline section
//////////////////////////////////////////////////////////////////////
#define WIN_SIZE 500
#define CAMERA_DISTANCE_FROM_AXIS_CENTER 10

typedef struct{
	GLfloat point3D[4];
	GLfloat normal[4];
	GLfloat point3DeyeCoordinates[4];
	GLfloat NormalEyeCoordinates[4];
	GLfloat pointScreen[4];
	GLfloat PixelValue;
} Vertex;

enum ProjectionTypeEnum{ ORTHOGRAPHIC = 1, PERSPECTIVE };
enum DisplayTypeEnum{ FACE_VERTEXES = 11, FACE_COLOR, LIGHTING_FLAT, LIGHTING_GOURARD, LIGHTING_PHONG };
enum DisplayNormalEnum{ DISPLAY_NORMAL_YES = 21, DISPLAY_NORMAL_NO };

typedef struct{
	GLfloat ModelMinVec[3]; //(left, bottom, near) of a model.
	GLfloat ModelMaxVec[3]; //(right, top, far) of a model.
	GLfloat CameraPos[3];
	GLfloat ModelScale;
	GLfloat ModelTranslateVector[3];
	enum ProjectionTypeEnum ProjectionType;
	enum DisplayTypeEnum DisplayType;
	enum DisplayNormalEnum DisplayNormals;
	GLfloat Lighting_Diffuse;
	GLfloat Lighting_Specular;
	GLfloat Lighting_Ambient;
	GLfloat Lighting_sHininess;
	GLfloat LightPosition[3];
} GuiParamsForYou;

GuiParamsForYou GlobalGuiParamsForYou;

//written for you
void setPixel(GLint x, GLint y, GLfloat r, GLfloat g, GLfloat b);

//you should write
void ModelProcessing();
void VertexProcessing(Vertex *v);
void FaceProcessing(Vertex *v1, Vertex *v2, Vertex *v3, GLfloat FaceColor[3]);
GLfloat LightingEquation(GLfloat point[3], GLfloat PointNormal[3], GLfloat LightPos[3], GLfloat Kd, GLfloat Ks, GLfloat Ka, GLfloat n);
void DrawLineBresenham(GLint x1, GLint y1, GLint x2, GLint y2, GLfloat r, GLfloat g, GLfloat b);
void linearEquation(GLfloat x1, GLfloat y1, GLfloat x2, GLfloat y2, GLfloat linear[3]);
void rectanglePoints(GLfloat p1[3], GLfloat p2[3], GLfloat p3[3], int rectangle[4]);
void barycentricCoordinatesDrawPixel(GLfloat Px, GLfloat Py, GLfloat P1x, GLfloat P1y, GLfloat P2x, GLfloat P2y, GLfloat P3x, GLfloat P3y, GLfloat faceColor[3]);
GLfloat distanceFromLinear(GLfloat linear[3], GLfloat Px, GLfloat Py);
GLfloat coefficientForBarycentricCoordinates(GLfloat Px, GLfloat Py, GLfloat P1x, GLfloat P1y, GLfloat P2x, GLfloat P2y, GLfloat P3x, GLfloat P3y);

GLMmodel *model_ptr;
void ClearColorBuffer();
void DisplayColorBuffer();

void GraphicsPipeline()
{
	static GLuint i;
	static GLMgroup* group;
	static GLMtriangle* triangle;
	Vertex v1, v2, v3;
	GLfloat FaceColor[3];

	//calling ModelProcessing every time refreshing screen
	ModelProcessing();

	//call VertexProcessing for every vertrx
	//and then call FaceProcessing for every face
	group = model_ptr->groups;
	srand(0);
	while (group) {
		for (i = 0; i < group->numtriangles; i++) {
			triangle = &(model_ptr->triangles[group->triangles[i]]);

			MatrixCopy(v1.point3D, &model_ptr->vertices[3 * triangle->vindices[0]], 3);
			v1.point3D[3] = 1;
			MatrixCopy(v1.normal, &model_ptr->normals[3 * triangle->nindices[0]], 3);
			v1.normal[3] = 1;
			VertexProcessing(&v1);

			MatrixCopy(v2.point3D, &model_ptr->vertices[3 * triangle->vindices[1]], 3);
			v2.point3D[3] = 1;
			MatrixCopy(v2.normal, &model_ptr->normals[3 * triangle->nindices[1]], 3);
			v2.normal[3] = 1;
			VertexProcessing(&v2);

			MatrixCopy(v3.point3D, &model_ptr->vertices[3 * triangle->vindices[2]], 3);
			v3.point3D[3] = 1;
			MatrixCopy(v3.normal, &model_ptr->normals[3 * triangle->nindices[2]], 3);
			v3.normal[3] = 1;
			VertexProcessing(&v3);

			FaceColor[0] = (GLfloat)rand() / ((GLfloat)RAND_MAX + 1);
			FaceColor[1] = (GLfloat)rand() / ((GLfloat)RAND_MAX + 1);
			FaceColor[2] = (GLfloat)rand() / ((GLfloat)RAND_MAX + 1);
			FaceProcessing(&v1, &v2, &v3, FaceColor);
		}
		group = group->next;
	}

	DisplayColorBuffer();
}

GLfloat Zbuffer[WIN_SIZE][WIN_SIZE];
GLfloat Mmodeling[16];
GLfloat Mlookat[16];
GLfloat Mprojection[16];
GLfloat Mviewport[16];

void ModelProcessing()
{
	// temp added : 
	GLfloat r=1.0, l=-1.0, t=1.0, b=-1.0, n= CAMERA_DISTANCE_FROM_AXIS_CENTER+1.0,f= CAMERA_DISTANCE_FROM_AXIS_CENTER-1.0;
	GLfloat  tempVector[3],wVector[3], uVector[3], vVector[3], center[3] = { 0,0,0 }, up[3] = { 0,1,0 };
	// end temp added.

	// ex2-3-extra: calculating model scaling and translating transformation matrix
	//////////////////////////////////////////////////////////////////////////////////
	M4x4identity(Mmodeling);



	// ex2-3: calculating translate transformation matrix
	//////////////////////////////////////////////////////////////////////////////////
	Mmodeling[12] = GlobalGuiParamsForYou.ModelTranslateVector[0];
	Mmodeling[13] = GlobalGuiParamsForYou.ModelTranslateVector[1];
	Mmodeling[14] = GlobalGuiParamsForYou.ModelTranslateVector[2];


	// ex2-3: calculating scale transformation matrix
	//////////////////////////////////////////////////////////////////////////////////
	Mmodeling[0] = GlobalGuiParamsForYou.ModelScale;
	Mmodeling[5] = GlobalGuiParamsForYou.ModelScale;
	Mmodeling[10] = GlobalGuiParamsForYou.ModelScale;


	// ex2-4: calculating lookat transformation matrix
	//////////////////////////////////////////////////////////////////////////////////
	M4x4identity(Mlookat);
	Vminus(wVector, GlobalGuiParamsForYou.CameraPos, center, 3); // camera - center
	MatrixCopy(tempVector, wVector, 3);
	VscalarMultiply(wVector, wVector, 1.0 / V3Normalize(tempVector), 3); // W = camera - center / | camera - center|
	V3cross(uVector, up, wVector); // up x W
	MatrixCopy(tempVector, uVector, 3);
	VscalarMultiply(uVector, uVector, 1.0 / V3Normalize(tempVector), 3); // U = up x W / | up x W |
	V3cross(vVector, wVector, uVector); // V= W x U
	Mlookat[2] = wVector[0]; //Wx
	Mlookat[0] = uVector[0]; //Ux
	Mlookat[1] = vVector[0]; //Vx

	Mlookat[6] = wVector[1]; //Wy
	Mlookat[4] = uVector[1]; //Uy
	Mlookat[5] = vVector[1]; //Vy

	Mlookat[10] = wVector[2]; //Wz
	Mlookat[8] = uVector[2];  //Uz
	Mlookat[9] = vVector[2];  //Vz

	Mlookat[12] = -GlobalGuiParamsForYou.CameraPos[0] * Mlookat[0] - GlobalGuiParamsForYou.CameraPos[1] * Mlookat[4] - GlobalGuiParamsForYou.CameraPos[2] * Mlookat[8];
	Mlookat[13] = -GlobalGuiParamsForYou.CameraPos[0] * Mlookat[1] - GlobalGuiParamsForYou.CameraPos[1] * Mlookat[5] - GlobalGuiParamsForYou.CameraPos[2] * Mlookat[9];
	Mlookat[14] = -GlobalGuiParamsForYou.CameraPos[0] * Mlookat[2] - GlobalGuiParamsForYou.CameraPos[1] * Mlookat[6] - GlobalGuiParamsForYou.CameraPos[2] * Mlookat[10];


	// ex2-2: calculating Orthographic or Perspective projection transformation matrix
	//////////////////////////////////////////////////////////////////////////////////
	M4x4identity(Mprojection);
	if (GlobalGuiParamsForYou.ProjectionType == ORTHOGRAPHIC) { // if the projection is Orthographic.

		Mprojection[0] = (GLfloat)2 / (r - l);
		Mprojection[5] = (GLfloat)2 / (t - b);
		Mprojection[10] = (GLfloat)-2 / (f - n);
		Mprojection[12] = (-r - l) / (r - l);
		Mprojection[13] = (-t - b) / (t - b);
		Mprojection[14] = (-f - n) / (f - n);
		Mprojection[15] = (GLfloat)1;
	}
	else { // if the projection is Perspective.
		Mprojection[0] = (((GLfloat)2.0) * n) / (r - l);
		Mprojection[5] = (((GLfloat)2.0) * n) / (t - b);
		Mprojection[8] = (r + l) / (r - l);
		Mprojection[9] = (t + b) / (t - b);
		Mprojection[10] = (-f - n) / (f - n);
		Mprojection[11] = ((GLfloat)(-1.0));
		Mprojection[14] = ((((GLfloat)-2.0) * f * n) / (f - n));
		Mprojection[15] = 0.0;
	}

	// ex2-2: calculating viewport transformation matrix
	//////////////////////////////////////////////////////////////////////////////////
	M4x4identity(Mviewport);
	Mviewport[0] = WIN_SIZE / 2.0;
	Mviewport[5] = WIN_SIZE / 2.0;
	Mviewport[10] = 0.5;
	Mviewport[12] = WIN_SIZE / 2.0;
	Mviewport[13] = WIN_SIZE / 2.0;
	Mviewport[14] = 0.5;
	Mviewport[15] = 1.0;

	// ex3: clearing color and Z-buffer
	//////////////////////////////////////////////////////////////////////////////////
	ClearColorBuffer(); // setting color buffer to background color
	//add here clearing z-buffer
	
}


void VertexProcessing(Vertex *v)
{
	GLfloat point3DafterModelingTrans[4];
	GLfloat temp1[4], temp2[4];
	GLfloat point3D_plusNormal_screen[4];
	GLfloat Mmodeling3x3[9], Mlookat3x3[9];
	//temp added :
	GLfloat tempMatrix[16];
	//end temp added!
	// ex2-3: modeling transformation v->point3D --> point3DafterModelingTrans
	//////////////////////////////////////////////////////////////////////////////////
	M4multiplyV4(point3DafterModelingTrans, Mmodeling, v->point3D);


	// ex2-4: lookat transformation point3DafterModelingTrans --> v->point3DeyeCoordinates
	//////////////////////////////////////////////////////////////////////////////////
	M4multiplyV4(v->point3DeyeCoordinates, Mlookat, point3DafterModelingTrans);


	// ex2-2: transformation from eye coordinates to screen coordinates v->point3DeyeCoordinates --> v->pointScreen
	//////////////////////////////////////////////////////////////////////////////////
	M4multiplyM4(tempMatrix, Mviewport, Mprojection);
	M4multiplyV4(v->pointScreen, tempMatrix, v->point3DeyeCoordinates);


	// ex2-5: transformation normal from object coordinates to eye coordinates v->normal --> v->NormalEyeCoordinates
	//////////////////////////////////////////////////////////////////////////////////
	M3fromM4(Mmodeling3x3, Mmodeling);
	M3fromM4(Mlookat3x3, Mlookat);
	M3multiplyV3(temp1, Mmodeling3x3, v->normal);
	M3multiplyV3(v->NormalEyeCoordinates, Mlookat3x3, temp1);
	V3Normalize(v->NormalEyeCoordinates);
	v->NormalEyeCoordinates[3] = 1;
	
	// ex2-5: drawing normals 
	//////////////////////////////////////////////////////////////////////////////////
	if (GlobalGuiParamsForYou.DisplayNormals == DISPLAY_NORMAL_YES){
		V4HomogeneousDivide(v->point3DeyeCoordinates);
		VscalarMultiply(temp1, v->NormalEyeCoordinates, 0.05, 3);
		Vplus(temp2, v->point3DeyeCoordinates, temp1, 4);
		temp2[3] = 1;
		M4multiplyV4(temp1, Mprojection, temp2);
		M4multiplyV4(point3D_plusNormal_screen, Mviewport, temp1);
		V4HomogeneousDivide(point3D_plusNormal_screen);
		V4HomogeneousDivide(v->pointScreen);
		DrawLineBresenham(round(v->pointScreen[0]), round(v->pointScreen[1]), round(point3D_plusNormal_screen[0]), round(point3D_plusNormal_screen[1]), 0, 0, 1);
	}
	
	// ex3: calculating lighting for vertex
	//////////////////////////////////////////////////////////////////////////////////

}


void FaceProcessing(Vertex *v1, Vertex *v2, Vertex *v3, GLfloat FaceColor[3])
{
	//local var addded:
	int rectangle[4];
	rectanglePoints(v1->pointScreen, v2->pointScreen, v3->pointScreen, rectangle);
	int i, j;
	//end local var added

	V4HomogeneousDivide(v1->pointScreen);
	V4HomogeneousDivide(v2->pointScreen);
	V4HomogeneousDivide(v3->pointScreen);

	if (GlobalGuiParamsForYou.DisplayType == FACE_VERTEXES)
	{

		DrawLineBresenham(round(v1->pointScreen[0]), round(v1->pointScreen[1]), round(v2->pointScreen[0]), round(v2->pointScreen[1]), 1, 1, 1);
		DrawLineBresenham(round(v2->pointScreen[0]), round(v2->pointScreen[1]), round(v3->pointScreen[0]), round(v3->pointScreen[1]), 1, 1, 1);
		DrawLineBresenham(round(v3->pointScreen[0]), round(v3->pointScreen[1]), round(v1->pointScreen[0]), round(v1->pointScreen[1]), 1, 1, 1);
	}
	else{
		//ex3: Barycentric Coordinates and lighting
		//////////////////////////////////////////////////////////////////////////////////
		
		for (i = rectangle[0]; i <= rectangle[1]; i++)
			for (j = rectangle[2]; j <= rectangle[3]; j++)
				//barycentricCoordinatesDrawPixel(i, j, v1->pointScreen[0]/ v1->pointScreen[3], v1->pointScreen[1]/v1->pointScreen[3], v2->pointScreen[0]/v2->pointScreen[3], v2->pointScreen[1] / v2->pointScreen[3], v3->pointScreen[0] / v3->pointScreen[3], v3->pointScreen[1] / v3->pointScreen[3], FaceColor);
		barycentricCoordinatesDrawPixel(i, j, v1->pointScreen[0] , v1->pointScreen[1] , v2->pointScreen[0], v2->pointScreen[1], v3->pointScreen[0], v3->pointScreen[1], FaceColor);




	}
}
void rectanglePoints(GLfloat p1[3], GLfloat p2[3], GLfloat p3[3], int rectangle[4]) {
	GLfloat x[3] = { p1[0],p2[0],p3[0] };
	GLfloat y[3] = { p1[1],p2[1],p3[1] };
	//GLfloat maxX = x[0], minX = x[0], maxY = y[0], minY = y[0];
	GLfloat maxX = 0, minX = WIN_SIZE, maxY = 0, minY = WIN_SIZE;
	int i = 0;

	for (i = 0; i <= 2; i++) {
		if (maxX < x[i]) maxX = x[i];
		if (minX > x[i]) minX = x[i];
		if (maxY < y[i]) maxY = y[i];
		if (minY > y[i]) minY = y[i];
	}
	rectangle[0] = floorf(minX);
	rectangle[1] = ceilf(maxX);
	rectangle[2] = floorf( minY);
	rectangle[3] = ceilf(maxY);
}
void linearEquation(GLfloat x1, GLfloat y1, GLfloat x2, GLfloat y2, GLfloat linear[3]) { // return 0=A , 1=B , 2=C
	GLfloat deltaX, delatY,  m,  n,  A,  B,  C;
	// y=mx+n
	//Ax+By+C=0
	//A=-m
	//B=1
	//C=-n
	// m= (y2 - y1) /(x2 - x1)
	// n = y1-mx1
	deltaX =x2 - x1;
	delatY = y2 - y1;
	m = delatY / deltaX;
	n = y1 - m * x1;
	A = -m;
	B = 1.0;
	C = -n;
	linear[0] = A;
	linear[1] = B;
	linear[2] = C;

}
GLfloat distanceFromLinear(GLfloat linear[3], GLfloat Px, GLfloat Py) { // return the distnace 
	GLfloat distance;
	distance = linear[0] * Px + linear[1] * Py + linear[2];
	//distance /= sqrtf(powf(linear[0], 2) + powf(linear[1], 2));
	return distance;
}
GLfloat coefficientForBarycentricCoordinates(GLfloat Px, GLfloat Py, GLfloat P1x, GLfloat P1y, GLfloat P2x, GLfloat P2y, GLfloat P3x, GLfloat P3y) {
	GLfloat linear[3];
	GLfloat bigDistance, smallDistance, coefficient;
	linearEquation(P2x, P2y, P3x, P3y, linear);
	bigDistance = distanceFromLinear(linear, P1x, P1y);
	smallDistance = distanceFromLinear(linear, Px, Py);
	coefficient = smallDistance / bigDistance;
	return coefficient;
}
void barycentricCoordinatesDrawPixel(GLfloat Px, GLfloat Py, GLfloat P1x, GLfloat P1y, GLfloat P2x, GLfloat P2y, GLfloat P3x, GLfloat P3y, GLfloat faceColor[3]) {
	GLfloat alpha, beta, gamma;
	alpha = coefficientForBarycentricCoordinates(Px, Py, P1x, P1y, P2x, P2y, P3x, P3y);
	beta= coefficientForBarycentricCoordinates(Px, Py, P2x, P2y, P3x, P3y, P1x, P1y);
	gamma = coefficientForBarycentricCoordinates(Px, Py, P3x, P3y, P2x, P2y, P1x, P1y);
	if (alpha>=0 && beta>=0 && gamma>=0 && alpha <=1 && beta<=1 && gamma<=1)
	setPixel(Px, Py, faceColor[0], faceColor[1], faceColor[2]);
}
void DrawLineDDA(GLfloat x1, GLfloat y1, GLfloat x2, GLfloat y2,GLfloat r, GLfloat g, GLfloat b)
{
	float dx, dy, x, y, a, x1_, y1_, x2_, y2_;

	if ((y2 - y1) > -(x2 - x1)) {
		x1_ = x1;
		y1_ = y1;
		x2_ = x2;
		y2_ = y2;
	}
	else
	{
		x1_ = x2;
		y1_ = y2;
		x2_ = x1;
		y2_ = y1;
	}

	dx = x2_ - x1_;
	dy = y2_ - y1_;
	if (fabs(dx) > fabs(dy)) {
		a = dy / dx;
		y = y1_;
		for (x = x1_; x < x2_; x++) {
			setPixel(x, round(y), r, g, b);
			y = y + a;
		}
	}
	else {
		a = dx / dy;
		x = x1_;
		for (y = y1_; y < y2_; y++) {
			setPixel(round(x), y, r, g, b);
			x = x + a;
		}
	}
}

void DrawLineBresenham(GLint x1, GLint y1, GLint x2, GLint y2, GLfloat r, GLfloat g, GLfloat b)
{
	//ex2.1: implement Bresenham line drawing algorithm
	DrawLineDDA(x1, y1, x2, y2,r,g,b);
	//////////////////////////////////////////////////////////////////////////////////
	setPixel(x1, y1, r, g, b);
	setPixel(x2, y2, r, g, b);
}




GLfloat LightingEquation(GLfloat point[3], GLfloat PointNormal[3], GLfloat LightPos[3], GLfloat Kd, GLfloat Ks, GLfloat Ka, GLfloat n)
{
	//ex3: calculate lighting equation
	//////////////////////////////////////////////////////////////////////////////////


	return 1;
}





//////////////////////////////////////////////////////////////////////
// GUI section
//////////////////////////////////////////////////////////////////////

//function declerations
void InitGuiGlobalParams();
void drawingCB(void);
void reshapeCB(int width, int height);
void keyboardCB(unsigned char key, int x, int y);
void keyboardSpecialCB(int key, int x, int y);
void MouseClickCB(int button, int state, int x, int y);
void MouseMotionCB(int x, int y);
void menuCB(int value);
void TerminationErrorFunc(char *ErrorString);
void LoadModelFile();
void DisplayColorBuffer();
void drawstr(char* FontName, int FontSize, GLuint x, GLuint y, char* format, ...);
void TerminationErrorFunc(char *ErrorString);

enum FileNumberEnum{ TEAPOT = 100, TEDDY, PUMPKIN, COW, SIMPLE_PYRAMID, FIRST_EXAMPLE, SIMPLE_3D_EXAMPLE, SPHERE, TRIANGLE, Z_BUFFER_EXAMPLE };

typedef struct{
	enum FileNumberEnum FileNum;
	GLfloat CameraRaduis;
	GLint   CameraAnleHorizontal;
	GLint   CameraAnleVertical;
	GLint   MouseLastPos[2];
} GuiCalculations;

GuiCalculations GlobalGuiCalculations;

GLuint ColorBuffer[WIN_SIZE][WIN_SIZE][3];

int main(int argc, char** argv)
{
	GLint submenu1_id, submenu2_id, submenu3_id, submenu4_id;

	//initizlizing GLUT
	glutInit(&argc, argv);

	//initizlizing GUI globals
	InitGuiGlobalParams();

	//initializing window
	glutInitWindowSize(WIN_SIZE, WIN_SIZE);
	glutInitWindowPosition(900, 100);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
	glutCreateWindow("Computer Graphics HW");

	//registering callbacks
	glutDisplayFunc(drawingCB);
	glutReshapeFunc(reshapeCB);
	glutKeyboardFunc(keyboardCB);
	glutSpecialFunc(keyboardSpecialCB);
	glutMouseFunc(MouseClickCB);
	glutMotionFunc(MouseMotionCB);

	//registering and creating menu
	submenu1_id = glutCreateMenu(menuCB);
		glutAddMenuEntry("open teapot.obj", TEAPOT);
		glutAddMenuEntry("open teddy.obj", TEDDY);
		glutAddMenuEntry("open pumpkin.obj", PUMPKIN);
		glutAddMenuEntry("open cow.obj", COW);
		glutAddMenuEntry("open Simple3Dexample.obj", SIMPLE_3D_EXAMPLE);
		glutAddMenuEntry("open SimplePyramid.obj", SIMPLE_PYRAMID);
		glutAddMenuEntry("open sphere.obj", SPHERE);
		glutAddMenuEntry("open triangle.obj", TRIANGLE);
		glutAddMenuEntry("open FirstExample.obj", FIRST_EXAMPLE);
		glutAddMenuEntry("open ZbufferExample.obj", Z_BUFFER_EXAMPLE);
		submenu2_id = glutCreateMenu(menuCB);
		glutAddMenuEntry("Orthographic", ORTHOGRAPHIC);
		glutAddMenuEntry("Perspective",  PERSPECTIVE);
	submenu3_id = glutCreateMenu(menuCB);
		glutAddMenuEntry("Face Vertexes", FACE_VERTEXES);
		glutAddMenuEntry("Face Color", FACE_COLOR);
		glutAddMenuEntry("Lighting Flat", LIGHTING_FLAT);
		glutAddMenuEntry("Lighting Gourard", LIGHTING_GOURARD);
		glutAddMenuEntry("Lighting Phong", LIGHTING_PHONG);
	submenu4_id = glutCreateMenu(menuCB);
		glutAddMenuEntry("Yes", DISPLAY_NORMAL_YES);
		glutAddMenuEntry("No", DISPLAY_NORMAL_NO);
	glutCreateMenu(menuCB);
		glutAddSubMenu("Open Model File", submenu1_id);
		glutAddSubMenu("Projection Type", submenu2_id);
		glutAddSubMenu("Display type",    submenu3_id);
		glutAddSubMenu("Display Normals", submenu4_id);
	glutAttachMenu(GLUT_RIGHT_BUTTON);

	LoadModelFile();

	//starting main loop
	glutMainLoop();
}

void drawingCB(void)
{
	GLenum er;

	char DisplayString1[200], DisplayString2[200];

	//clearing the background
	glClearColor(0, 0, 0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//initializing modelview transformation matrix
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	GraphicsPipeline();

	glColor3f(0, 1, 0);
	sprintf(DisplayString1, "Scale:%.1f , Translate: (%.1f,%.1f,%.1f), Camera angles:(%d,%d) position:(%.1f,%.1f,%.1f) ", GlobalGuiParamsForYou.ModelScale, GlobalGuiParamsForYou.ModelTranslateVector[0], GlobalGuiParamsForYou.ModelTranslateVector[1], GlobalGuiParamsForYou.ModelTranslateVector[2], GlobalGuiCalculations.CameraAnleHorizontal, GlobalGuiCalculations.CameraAnleVertical, GlobalGuiParamsForYou.CameraPos[0], GlobalGuiParamsForYou.CameraPos[1], GlobalGuiParamsForYou.CameraPos[2]);
	drawstr("helvetica", 12, 15, 25, DisplayString1);
	sprintf(DisplayString2, "Lighting reflection - Diffuse:%1.2f, Specular:%1.2f, Ambient:%1.2f, sHininess:%1.2f", GlobalGuiParamsForYou.Lighting_Diffuse, GlobalGuiParamsForYou.Lighting_Specular, GlobalGuiParamsForYou.Lighting_Ambient, GlobalGuiParamsForYou.Lighting_sHininess);
	drawstr("helvetica", 12, 15, 10, DisplayString2);

	//swapping buffers and displaying
	glutSwapBuffers();

	//check for errors
	er = glGetError();  //get errors. 0 for no error, find the error codes in: https://www.opengl.org/wiki/OpenGL_Error
	if (er) printf("error: %d\n", er);
}


void LoadModelFile()
{
	if (model_ptr){
		glmDelete(model_ptr);
		model_ptr = 0;
	}

	switch (GlobalGuiCalculations.FileNum){
	case TEAPOT:
		model_ptr = glmReadOBJ("teapot.obj");
		break;
	case TEDDY:
		model_ptr = glmReadOBJ("teddy.obj");
		break;
	case PUMPKIN:
		model_ptr = glmReadOBJ("pumpkin.obj");
		break;
	case COW:
		model_ptr = glmReadOBJ("cow.obj");
		break;
	case SIMPLE_PYRAMID:
		model_ptr = glmReadOBJ("SimplePyramid.obj");
		break;
	case FIRST_EXAMPLE:
		model_ptr = glmReadOBJ("FirstExample.obj");
		break;
	case SIMPLE_3D_EXAMPLE:
		model_ptr = glmReadOBJ("Simple3Dexample.obj");
		break;
	case SPHERE:
		model_ptr = glmReadOBJ("sphere.obj");
		break;
	case TRIANGLE:
		model_ptr = glmReadOBJ("triangle.obj");
		break;
	case Z_BUFFER_EXAMPLE:
		model_ptr = glmReadOBJ("ZbufferExample.obj");
		break;
	default:
		TerminationErrorFunc("File number not valid");
		break;
	}

	if (!model_ptr)
		TerminationErrorFunc("can not load 3D model");
	//glmUnitize(model_ptr);  //"unitize" a model by translating it

	//to the origin and scaling it to fit in a unit cube around
	//the origin
	glmFacetNormals(model_ptr);  //adding facet normals
	glmVertexNormals(model_ptr, 90.0);  //adding vertex normals

	glmBoundingBox(model_ptr, GlobalGuiParamsForYou.ModelMinVec, GlobalGuiParamsForYou.ModelMaxVec);
}

void ClearColorBuffer()
{
	GLuint x, y;
	for (y = 0; y <WIN_SIZE; y++){
		for (x = 0; x < WIN_SIZE; x++){
			ColorBuffer[y][x][0] = 0;
			ColorBuffer[y][x][1] = 0;
			ColorBuffer[y][x][2] = 0;
		}
	}
}

void setPixel(GLint x, GLint y, GLfloat r, GLfloat g, GLfloat b)
{
	if (x >= 0 && x < WIN_SIZE && y >= 0 && y < WIN_SIZE){
		ColorBuffer[y][x][0] = round(r * 255);
		ColorBuffer[y][x][1] = round(g * 255);
		ColorBuffer[y][x][2] = round(b * 255);
	}
}

void DisplayColorBuffer()
{
	GLuint x, y;
	glBegin(GL_POINTS);
	for (y = 0; y < WIN_SIZE; y++){
		for (x = 0; x <WIN_SIZE; x++){
			glColor3ub(min(255, ColorBuffer[y][x][0]), min(255, ColorBuffer[y][x][1]), min(255, ColorBuffer[y][x][2]));
				glVertex2f(x + 0.5, y + 0.5);   // The 0.5 is to target pixel
			}
		}
	glEnd();
}


void InitGuiGlobalParams()
{
	GlobalGuiCalculations.FileNum = TEAPOT;
	GlobalGuiCalculations.CameraRaduis = CAMERA_DISTANCE_FROM_AXIS_CENTER;
	GlobalGuiCalculations.CameraAnleHorizontal = 0;
	GlobalGuiCalculations.CameraAnleVertical = 0;

	GlobalGuiParamsForYou.CameraPos[0] = 0;
	GlobalGuiParamsForYou.CameraPos[1] = 0;
	GlobalGuiParamsForYou.CameraPos[2] = GlobalGuiCalculations.CameraRaduis;

	GlobalGuiParamsForYou.ModelScale = 1;

	GlobalGuiParamsForYou.ModelTranslateVector[0] = 0;
	GlobalGuiParamsForYou.ModelTranslateVector[1] = 0;
	GlobalGuiParamsForYou.ModelTranslateVector[2] = 0;
	GlobalGuiParamsForYou.DisplayType = FACE_VERTEXES;
	GlobalGuiParamsForYou.ProjectionType = ORTHOGRAPHIC;
	GlobalGuiParamsForYou.DisplayNormals = DISPLAY_NORMAL_NO;
	GlobalGuiParamsForYou.Lighting_Diffuse   = 0.75;
	GlobalGuiParamsForYou.Lighting_Specular  = 0.2;
	GlobalGuiParamsForYou.Lighting_Ambient   = 0.2;
	GlobalGuiParamsForYou.Lighting_sHininess = 40;
	GlobalGuiParamsForYou.LightPosition[0] = 10;
	GlobalGuiParamsForYou.LightPosition[1] = 5;
	GlobalGuiParamsForYou.LightPosition[2] = 0;

}


void reshapeCB(int width, int height)
{
	if (width != WIN_SIZE || height != WIN_SIZE)
	{
		glutReshapeWindow(WIN_SIZE, WIN_SIZE);
	}

	//update viewport
	glViewport(0, 0, width, height);

	//clear the transformation matrices (load identity)
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	//projection
	gluOrtho2D(0, WIN_SIZE, 0, WIN_SIZE);
}


void keyboardCB(unsigned char key, int x, int y){
	switch (key) {
	case 27:
		exit(0);
		break;
	case '+':
		GlobalGuiParamsForYou.ModelScale += 0.1;
		glutPostRedisplay();
		break;
	case '-':
		GlobalGuiParamsForYou.ModelScale -= 0.1;
		glutPostRedisplay();
		break;
	case 'd':
	case 'D':
		GlobalGuiParamsForYou.Lighting_Diffuse += 0.05;
		glutPostRedisplay();
		break;
	case 'c':
	case 'C':
		GlobalGuiParamsForYou.Lighting_Diffuse -= 0.05;
		glutPostRedisplay();
		break;
	case 's':
	case 'S':
		GlobalGuiParamsForYou.Lighting_Specular += 0.05;
		glutPostRedisplay();
		break;
	case 'x':
	case 'X':
		GlobalGuiParamsForYou.Lighting_Specular -= 0.05;
		glutPostRedisplay();
		break;
	case 'a':
	case 'A':
		GlobalGuiParamsForYou.Lighting_Ambient += 0.05;
		glutPostRedisplay();
		break;
	case 'z':
	case 'Z':
		GlobalGuiParamsForYou.Lighting_Ambient -= 0.05;
		glutPostRedisplay();
		break;
	case 'h':
	case 'H':
		GlobalGuiParamsForYou.Lighting_sHininess += 1;
		glutPostRedisplay();
		break;
	case 'n':
	case 'N':
		GlobalGuiParamsForYou.Lighting_sHininess -= 1;
		glutPostRedisplay();
		break;
	default:
		printf("Key not valid (language shold be english)\n");
	}
}


void keyboardSpecialCB(int key, int x, int y)
{
	switch (key) {
	case GLUT_KEY_LEFT:
		GlobalGuiParamsForYou.ModelTranslateVector[0] -= 0.1;
		glutPostRedisplay();
		break;
	case GLUT_KEY_RIGHT:
		GlobalGuiParamsForYou.ModelTranslateVector[0] += 0.1;
		glutPostRedisplay();
		break;
	case GLUT_KEY_DOWN:
		GlobalGuiParamsForYou.ModelTranslateVector[2] -= 0.1;
		glutPostRedisplay();
		break;
	case GLUT_KEY_UP:
		GlobalGuiParamsForYou.ModelTranslateVector[2] += 0.1;
		glutPostRedisplay();
		break;
	}
}


void MouseClickCB(int button, int state, int x, int y)
{
	GlobalGuiCalculations.MouseLastPos[0] = x;
	GlobalGuiCalculations.MouseLastPos[1] = y;
}

void MouseMotionCB(int x, int y)
{
	GlobalGuiCalculations.CameraAnleHorizontal += (x - GlobalGuiCalculations.MouseLastPos[0])/40;
	GlobalGuiCalculations.CameraAnleVertical   -= (y - GlobalGuiCalculations.MouseLastPos[1])/40;

	if (GlobalGuiCalculations.CameraAnleVertical > 30)
		GlobalGuiCalculations.CameraAnleVertical = 30;
	if (GlobalGuiCalculations.CameraAnleVertical < -30)
		GlobalGuiCalculations.CameraAnleVertical = -30;

	GlobalGuiCalculations.CameraAnleHorizontal = (GlobalGuiCalculations.CameraAnleHorizontal + 360) % 360;
//	GlobalGuiCalculations.CameraAnleVertical   = (GlobalGuiCalculations.CameraAnleVertical   + 360) % 360;

	GlobalGuiParamsForYou.CameraPos[0] = GlobalGuiCalculations.CameraRaduis * sin((float)(GlobalGuiCalculations.CameraAnleVertical+90)*M_PI / 180) * cos((float)(GlobalGuiCalculations.CameraAnleHorizontal+90)*M_PI / 180);
	GlobalGuiParamsForYou.CameraPos[2] = GlobalGuiCalculations.CameraRaduis * sin((float)(GlobalGuiCalculations.CameraAnleVertical+90)*M_PI / 180) * sin((float)(GlobalGuiCalculations.CameraAnleHorizontal+90)*M_PI / 180);
	GlobalGuiParamsForYou.CameraPos[1] = GlobalGuiCalculations.CameraRaduis * cos((float)(GlobalGuiCalculations.CameraAnleVertical+90)*M_PI / 180);
	glutPostRedisplay();
}

void menuCB(int value)
{
	switch (value){
	case ORTHOGRAPHIC:
	case PERSPECTIVE:
		GlobalGuiParamsForYou.ProjectionType = value;
		glutPostRedisplay();
		break;
	case FACE_VERTEXES:
	case FACE_COLOR:
	case LIGHTING_FLAT:
	case LIGHTING_GOURARD:
	case LIGHTING_PHONG:
		GlobalGuiParamsForYou.DisplayType = value;
		glutPostRedisplay();
		break;
	case DISPLAY_NORMAL_YES:
	case DISPLAY_NORMAL_NO:
		GlobalGuiParamsForYou.DisplayNormals = value;
		glutPostRedisplay();
		break;
	default:
		GlobalGuiCalculations.FileNum = value;
		LoadModelFile();
		glutPostRedisplay();
	}
}



void drawstr(char* FontName, int FontSize, GLuint x, GLuint y, char* format, ...)
{
	va_list args;
	char buffer[255], *s;

	GLvoid *font_style = GLUT_BITMAP_TIMES_ROMAN_10;

	font_style = GLUT_BITMAP_HELVETICA_10;
	if (strcmp(FontName, "helvetica") == 0) {
		if (FontSize == 12)
			font_style = GLUT_BITMAP_HELVETICA_12;
		else if (FontSize == 18)
			font_style = GLUT_BITMAP_HELVETICA_18;
	}
	else if (strcmp(FontName, "times roman") == 0) {
		font_style = GLUT_BITMAP_TIMES_ROMAN_10;
		if (FontSize == 24)
			font_style = GLUT_BITMAP_TIMES_ROMAN_24;
	}
	else if (strcmp(FontName, "8x13") == 0) {
		font_style = GLUT_BITMAP_8_BY_13;
	}
	else if (strcmp(FontName, "9x15") == 0) {
		font_style = GLUT_BITMAP_9_BY_15;
	}

	va_start(args, format);
	vsprintf(buffer, format, args);
	va_end(args);

	glRasterPos2i(x, y);
	for (s = buffer; *s; s++)
		glutBitmapCharacter(font_style, *s);
}


void TerminationErrorFunc(char *ErrorString)
{
	char string[256];
	printf(ErrorString);
	fgets(string, 256, stdin);

	exit(0);
}

