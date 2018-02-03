//Imran Adham
//This program is a custom 3D graphics library. It implements most of the OpenGL functions
//to allow you to create basic 3D scenes and render them as JPG images.
//Type "make" in terminal to compile, and "./images" to produce a set of 300 frames
//that animate rotating boxes and moving numbers.
//needs libjpeg version 90 to run.

#include <stdio.h>
#include<iostream>
#include<math.h>
#include <stack>
#include "image.h"
#include<vector>
#include<string>
using namespace std;
static stack<double*> MatStack;
static double* I=new double[16]{1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};
static Image frame(400,400);
static double zvals[400][400];
static int polycount=0;
static int GL_QUADS=1;
static int GL_QUAD_STRIP=2;
static int mode=0;
static int vcount=0;
static bool flag=false;
#define MAXSIZE  10
#define ERROR 	-1

void drawScanLine(int x1, int x2, int y1, double a, double b, double d, double e, bool m){
	int x=x1;
	int z=x2;
	int y=y1;
	if(x1<0)
		x=0;
	if(x1>399)
		x=399;
	if(x2<0)
		z=0;
	if(x2>399)
		z=399;
	if(y1<0)
		y=0;
	if(y1>399)
		y=399;
	if(d!=0){
	if(x==z){
		if(-1*(a*(1.0*x)+b*(399-1.0*y)+e)/d>zvals[x][399-y]){
		frame.setPixel(x,399-y,0,0);
		frame.setPixel(x,399-y,1,0);
		frame.setPixel(x,399-y,2,0);
		zvals[x][399-y]=-1*(a*(1.0*x)+b*(399-1.0*y)+e)/d;
		}
		if(-1*(a*(1.0*x)+b*(399-1.0*y)+e)/d==zvals[x][399-y]){
		frame.setPixel(x,399-y,0,1);
		frame.setPixel(x,399-y,1,1);
		frame.setPixel(x,399-y,2,1);
		}
	}
	else if(x>z){
		for (int i=1;i<x-z;i++){
			if(-1*(a*(1.0*x-1.0*i)+b*(399-1.0*y)+e)/d>zvals[x-i][399-y]){
			frame.setPixel(x-i,399-y,0,1);
			frame.setPixel(x-i,399-y,1,1);
			frame.setPixel(x-i,399-y,2,1);
			zvals[x-i][399-y]=-1*(a*(1.0*x-1.0*i)+b*(399-1.0*y)+e)/d;
			}
		}
		if(-1*(a*(1.0*x)+b*(399-1.0*y)+e)/d>zvals[x][399-y]){
		frame.setPixel(x,399-y,0,0);
		frame.setPixel(x,399-y,1,0);
		frame.setPixel(x,399-y,2,0);
		zvals[x][399-y]=-1*(a*(1.0*x)+b*(399-1.0*y)+e)/d;
		}
		else if(-1*(a*(1.0*x)+b*(399-1.0*y)+e)/d==zvals[x][399-y]){
		frame.setPixel(x,399-y,0,1);
		frame.setPixel(x,399-y,1,1);
		frame.setPixel(x,399-y,2,1);
		}
		if(-1*(a*(1.0*z)+b*(399-1.0*y)+e)/d>zvals[z][399-y]){
		frame.setPixel(z,399-y,0,0);
		frame.setPixel(z,399-y,1,0);
		frame.setPixel(z,399-y,2,0);
		zvals[z][399-y]=-1*(a*(1.0*z)+b*(399-1.0*y)+e)/d;
		}
		else if(-1*(a*(1.0*z)+b*(399-1.0*y)+e)/d==zvals[z][399-y]){
		frame.setPixel(z,399-y,0,1);
		frame.setPixel(z,399-y,1,1);
		frame.setPixel(z,399-y,2,1);
		}
	}
	else if(z>x){
		for (int i=1;i<z-x;i++){
			if(-1*(a*(1.0*x+1.0*i)+b*(399-1.0*y)+e)/d>zvals[x+i][399-y]){
			frame.setPixel(x+i,399-y,0,1);
			frame.setPixel(x+i,399-y,1,1);
			frame.setPixel(x+i,399-y,2,1);
			zvals[x+i][399-y]=-1*(a*(1.0*x+1.0*i)+b*(399-1.0*y)+e)/d;
			}
		}
		if(-1*(a*(1.0*x)+b*(399-1.0*y)+e)/d>zvals[x][399-y]){
		frame.setPixel(x,399-y,0,0);
		frame.setPixel(x,399-y,1,0);
		frame.setPixel(x,399-y,2,0);
		zvals[x][399-y]=-1*(a*(1.0*x)+b*(399-1.0*y)+e)/d;
		}
		else if(-1*(a*(1.0*x)+b*(399-1.0*y)+e)/d==zvals[x][399-y]){
		frame.setPixel(x,399-y,0,1);
		frame.setPixel(x,399-y,1,1);
		frame.setPixel(x,399-y,2,1);
		}
		if(-1*(a*(1.0*z)+b*(399-1.0*y)+e)/d>zvals[z][399-y]){
		frame.setPixel(z,399-y,0,0);
		frame.setPixel(z,399-y,1,0);
		frame.setPixel(z,399-y,2,0);
		zvals[z][399-y]=-1*(a*(1.0*z)+b*(399-1.0*y)+e)/d;
		}
		else if(-1*(a*(1.0*z)+b*(399-1.0*y)+e)/d==zvals[z][399-y]){
		frame.setPixel(z,399-y,0,1);
		frame.setPixel(z,399-y,1,1);
		frame.setPixel(z,399-y,2,1);
		}
	}
	}
	/*else{
		if(x==z){
			frame.setPixel(x,399-y,0,0);
			frame.setPixel(x,399-y,1,0);
			frame.setPixel(x,399-y,2,0);
			}
		else if(x>z){
			for (int i=1;i<x-z;i++){
				frame.setPixel(x-i,399-y,0,0);
				frame.setPixel(x-i,399-y,1,0);
				frame.setPixel(x-i,399-y,2,0);
				}
			frame.setPixel(x,399-y,0,0);
			frame.setPixel(x,399-y,1,0);
			frame.setPixel(x,399-y,2,0);
			frame.setPixel(z,399-y,0,0);
			frame.setPixel(z,399-y,1,0);
			frame.setPixel(z,399-y,2,0);
			}
		else if(z>x){
			for (int i=1;i<z-x;i++){
				frame.setPixel(x+i,399-y,0,1);
				frame.setPixel(x+i,399-y,1,1);
				frame.setPixel(x+i,399-y,2,1);
				}
			frame.setPixel(x,399-y,0,0);
			frame.setPixel(x,399-y,1,0);
			frame.setPixel(x,399-y,2,0);
			frame.setPixel(z,399-y,0,0);
			frame.setPixel(z,399-y,1,0);
			frame.setPixel(z,399-y,2,0);
			}
	}*/
	/*if(m){
		if(x==z){
			frame.setPixel(x,399-y,0,0);
			frame.setPixel(x,399-y,1,0);
			frame.setPixel(x,399-y,2,0);
		}
		else if(x>z){
			for (int i=0;i<x-z+1;i++){
				frame.setPixel(x-i,399-y,0,0);
				frame.setPixel(x-i,399-y,1,0);
				frame.setPixel(x-i,399-y,2,0);
				}
			}
		else if(z>x){
			for (int i=0;i<z-x+1;i++){
				frame.setPixel(x+i,399-y,0,0);
				frame.setPixel(x+i,399-y,1,0);
				frame.setPixel(x+i,399-y,2,0);
				}
			}
	}*/
}

class Point{
public:
	double* p;
};

class Quad{
public:
	vector<Point> points;
	int inter[4];
	int x, y, xmin, ymin, xmax, ymax, c;
	double a, b, d, e;
	void calcVals();
	void display();
	void inters(double);
	void fill(int);
};

static vector< Quad > quads;

void Quad::calcVals(){
	xmin=points[0].p[0];
	xmax=points[0].p[0];
	ymin=points[0].p[1];
	ymax=points[0].p[1];
	for(int i=0;i<4;i++){
		if(xmin>points[i].p[0])
			xmin=points[i].p[0];
		if(xmax<points[i].p[0])
			xmax=points[i].p[0];
		if(ymin>points[i].p[1])
			ymin=points[i].p[1];
		if(ymax<points[i].p[1])
			ymax=points[i].p[1];
	}
	//find plane equation
	a=((points[0].p[1]-points[2].p[1])*(points[1].p[2]-points[0].p[2]))-((points[0].p[2]-points[2].p[2])*(points[1].p[1]-points[0].p[1]));
	b=((points[0].p[2]-points[2].p[2])*(points[1].p[0]-points[0].p[0]))-((points[0].p[0]-points[2].p[0])*(points[1].p[2]-points[0].p[2]));
	d=((points[0].p[0]-points[2].p[0])*(points[1].p[1]-points[0].p[1]))-((points[0].p[1]-points[2].p[1])*(points[1].p[0]-points[0].p[0]));
	e=-a*points[0].p[0]-b*points[0].p[1]-d*points[0].p[2];
}

void Quad::display(){
	calcVals();
	double ytemp=ymin+0.01;
	//cout<<ymin<<" "<<ytemp<<" "<<ymax<<endl;
	while(ytemp<=ymax){
		//cout<<"t"<<endl;
		inters(ytemp);
		fill(ytemp);
		ytemp++;
	}
}

void Quad::inters(double ytemp){
	int x1, x2, y1, y2, temp;
	c=0;
	for(int i=0;i<4;i++){
		x1=points[i].p[0];
		y1=points[i].p[1];
		x2=points[(i+1)%4].p[0];
		y2=points[(i+1)%4].p[1];
		if(y2<y1){
			temp=x1;
			x1=x2;
			x2=temp;
			temp=y1;
			y1=y2;
			y2=temp;
		}
		if(ytemp<=y2 && ytemp>=y1){
			if(y1==y2){
				x=x1;
			}
			else{
				x=((x2-x1)*(ytemp-y1))/(y2-y1);
				x+=x1;
			}
			if(x<=xmax && x>=xmin){
				inter[c++]=x;
			}
		}
	}
}

void Quad::fill(int ytemp){
	for(int i=0; i<c; i+=2){
		if(abs(ytemp-ymax)<=0.2 || abs(ytemp-ymin)<=0.2){
		drawScanLine(inter[i], inter[i+1], ytemp,a,b,d,e,true);
		}
		else
			drawScanLine(inter[i], inter[i+1], ytemp,a,b,d,e,false);
	}
}

//CONVERT VERTICES TO INTS BEFORE DISPLAYING
void printMatrix(double* mat){
	for(int i=0;i<16;i++){
		cout<<mat[i]<<"  ";
		if (i%4==3)
			cout<<endl;
	}
}
//MATRIX OPERATIONS
void multiplication(int row, int col1, int col2, 
double *imat1, double *imat2, double *omat) {
int i, j, k;
double mat[2*MAXSIZE*MAXSIZE];

	for (i = 0; i < row; i++)
	   for ( j = 0; j < col2; j++) {
	      mat[i*col2+j] = 0;
	      for (k = 0; k < col1; k++)
		 mat[i*col2+j] += imat1[i*col1+k]*imat2[k*col2+j];
	   }
	for (i = 0; i < row; i++)
	   for ( j = 0; j < col2; j++)
	      omat[i*col2+j] = mat[i*col2+j];
}
void copymatrix(int row, int column, double *mat1, double *mat2) {
int i, j;

	for (i = 0; i < row; i++)
	   for (j = 0; j < column; j++)
	      *mat1++ = *mat2++;
}
void rotate(char direction, double angle, double *mat) {
#define CUTOFF 1e-10
double cosa, sina;
	
	cosa = cos(angle);
	sina = sin(angle);
	if (fabs(cosa) < CUTOFF) cosa = 0;
	if (fabs(sina) < CUTOFF) sina = 0;
	switch (direction) {
	   case 'x': mat[5] = cosa;
		     mat[6] = sina;
		     mat[9] = -mat[6];
		     mat[10] = mat[5];
		     break;
	   case 'y': mat[0] = cosa;
		     mat[2] = -sina;
		     mat[8] = -mat[2];
		     mat[10] = mat[0];
		     break;
	   case 'z': mat[0] = cosa;
		     mat[1] = sina;
		     mat[4] = -mat[1];
		     mat[5] = mat[0];
		     break;
	   default: break;
	}
}
//END OF MATRIX OPERATIONS


//OPENGL OPERATIONS
void glPushMatrix(){
	double* temp=new double[16];
	copymatrix(4,4,temp,MatStack.top());
	MatStack.push(temp);
}
void glPopMatrix(){
	delete[] MatStack.top();
	MatStack.pop();
}
void glTranslatef(double x, double y, double z){
	double* temp=new double[16]{1,0,0,x,0,1,0,y,0,0,1,z,0,0,0,1};
	double* temp2=new double[16];
	copymatrix(4,4, temp2, MatStack.top());
	multiplication(4,4,4,temp2, temp, MatStack.top());
	delete [] temp;
	delete [] temp2;
}
void glRotatef(double angle, int x, int y, int z){
	double* temp=new double[16]{1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};
	if(x==1)
		rotate('x', angle*3.14159/180, temp);
	if(y==1)
			rotate('y', angle*3.14159/180, temp);
	if(z==1)
			rotate('z', -1*angle*3.14159/180, temp);
	if(x==-1)
			rotate('x', -1*angle*3.14159/180, temp);
	if(y==-1)
			rotate('y', -1*angle*3.14159/180, temp);
	if(z==-1)
			rotate('z', angle*3.14159/180, temp);
	double* temp2=new double[16];
	copymatrix(4,4, temp2, MatStack.top());
	multiplication(4,4,4,temp2, temp, MatStack.top());
	delete [] temp;
	delete [] temp2;
}

void glScalef(double x, double y, double z){
	double* temp=new double[16]{x,0,0,0,0,y,0,0,0,0,z,0,0,0,0,1};
	double* temp2=new double[16]{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	copymatrix(4,4, temp2, MatStack.top());
	multiplication(4,4,4,temp2, temp, MatStack.top());
	delete [] temp;
	delete [] temp2;
}
void glBegin(int type){
	if (type==GL_QUADS){
		Quad temp;
		quads.push_back(temp);
		polycount++;
		mode=0;
	}
	else if(type==GL_QUAD_STRIP){
		Quad temp;
		quads.push_back(temp);
		//polycount++;
		mode=1;
	}
}

void glVertex3f(double x, double y, double z){
	Point temp;
	temp.p=new double[4]{x,y,z,1};
	if(mode==0)
	quads[quads.size()-1].points.push_back(temp);
	if(mode==1){
		if(!flag){
		if(vcount==0 || vcount==1){
		quads[quads.size()-1].points.push_back(temp);
		}
		else if(vcount==2){
			Point temp2;
			temp2.p=new double[4]{0,0,0,1};
			quads[quads.size()-1].points.push_back(temp2);
			quads[quads.size()-1].points.push_back(temp);
		}
		else if(vcount==3){
			quads[quads.size()-1].points[2].p[0]=x;
			quads[quads.size()-1].points[2].p[1]=y;
			quads[quads.size()-1].points[2].p[2]=z;
		}
		vcount++;
		if(vcount==4){
			vcount=0;
			flag=true;
			polycount++;
		}
		}
		else{
			if(vcount==0){
				Quad temp2;
				Point temp3;
				Point temp4;
				temp3.p=new double[4];
				temp4.p=new double[4];
				temp3.p[0]=quads[quads.size()-1].points[3].p[0];
				temp3.p[1]=quads[quads.size()-1].points[3].p[1];
				temp3.p[2]=quads[quads.size()-1].points[3].p[2];
				temp3.p[3]=1;
				quads.push_back(temp2);
				quads[quads.size()-1].points.push_back(temp3);
				temp4.p[0]=quads[quads.size()-2].points[2].p[0];
				temp4.p[1]=quads[quads.size()-2].points[2].p[1];
				temp4.p[2]=quads[quads.size()-2].points[2].p[2];
				temp4.p[3]=1;
				quads[quads.size()-1].points.push_back(temp4);
				Point temp5;
				temp5.p=new double[4]{0,0,0,1};
				quads[quads.size()-1].points.push_back(temp5);
				quads[quads.size()-1].points.push_back(temp);
				vcount++;
			}
			else if(vcount==1){
				quads[quads.size()-1].points[2].p[0]=x;
				quads[quads.size()-1].points[2].p[1]=y;
				quads[quads.size()-1].points[2].p[2]=z;
				vcount=0;
				polycount++;
			}
		}
	}
}
void glEnd(){
	double n=75;
	for(int k=0;k<polycount;k++){
	for(int i=0;i<quads[quads.size()-1-k].points.size();i++){
		double* temp=new double[4];
		for(int j=0;j<4;j++)
			temp[j]=quads[quads.size()-1-k].points[i].p[j];
		multiplication(4,4,1,MatStack.top(),temp,quads[quads.size()-1-k].points[i].p);
		quads[quads.size()-1-k].points[i].p[0]=quads[quads.size()-1-k].points[i].p[0]*
				(100-n/3+n/15*quads[quads.size()-1-k].points[i].p[2])+200;
		quads[quads.size()-1-k].points[i].p[1]=quads[quads.size()-1-k].points[i].p[1]*
				(100-n/3+n/15*quads[quads.size()-1-k].points[i].p[2])+200;
		delete[] temp;
	}
	}
	mode=0;
	flag=false;
	vcount=0;
	polycount=0;
}
//END OF OPENGL OPERATIONS

bool testpixel(int x, int y, vector<double*> points){
	for (int i=0;i<points.size();i++){
	if(0<(points[(i+1)%4][1]-points[i][1])*(x-points[i][0])+(points[(i+1)%4][0]-points[i][0])*(points[i][1]-y))
		return false;
	}
	return true;
}

void clearPoints(){
	for(int i=0;i<quads.size();i++){
		delete[] quads[i].points[0].p;
		delete[] quads[i].points[1].p;
		delete[] quads[i].points[2].p;
		delete[] quads[i].points[3].p;
	}
	quads.clear();
}


void writeFrame(){
	for (int i=0;i<400;i++){
		for (int j=0;j<400;j++){
			for (int k=0;k<3;k++){
				/*if(j%2==0){
				frame.setPixel(i,399-j,k,1);
				int l=0;
				while(l<quads.size() && frame.getPixel(i, 399-j, k)==1){
					if(testpixel(i,j,quads[l]))
						frame.setPixel(i,399-j,k,0);
					l++;
				}
				}
				else{
					frame.setPixel(i,399-j,k,frame.getPixel(i, 400-j, k));
				}*/
				frame.setPixel(i,399-j,k,1);
			}
		}
	}
	for(int i=0;i<quads.size();i++)
		quads[i].display();
	
}


void convertToInt(){
	for (int i=0;i<quads.size();i++){
		for(int j=0;j<4;j++){
			for(int k=0;k<4;k++){
				int temp=quads[i].points[j].p[k];
				quads[i].points[j].p[k]=temp;
			}
		}
	}
}

void setZvals(){
	for(int i=0;i<400;i++){
		for(int j=0;j<400;j++){
			zvals[i][j]=-20.0;
		}
	}
}

void
drawCube(double length, double height, double depth, double r, double g, double b)
{
    glBegin(GL_QUAD_STRIP);
    //front face
    //glColor3f(r,g,b);
    glVertex3f(-1*(length/2.0),-1*(height/2.0),depth/2.0);
    glVertex3f(length/2.0,-1*(height/2.0),depth/2.0);
    glVertex3f(-1*(length/2.0),height/2.0,depth/2.0);
    glVertex3f(length/2.0,height/2.0,depth/2.0);
    //top face
    //glColor3f(.80*r,.80*g,0.80*b);
    glVertex3f(-1*(length/2.0),height/2.0,-1*(depth/2.0));
    glVertex3f(length/2.0,height/2.0,-1*(depth/2.0));
    //back face
    //glColor3f(.60*r,.60*g,0.60*b);
    glVertex3f(-1*(length/2.0),-1*(height/2.0),-1*depth/2.0);
    glVertex3f(length/2.0,-1*(height/2.0),-1*depth/2.0);
    //bottom face
    //glColor3f(.50*r,.50*g,0.50*b);
    glVertex3f(-1*(length/2.0),-1*(height/2.0),depth/2.0);
    glVertex3f(length/2.0,-1*(height/2.0),depth/2.0);
    glEnd();
    //left face
    glBegin(GL_QUADS);
    //glColor3f(.70*r,.70*g,0.70*b);
    glVertex3f(-1*(length/2.0),-1*(height/2.0),depth/2.0);
    glVertex3f(-1*(length/2.0),height/2.0,depth/2.0);
    glVertex3f(-1*(length/2.0),height/2.0,-1*(depth/2.0));
    glVertex3f(-1*(length/2.0),-1*(height/2.0),-1*depth/2.0);
    glEnd();
    //right face
    glBegin(GL_QUADS);
    //glColor3f(.70*r,.70*g,0.70*b);
    glVertex3f((length/2.0),-1*(height/2.0),depth/2.0);
    glVertex3f((length/2.0),height/2.0,depth/2.0);
    glVertex3f((length/2.0),height/2.0,-1*(depth/2.0));
    glVertex3f((length/2.0),-1*(height/2.0),-1*depth/2.0);
    glEnd();

}

void
drawFullCube(double length, double height, double depth, float r, float g, float b, float r2, float g2, float b2)
{
	//main cube
	drawCube(length, height, depth, r, g, b);

	//front top
	glPushMatrix();
	glTranslatef(0,0.5*height,0.5*depth);
	drawCube(length*1.2, height*0.2, depth*0.2, r2, g2, b2);
	glPopMatrix();

	//front bottom
	glPushMatrix();
	glTranslatef(0,-0.5*height,0.5*depth);
	drawCube(length*1.2, height*0.2, depth*0.2, r2, g2, b2);
	glPopMatrix();

	//front right
	glPushMatrix();
	glTranslatef(0.5*length,0,0.5*depth);
	drawCube(length*0.2, height*1.2, depth*0.2, r2, g2, b2);
	glPopMatrix();

	//front left
	glPushMatrix();
	glTranslatef(-0.5*length,0,0.5*depth);
	drawCube(length*0.2, height*1.2, depth*0.2, r2, g2, b2);
	glPopMatrix();

	//right top
	glPushMatrix();
	glTranslatef(0.5*length,0.5*height,0);
	drawCube(length*0.2, height*0.2, depth*1.2, r2, g2, b2);
	glPopMatrix();

	//right bottom
	glPushMatrix();
	glTranslatef(0.5*length,-0.5*height,0);
	drawCube(length*0.2, height*0.2, depth*1.2, r2, g2, b2);
	glPopMatrix();

	//right back
	glPushMatrix();
	glTranslatef(0.5*length,0,-0.5*depth);
	drawCube(length*0.2, height*1.2, depth*0.2, r2, g2, b2);
	glPopMatrix();

	//left top
	glPushMatrix();
	glTranslatef(-0.5*length,0.5*height,0);
	drawCube(length*0.2, height*0.2, depth*1.2, r2, g2, b2);
	glPopMatrix();

	//left bottom
	glPushMatrix();
	glTranslatef(-0.5*length,-0.5*height,0);
	drawCube(length*0.2, height*0.2, depth*1.2, r2, g2, b2);
	glPopMatrix();

	//left back
	glPushMatrix();
	glTranslatef(-0.5*length,0,-0.5*depth);
	drawCube(length*0.2, height*1.2, depth*0.2, r2, g2, b2);
	glPopMatrix();

	//top back
	glPushMatrix();
	glTranslatef(0,0.5*height,-0.5*depth);
	drawCube(length*1.2, height*0.2, depth*0.2, r2, g2, b2);
	glPopMatrix();

	//bottom back
	glPushMatrix();
	glTranslatef(0,-0.5*height,-0.5*depth);
	drawCube(length*1.2, height*0.2, depth*0.2, r2, g2, b2);
	glPopMatrix();

}

void
drawArc(double angle, float r, float g, float b){
	//front face;
	glBegin(GL_QUAD_STRIP);
	for(int i=0;i<angle+1;i++){
		//if(i>1)
		if(i%5==0){
		glVertex3f(0.5*cos((double)i/180*3.14159),0.5*sin((double)i/180*3.14159),0.5);
		glVertex3f(cos((double)i/180*3.14159),sin((double)i/180*3.14159),0.5);
		}
	}
	glEnd();
	//back face
	glBegin(GL_QUAD_STRIP);
	for(int i=0;i<angle+1;i++){
		if(i%5==0){
		glVertex3f(0.5*cos((double)i/180*3.14159),0.5*sin((double)i/180*3.14159),-0.5);
		glVertex3f(cos((double)i/180*3.14159),sin((double)i/180*3.14159),-0.5);
		}
	}
	glEnd();
	//inner ring
	glBegin(GL_QUAD_STRIP);
	for(int i=0;i<angle+1;i++){
		if(i%5==0){
		glVertex3f(0.5*cos((double)i/180*3.14159),0.5*sin((double)i/180*3.14159),-0.5);
		glVertex3f(0.5*cos((double)i/180*3.14159),0.5*sin((double)i/180*3.14159),0.5);
		}
	}
	glEnd();
	//outer ring
	glBegin(GL_QUAD_STRIP);
	for(int i=0;i<angle+1;i++){
		if(i%5==0){
		glVertex3f(cos((double)i/180*3.14159),sin((double)i/180*3.14159),-0.5);
		glVertex3f(cos((double)i/180*3.14159),sin((double)i/180*3.14159),0.5);
		}
	}
	glEnd();

	//first cap
	//glColor3f(0.85*r,0.85*g,0.85*b);
	glBegin(GL_QUADS);;
	glVertex3f(0.5,0,0.5);
	glVertex3f(1,0,0.5);
	glVertex3f(1,0,-0.5);
	glVertex3f(0.5,0,-0.5);
	glEnd();

	//last cap
	//glColor3f(0.85*r,0.85*g,0.85*b);
	glBegin(GL_QUADS);
	glVertex3f((float)0.5*cos((double)angle/180*3.14159),(float)0.5*sin((double)angle/180*3.14159),0.5);
	glVertex3f((float)cos((double)angle/180*3.14159),(float)sin((double)angle/180*3.14159),0.5);
	glVertex3f((float)cos((double)angle/180*3.14159),(float)sin((double)angle/180*3.14159),-0.5);
	glVertex3f((float)0.5*cos((double)angle/180*3.14159),(float)0.5*sin((double)angle/180*3.14159),-0.5);
	glEnd();
}

void drawZero(){
	glPushMatrix();
	//glTranslatef(-1.8,1.5,0);
	glScalef(0.25,0.4,0.1);
	drawArc(360,0.2,1,1);
	glPopMatrix();
}

void drawOne(){
	glPushMatrix();
	//glTranslatef(-0.9,1.5,0);
	drawCube(0.125,0.75,0.1,.1,.8,0.4);
	glPopMatrix();
}

void drawTwo(){
	glPushMatrix();
	//glTranslatef(0,1.55,0);
	glScalef(0.8,0.85,1);
	glPushMatrix();
	glTranslatef(0, 0.1, 0);
	glScalef(0.3,0.3,0.1);
	drawArc(180,1,.4,.3);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(0.017, -0.171, 0);
	glRotatef(325,0,0,1);
	drawCube(0.125,0.77,0.1,1,.4,.3);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(0.035, -0.45, 0);
	drawCube(0.6,0.125,0.125,1,.4,.3);
	glPopMatrix();
	glPopMatrix();
}

void drawThree(){
	glPushMatrix();
	//glTranslatef(0.95,1.5,0);
	glScalef(0.3,0.25,0.1);
	glPushMatrix();
	glTranslatef(0,0.7,0);
	glRotatef(270,0,0,1);
	drawArc(270,.3,.1,0.5);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(0,-0.7,0);
	glRotatef(180,0,0,1);
	drawArc(270,.3,.1,0.5);
	glPopMatrix();
	glPopMatrix();
}

void drawFour(){
	glPushMatrix();
	//glTranslatef(1.8,1.5,0);
	glScalef(0.4,0.45,1);
	drawCube(1,0.2,.1,.6,.7,.8);
	glPushMatrix();
	glTranslatef(-0.4,0.4,0);
	drawCube(0.2,1,0.1,.6,.7,.8);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(0.4,0.4,0);
	drawCube(0.2,1,0.1,.6,.7,.8);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(0.4,-0.4,0);
	drawCube(0.2,1,0.1,.6,.7,.8);
	glPopMatrix();
	glPopMatrix();
}

void drawFive(){
	glPushMatrix();
	//glTranslatef(-1.8,0,0);
	glPushMatrix();
	glScalef(0.3,0.3,0.1);
	glRotatef(205,0,0,1);
	drawArc(320,0.6,0.4,0.5);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(-0.12,0.35,0);
	glRotatef(345,0,0,1);
	drawCube(0.15,0.4,0.1,0.6,0.4,0.5);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(0.06,0.5,0);
	drawCube(0.35,0.15,0.1,0.6,0.4,0.5);
	glPopMatrix();
	glPopMatrix();
}

void drawSix(){
	glPushMatrix();
	//glTranslatef(-0.9,0,0);
	glPushMatrix();
	glScalef(0.3,0.3,0.1);
	glRotatef(205,0,0,1);
	drawArc(360,0.5,0.3,1);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(-0.1,0.3,0);
	glRotatef(335,0,0,1);
	drawCube(0.13,0.75,0.1,0.5,0.3,1);
	glPopMatrix();
	glPopMatrix();
}

void drawSeven(){
	glPushMatrix();
	//glTranslatef(0,0.15,0);
	glPushMatrix();
	glRotatef(-30,0,0,1);
	drawCube(0.2,0.9,0.1,1,0.5,1);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(-0.02,0.39,0);
	drawCube(0.65,0.2,0.1,1,0.5,1);
	glPopMatrix();
	glPopMatrix();
}

void drawEight(){
	glPushMatrix();
	//glTranslatef(0.9,-0.13,0);
	glScalef(1,0.9,1);
	glPushMatrix();
	glScalef(0.3,0.3,0.1);
	drawArc(360,1,1,0.5);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(0,0.45,0);
	glScalef(0.25,0.25,0.1);
	drawArc(360,1,1,0.5);
	glPopMatrix();
	glPopMatrix();
}

void drawNine(){
	glPushMatrix();
	//glTranslatef(1.8,0.25,0);
	glRotatef(180,0,0,1);
	glPushMatrix();
	glScalef(0.3,0.3,0.1);
	glRotatef(205,0,0,1);
	drawArc(360,0.4,1,0.6);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(-0.1,0.3,0);
	glRotatef(335,0,0,1);
	drawCube(0.13,0.75,0.1,0.4,1,0.6);
	glPopMatrix();
	glPopMatrix();
}

void
drawCubeOne(void){
	drawFullCube(.5,.5,.5,.4,.2,.6,0,1,0);

	//draw triangle front
	glPushMatrix();
	glTranslatef(0, -0.03125, 0.3);
	//drawTriangle(1, 1, 1);
	glPopMatrix();

	//draw donut top
	glPushMatrix();
	glTranslatef(0, .25, 0);
	glRotatef(270,1,0,0);
	glScalef(0.15,0.15,0.1);
	drawArc(360, 0.7, 0.5, 0.9);
	glPopMatrix();

	//draw plus sign left
	glPushMatrix();
	glTranslatef(-0.25, 0, 0);
	glRotatef(45,1,0,0);
	glRotatef(270,0,1,0);
	drawCube(0.1, 0.3, 0.1, .4,.6,.4);
	glPushMatrix();
	glRotatef(90, 0, 0, 1);
	drawCube(0.1, 0.3, 0.1, .4,.6,.4);
	glPopMatrix();
	glPopMatrix();

	//draw square bottom
	glPushMatrix();
	glTranslatef(0, -0.25, 0);
	glRotatef(270,1,0,0);
	drawCube(.3,.3,.1,1,.1,.3);
	glPopMatrix();

	//draw minus right
	glPushMatrix();
	glTranslatef(0.25, 0, 0);
	glRotatef(90,0,1,0);
	glPushMatrix();
	glRotatef(90, 0, 0, 1);
	drawCube(0.1, 0.3, 0.1, .6,.6,.3);
	glPopMatrix();
	glPopMatrix();

	//draw vertical back
	glPushMatrix();
	glTranslatef(0, 0, -0.25);
	glRotatef(90,0,1,0);
	drawCube(0.1, 0.3, 0.1, .6,.6,.3);
	glPopMatrix();
}

void
drawCubeTwo(void){
	drawFullCube(.5,.5,.5,.7,.5,.2,1,0,1);
	//draw plus front
	glPushMatrix();
	glTranslatef(0, 0, 0.25);
	drawCube(0.1, 0.3, 0.1, .2,.4,.6);
	glPushMatrix();
	glRotatef(90, 0, 0, 1);
	drawCube(0.1, 0.3, 0.1, .2,.4,.6);
	glPopMatrix();
	glPopMatrix();

	//draw arrow top
	/*glPushMatrix();
	glTranslatef(0, 0.25, 0.05);
	glRotatef(270,1,0,0);
	drawCube(0.1, 0.2, 0.1, .3,.6,.6);
	glPushMatrix();
	glTranslatef(0,0.1,0.05);
	glScalef(0.7,0.7,1);
	drawTriangle(.3,.6,.6);
	glPopMatrix();
	glPopMatrix();*/

	//draw xor back
	glPushMatrix();
	glTranslatef(0, 0, -0.25);
	glPushMatrix();
	glScalef(0.175,0.175,0.1);
	drawArc(360, 0.7, 0.8, 0.4);
	glPopMatrix();
	drawCube(0.03, 0.2, 0.1, .7,.8,.4);
	glPushMatrix();
	glRotatef(90, 0, 0, 1);
	drawCube(0.03, 0.2, 0.1, .7,.8,.4);
	glPopMatrix();
	glPopMatrix();

	//draw fat diamond bottom
	glPushMatrix();
	glTranslatef(0, -0.25, 0);
	glRotatef(45,0,1,0);
	glRotatef(270,1,0,0);
	drawCube(.25,.25,.1,0.9,.4,.6);
	glPopMatrix();

	//draw L left
	glPushMatrix();
	glTranslatef(-0.25, 0, 0);
	glRotatef(270,0,1,0);
	glPushMatrix();
	glTranslatef(0,-0.1,0);
	drawCube(0.3, 0.1, 0.1, .9,.6,.6);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(-0.1,0,0);
	drawCube(0.1, 0.3, 0.1, .9,.6,.6);
	glPopMatrix();
	glPopMatrix();

	//draw hollow square right
	glPushMatrix();
	glTranslatef(0.25, 0, 0);
	glRotatef(90,0,1,0);
	glPushMatrix();
	glTranslatef(0,-0.1,0);
	drawCube(0.3, 0.1, 0.1, .6,.9,.6);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(-0.1,-0.005,0);
	drawCube(0.1, 0.28, 0.1, .6,.9,.6);
	glPopMatrix();
	glPushMatrix();
	glRotatef(180,0,0,1);
	glPushMatrix();
	glTranslatef(0,-0.1,0);
	drawCube(0.3, 0.1, 0.1, .6,.9,.6);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(-0.1,0,0);
	drawCube(0.1, 0.3, 0.1, .6,.9,.6);
	glPopMatrix();
	glPopMatrix();
	glPopMatrix();

}

void
drawCubeThree(void){
	drawFullCube(.5,.5,.5,.7,.8,.6,0,0,1);

	//draw C front
	glPushMatrix();
	glTranslatef(0, 0, 0.25);
	glScalef(0.1,0.15,0.1);
	glRotatef(30,0,0,1);
	drawArc(300,0.7,0.3,0.4);
	glPopMatrix();

	//draw T top
	glPushMatrix();
	glTranslatef(0, 0.25, 0);
	glRotatef(270,1,0,0);
	glPushMatrix();
	glTranslatef(0,0.1,0);
	drawCube(0.3, 0.1, 0.1, .9,.7,.8);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(0,0,0);
	drawCube(0.1, 0.3, 0.1, .9,.7,.8);
	glPopMatrix();
	glPopMatrix();

	//draw D left
	glPushMatrix();
	glTranslatef(-0.25, 0, -0.03);
	glRotatef(270,0,1,0);
	glPushMatrix();
	glScalef(0.15,0.15,0.1);
	glRotatef(270,0,0,1);
	drawArc(180,0.7,0.9,0.3);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(-0.03,0,0);
	drawCube(0.05,0.3,0.1,0.7,0.9,0.3);
	glPopMatrix();
	glPopMatrix();

	//draw P right
	glPushMatrix();
	glTranslatef(0.25, 0, 0.03);
	glRotatef(90,0,1,0);
	glPushMatrix();
	glTranslatef(0, 0.07, 0);
	glScalef(0.15,0.075,0.1);
	glRotatef(270,0,0,1);
	drawArc(180,0.3,0.9,0.7);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(-0.03,0,0);
	drawCube(0.05,0.3,0.1,0.3,0.9,0.7);
	glPopMatrix();
	glPopMatrix();

	//draw B back
	glPushMatrix();
	glTranslatef(0.05, 0, -0.25);
	glRotatef(180,0,1,0);
	glPushMatrix();
	glTranslatef(0, 0.07, 0);
	glScalef(0.15,0.075,0.1);
	glRotatef(270,0,0,1);
	drawArc(180,0.3,0.4,0.7);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(0, -0.07, 0);
	glScalef(0.15,0.075,0.1);
	glRotatef(270,0,0,1);
	drawArc(180,0.3,0.4,0.7);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(-0.03,0,0);
	drawCube(0.08,0.3,0.1,0.3,0.4,0.7);
	glPopMatrix();
	glPopMatrix();

	//draw S bottom
	glPushMatrix();
	glTranslatef(0,-0.25,0);
	glRotatef(90, 1, 0, 0);
	glPushMatrix();
	glTranslatef(0,0.06,0);
	glRotatef(90,0,0,1);
	glScalef(0.075,0.075,0.1);
	drawArc(180,0.9,0.9,0.5);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(0,-0.06,0);
	glScalef(0.075,0.075,0.1);
	glRotatef(270,0,0,1);
	drawArc(180,0.9,0.9,0.5);
	glPopMatrix();
	glPopMatrix();
}

void
drawCubeFour(void){
	drawFullCube(.5,.5,.5,.6,.8,.2,1,0,0);

	//draw diamond front
	/*glPushMatrix();
	glTranslatef(0,0.05,0.3);
	glPushMatrix();
	glTranslatef(0, 0, 0);
	glScalef(0.7,0.7,1);
	drawTriangle(0.6, 0.4, 0.3);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(0, -0.08, 0);
	glRotatef(180,0,0,1);
	glScalef(0.7,0.7,1);
	drawTriangle(0.6, 0.4, 0.3);
	glPopMatrix();
	glPopMatrix();*/

	//draw equal sign top
	glPushMatrix();
	glTranslatef(0,0.25,0);
	glRotatef(270,1,0,0);
	glPushMatrix();
	glTranslatef(-0.08,0,0);
	drawCube(0.08,0.3,0.1,0.4,0.2,0.8);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(0.08,0,0);
	drawCube(0.08,0.3,0.1,0.4,0.2,0.8);
	glPopMatrix();
	glPopMatrix();

	//draw U left
	glPushMatrix();
	glTranslatef(-0.25,0,0);
	glRotatef(270,0,1,0);
	glPushMatrix();
	glTranslatef(-0.08,0.03,0);
	drawCube(0.05,0.2,0.1,0.7,0.3,0.8);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(0.08,0.03,0);
	drawCube(0.05,0.2,0.1,0.7,0.3,0.8);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(0,-0.07,0);
	glRotatef(180,0,0,1);
	glScalef(0.1,0.08,0.1);
	drawArc(180,0.7,0.3,0.8);
	glPopMatrix();
	glPopMatrix();

	//draw E right
	glPushMatrix();
	glTranslatef(0.25,0.02,0);
	glRotatef(90,0,1,0);
	glPushMatrix();
	glTranslatef(0,0.1,0);
	drawCube(0.3,0.06,0.1,0.6,0.4,0.7);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(0,-0.02,0);
	drawCube(0.3,0.06,0.1,0.6,0.4,0.7);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(0,-0.14,0);
	drawCube(0.3,0.06,0.1,0.6,0.4,0.7);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(-0.12,-0.01,0);
	drawCube(0.06,0.28,0.1,0.6,0.4,0.7);
	glPopMatrix();
	glPopMatrix();

	//draw F bottom
	glPushMatrix();
	glTranslatef(0,-0.25,0.03);
	glRotatef(90,1,0,0);
	glPushMatrix();
	glTranslatef(0,0.1,0);
	drawCube(0.3,0.06,0.1,0.8,0.4,0.65);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(0,-0.02,0);
	drawCube(0.3,0.06,0.1,0.8,0.4,0.65);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(-0.12,-0.02,0);
	drawCube(0.06,0.3,0.1,0.8,0.4,0.65);
	glPopMatrix();
	glPopMatrix();

	//draw H back
	glPushMatrix();
	glTranslatef(0,0.02,-0.25);
	glRotatef(180,0,1,0);
	glPushMatrix();
	glTranslatef(0,-0.02,0);
	drawCube(0.3,0.06,0.1,0.8,0.5,0.9);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(-0.12,-0.01,0);
	drawCube(0.06,0.28,0.1,0.8,0.5,0.9);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(0.12,-0.01,0);
	drawCube(0.06,0.28,0.1,0.8,0.5,0.9);
	glPopMatrix();
	glPopMatrix();
}

void
drawCubeFive(void){
	drawFullCube(.5,.5,.5,.1,.5,.9,1,1,0);

	//draw cancel sign front
	glPushMatrix();
	glTranslatef(0,0,0.25);
	glPushMatrix();
	glScalef(0.17,0.17,0.1);
	drawArc(360,1,0,0);
	glPopMatrix();
	glPushMatrix();
	glRotatef(45,0,0,1);
	drawCube(0.05,0.3,0.1,1,0,0);
	glPopMatrix();
	glPopMatrix();

	//draw V top
	glPushMatrix();
	glTranslatef(0,0.25,0);
	glRotatef(270,1,0,0);
	glPushMatrix();
	glTranslatef(-.05,0,0);
	glRotatef(15,0,0,1);
	drawCube(.05,.3,.1,.7,.8,.6);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(.05,0,0);
	glRotatef(-15,0,0,1);
	drawCube(.05,.3,.1,.7,.8,.6);
	glPopMatrix();
	glPopMatrix();

	//draw W left
	glPushMatrix();
	glTranslatef(-0.25,0,0);
	glRotatef(270,0,1,0);
	glPushMatrix();
	glTranslatef(0.065,0,0);
	glScalef(0.6,1,1);
	glPushMatrix();
	glTranslatef(-.05,0,0);
	glRotatef(15,0,0,1);
	drawCube(.05,.3,.1,.6,.8,.8);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(.05,0,0);
	glRotatef(-15,0,0,1);
	drawCube(.05,.3,.1,.6,.8,.8);
	glPopMatrix();
	glPopMatrix();
	glPushMatrix();
	glTranslatef(-0.065,0,0);
	glScalef(0.6,1,1);
	glPushMatrix();
	glTranslatef(-.05,0,0);
	glRotatef(15,0,0,1);
	drawCube(.05,.3,.1,.6,.8,.8);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(.05,0,0);
	glRotatef(-15,0,0,1);
	drawCube(.05,.3,.1,.6,.8,.8);
	glPopMatrix();
	glPopMatrix();
	glPopMatrix();

	//draw Q right
	glPushMatrix();
	glTranslatef(0.25,0,0);
	glPushMatrix();
	glRotatef(90,0,1,0);
	glPushMatrix();
	glScalef(0.15,0.15,0.1);
	drawArc(360,0.5,0.8,0.2);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(0.07,-0.07,0);
	glRotatef(45,0,0,1);
	drawCube(0.05,0.15,0.1,0.5,0.8,0);
	glPopMatrix();
	glPopMatrix();
	glPopMatrix();

	//draw star thing bottom
	glPushMatrix();
	glTranslatef(0,-0.25,0);
	glRotatef(90,1,0,0);
	glPushMatrix();
	glRotatef(45,0,0,1);
	drawCube(0.17,0.17,0.1,0.9,0.7,0.8);
	glPopMatrix();
	glPushMatrix();
	drawCube(0.17,0.17,0.1,0.9,0.7,0.8);
	glPopMatrix();
	glPopMatrix();

	//draw triangle star back
	/*glPushMatrix();
	glTranslatef(0,-0.02,-0.3);
	glRotatef(180,0,1,0);
	glPushMatrix();
	glScalef(0.8,0.8,1);
	drawTriangle(0.4,0.9,1);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(0,0.03,0);
	glRotatef(180,0,0,1);
	glScalef(0.8,0.8,1);
	drawTriangle(0.4,0.9,1);
	glPopMatrix();
	glPopMatrix();*/
}

int main(int argc, char **argv)
{
	cout<<"Drawing Frames."<<endl;
	for(int i=1;i<301;i++){
		setZvals();
		double* T=new double[16]{1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};
		int tick=i;
		MatStack.push(T);
		 glPushMatrix();
		  //draw cubes

		  //center
		  glPushMatrix();
		  glTranslatef(0,1.5,-1.5);
		  //glTranslatef(0,-4,0);
		  glRotatef(1.5*tick,0,0,1);
		  glRotatef(-3*tick,0,1,0);
		  //drawFullCube(0.5,0.5,0.5,0,0,0,0,0,0);
		  drawCubeThree();
		  glPopMatrix();

		  //most bottom left
		  glPushMatrix();
		  glTranslatef(-1.5,-1.5,0.5);
		  //glTranslatef(0,-4,0);
		  glRotatef(2*tick,1,0,0);
		  glRotatef(-3.5*tick,0,1,0);
		  glRotatef(40,0,1,0);
		  //drawFullCube(0.5,0.5,0.5,0,0,0,0,0,0);
		  drawCubeOne();
		  glPopMatrix();

		  //bottom right
		  glPushMatrix();
		  glTranslatef(0.8,0,-0.5);
		  //glTranslatef(0,-4,0);
		  glRotatef(2*tick,0,0,1);
		  glRotatef(2.5*tick,1,0,0);
		  glRotatef(-20,0,1,0);
		  //drawFullCube(0.5,0.5,0.5,0,0,0,0,0,0);
		  drawCubeFour();
		  glPopMatrix();

		  //bottom left
		  glPushMatrix();
		  glTranslatef(-0.8,0,-0.5);
		  //glTranslatef(0,-4,0);
		  glRotatef(-2*tick,1,0,0);
		  glRotatef(3*tick,0,0,1);
		  glRotatef(20,0,1,0);
		  //drawFullCube(0.5,0.5,0.5,0,0,0,0,0,0);
		  drawCubeTwo();
		  glPopMatrix();

		  //most bottom right
		  glPushMatrix();
		  glTranslatef(1.5,-1.5,0.5);
		  //glTranslatef(0,-4,0);
		  glRotatef(2.5*tick,0,0,1);
		  glRotatef(-3.5*tick,0,1,0);
		  glRotatef(-40,0,1,0);
		  //drawFullCube(0.5,0.5,0.5,0,0,0,0,0,0);
		  drawCubeFive();
		  glPopMatrix();
		  //glTranslatef(0,-4,0);
		  //drawNumbers();
		 if(tick<36){
		  glPushMatrix();
		  glTranslatef(-4+0.2*tick,-1,0);
		  drawZero();
		  glPopMatrix();
	}
	if(tick<45){
		  glPushMatrix();
		  glTranslatef(-6+0.2*tick,-1,0);
		  drawOne();
		  glPopMatrix();
	}
		  if(tick>25 && tick<55){
		  glPushMatrix();
		  glTranslatef(-8+0.2*tick,-1,0);
		  drawTwo();
		  glPopMatrix();
		  }
		  if(tick>35 && tick<100){
		  glPushMatrix();
		  glTranslatef(-10+0.2*tick,-1,0);
		  drawThree();
		  glPopMatrix();
		  }
		  if(tick>45 && tick<75){
		  glPushMatrix();
		  glTranslatef(-12+0.2*tick,-1,0);
		  drawFour();
		  glPopMatrix();
		  }
		  if(tick>55 && tick<85){
		  glPushMatrix();
		  glTranslatef(-14+0.2*tick,-1,0);
		  drawFive();
		  glPopMatrix();
		  }
		  if(tick>65 && tick<95){
		  glPushMatrix();
		  glTranslatef(-16+0.2*tick,-1,0);
		  drawSix();
		  glPopMatrix();
		  }
		  if(tick>75 && tick<105){
		  glPushMatrix();
		  glTranslatef(-18+0.2*tick,-1,0);
		  drawSeven();
		  glPopMatrix();
		  }
		  if(tick>85 && tick<115){
		  glPushMatrix();
		  glTranslatef(-20+0.2*tick,-1.15,0);
		  drawEight();
		  glPopMatrix();
		  }
		  if(tick>95 && tick<125){
		  glPushMatrix();
		  glTranslatef(-22+0.2*tick,-0.8,0);
		  drawNine();
		  glPopMatrix();
		  }
		  glPopMatrix();
	
	//convertToInt();
	writeFrame();
	/*cout<<zvals[0][0]<<endl;
	cout<<-1*(quads[0].a*(1.0*0)+quads[0].b*(399-1.0*1.5)+quads[0].e)/quads[0].d<<endl;
	cout<<quads[0].a<<endl;
	cout<<quads[0].b<<endl;
	cout<<quads[0].d<<endl;
	cout<<quads[0].e<<endl;*/
	clearPoints();
	//cout<<"Frame "+std::to_string(i)<<endl;
	delete[] MatStack.top();
	MatStack.pop();
	//drawScanLine(10,10+2*i,50);
	if(i<10){
		frame.writeJPG(("000"+std::to_string(i)).c_str());
	}
	else if(i<99){
		frame.writeJPG(("00"+std::to_string(i)).c_str());
	}
	else if(i<999){
		frame.writeJPG(("0"+std::to_string(i)).c_str());
	}
	else{
		frame.writeJPG((std::to_string(i)).c_str());
	}
	}
	cout<<"Done."<<endl;
  return 0;             /* ANSI C requires main to return int. */
}
