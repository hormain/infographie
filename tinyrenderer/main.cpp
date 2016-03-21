#include "tgaimage.h"
#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <string.h>
#include <cstdlib>
#include <cmath>
#include <unistd.h>
using namespace std;

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const double g_Pi = 3.14159265358979323846;

void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color) { 
  bool steep = std::abs(y1-y0) > std::abs(x1-x0) ;
    if(steep){
      std::swap(x0,y0);
      std::swap(x1,y1);
    }
  if(x1<x0){
    std::swap(x0,x1);
    std::swap(y0,y1);
  }
  int dx = x1-x0;
  int dy = y1-y0;
  int derror = std::abs(dy)*2;
  int error = 0;
  int y = y0;
  for (int x = x0; x<x1; x++){ 
  if(steep){
     image.set(y, x, color); 
  }else{
     image.set(x, y, color); 
  }
  error += derror;
  if(error > dx){
    y += (y1>y0?1:-1);
    error -= 2*dx;
  }
  } 
}
  
void box(float x1, float y1, float x2, float y2, float x3, float y3,float *xmin,float *xmax, float *ymin, float *ymax){
  *xmin=min(min(x1,x2),x3);
  *ymin=min(min(y1,y2),y3);
  *xmax=max(max(x1,x2),x3);
  *ymax=max(max(y1,y2),y3);
}

void produitVec(float x1, float y1, float z1, float x2, float y2, float z2,float* resX, float* resY, float* resZ){
  *resX=y1*z2-z1*y2;
  *resY=z1*x2-x1*z2;
  *resZ=x1*y2-y1*x2;
}

void coordonneVec(float x1, float y1,float z1, float x2, float y2,float z2, float* coorX,float* coorY,float* coorZ){
  *coorX=x2-x1;
  *coorY=y2-y1;
  *coorZ=z2-z1;
}

void coordonneVec2(float x1, float y1, float x2, float y2, float* coorX,float* coorY){
  *coorX=x2-x1;
  *coorY=y2-y1;
}

void produitScal(float x1,float y1, float z1, float lum1, float lum2, float lum3, float *lum){
  float norm=sqrt(x1*x1+y1*y1+z1*z1);
  x1=x1/norm;
  y1=y1/norm;
  z1=z1/norm;
  *lum=x1*lum1+y1*lum2+z1*lum3;
} 

void addMat44(double* mat1, float* mat2,float* res){
   for(int i =0;i<16;i++)res[i]=0;
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){

	res[i*4+j]=mat2[i*4+j]+mat1[i*4+j];
      
    }
  }
}
void normalize2(float x, float y, float z, float* X, float* Y, float *Z){
  float tmp=sqrt(x*x+y*y+z*z);
  *X=(x/tmp)*2.f-1.f;
  *Y=(y/tmp)*2.f-1.f;
  *Z=(z/tmp)*2.f-1.f;
}
void normalize(float x, float y, float z, float* X, float* Y, float *Z){
  float tmp=sqrt(x*x+y*y+z*z);
  *X=(x/tmp);
  *Y=(y/tmp);
  *Z=(z/tmp);
}
void matId4(float * mat){
  for(int i =0;i<16;i++)mat[i]=0;
  mat[0]=1;
  mat[5]=1;
  mat[10]=1;
  mat[15]=1;
}
void matId4(double * mat){
  for(int i =0;i<16;i++)mat[i]=0;
  mat[0]=1;
  mat[5]=1;
  mat[10]=1;
  mat[15]=1;
}
void multMatSimpl4(float* mat1, float* mat2, float* res){
  for(int i =0;i<16;i++)res[i]=0;
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      for(int k=0;k<4;k++){
	res[i*4+j]+=mat2[j+k*4]*mat1[i*4+k];
      }
    }
  }
} 

void multMatSimpl4(double* mat1, double* mat2, double* res){
  for(int i =0;i<16;i++)res[i]=0;
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      for(int k=0;k<4;k++){
	res[i*4+j]+=mat2[j+k*4]*mat1[i*4+k];
      }
    }
  }
} 

void multMatSimpl4(double* mat1, float* mat2, float* res){
  for(int i =0;i<16;i++)res[i]=0;
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      for(int k=0;k<4;k++){
	res[i*4+j]+=mat2[j+k*4]*mat1[i*4+k];
      }
    }
  }
} 

void transformationVec(float *x1, float *y1, float *z1,float *mat,float *X1, float *Y1, float *Z1){
  //calcule de la matrice centrale 3D
  float x_1,y_1,z_1;
  float caj=1.f;//coordonne que l'on rajoute a chaque vecteur pour avoir la 4D
  x_1=*x1*mat[0]+*y1*mat[1]+*z1*mat[2]+caj*mat[3];
  y_1=*x1*mat[4]+*y1*mat[5]+*z1*mat[6]+caj*mat[7];
  z_1=*x1*mat[8]+*y1*mat[9]+*z1*mat[10]+caj*mat[11];
  // cout << x_1 <<"/"<< y_1 <<"/"<< z_1 << endl;
  //projection 4D -> 3D
  float rapport1;
  rapport1=*x1*mat[12]+*y1*mat[13]+*z1*mat[14]+caj*mat[15];
   if(rapport1!=0){
    *X1=x_1/rapport1;
    *Y1=y_1/rapport1;
    *Z1=z_1/rapport1;
    }
   //possible simplification avec mat1 le vecteur
   /*
float mat1[3];
mat1[0]=x1;
mat1[1]=y1;
mat1[2]=z1;
     for(int i=0;i<4;i++){
     for(int j=0;j<4;j++){
     res[i]+=mat[i*4+j]*mat1[j];
     }
     }
    */
}

void coordText(float x1,float y1,float x2,float y2,float x3,float y3,float p1,float p2,float p3,float *u,float *v){
  *u=x1*p1+x2*p2+x3*p3;
  *v=y1*p1+y2*p2+y3*p3;
}

void triangleBarycentre(float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3, TGAImage &image, TGAImage &text,TGAImage &normal,TGAImage &eclairage,float vtx1,float vtx2,float vty1,float vty2,float vtz1,float vtz2,float* tabZ,float Xn1,float Xn2,float Xn3, float Yn1, float Yn2, float Yn3, float Zn1, float Zn2, float Zn3,float lumX,float lumY, float lumZ){
  float xmin=0,xmax=0,ymin=0,ymax=0;
  float undeuxX=0,untroisX=0,undeuxY=0,untroisY=0,undeuxZ=0,untroisZ=0,pointunX=0,pointunY=0;
  float coordX=0,coordY=0,coordZ=0;
  int tailleX=image.get_width();
  float lum2=0;
  float planX=0,planY=0,planZ=0;
 
  float z;
  normalize2(lumX,lumY,lumZ,&lumX,&lumY,&lumZ);
  box(x1,y1,x2,y2,x3,y3,&xmin,&xmax,&ymin,&ymax);
  coordonneVec(x1,y1,z1,x2,y2,z2,&undeuxX,&undeuxY,&undeuxZ);
  coordonneVec(x1,y1,z1,x3,y3,z3,&untroisX,&untroisY,&untroisZ);
   produitVec(undeuxX,undeuxY,undeuxZ,untroisX,untroisY,untroisZ,&planX,&planY,&planZ);
  for(int x=xmin;x<xmax;x++){
    for(int y=ymin;y<ymax;y++){
      float X=0,Y=0,Z=0,u,v;
      coordonneVec2(x,y,x1,y1,&pointunX,&pointunY);
      produitVec(undeuxX,untroisX,pointunX,undeuxY,untroisY,pointunY,&coordX,&coordY,&coordZ);
      X=(1.-(coordX+coordY)/coordZ);
      Y=coordY/coordZ;
      Z=coordX/coordZ;
      z=X*z1+Y*z2+Z*z3;
      coordText(vtx1,vtx2,vty1,vty2,vtz1,vtz2,X,Z,Y,&u,&v);
      float n[3];
      // normal mapping
        TGAColor color2=normal.get(u,v);
      n[0]=color2.r;
      n[1]=color2.g;
      n[2]=color2.b;
      normalize2(n[0],n[1],n[2],&n[0],&n[1],&n[2]);
      produitScal(lumX,lumY,lumZ,n[0],n[1],n[2],&lum2);
      TGAColor spec=eclairage.get(u,v);
      
      //      lum2=abs(lum2);
      if(X>=0 && Y>=0 && Z>=0){
	if(tabZ[tailleX*y+x]<z){
	  TGAColor color=text.get(u,v);
	  color.r=max(min(color.r*lum2,255.f),0.f)+10;
	  color.g=max(min(color.g*lum2,255.f),0.f)+10; 
	  color.b=max(min(color.b*lum2,255.f),0.f)+10;
	  image.set(x,y,color);
	  tabZ[tailleX*y+x]=z;
	}
      }
    }
  }
}

void printmat44(float *m) {
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      std::cout << m[j+i*4] << "\t";
    }
    std::cout << std::endl;
  }
}

void printmat44(double *m) {
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      std::cout << m[j+i*4] << "\t";
    }
    std::cout << std::endl;
  }
}

void printVec(float* v){
  for(int i=0;i<3;i++){
    cout<< v[i] << endl;
  }
}

double rad(double angle){
  return angle*g_Pi/180;
}

double deg(double rad){
  return rad*180/g_Pi;
}
float det(float* mat33){
  return mat33[0]*mat33[4]*mat33[8]+mat33[3]*mat33[7]*mat33[2]+mat33[6]*mat33[5]*mat33[1]-mat33[2]*mat33[4]*mat33[6]-mat33[0]*mat33[5]*mat33[7]-mat33[1]*mat33[4]*mat33[8];
}

void inv33(float* mat,float* res){
  float d=det(mat);
  if(d>0){
  res[0]=(mat[4]*mat[8]-mat[5]*mat[7]);
  res[1]=(mat[5]*mat[6]-mat[3]*mat[8]);
  res[2]=(mat[3]*mat[7]-mat[4]*mat[6]);
  res[3]=(mat[2]*mat[7]-mat[1]*mat[8]);
  res[4]=(mat[0]*mat[8]-mat[2]*mat[6]);
  res[5]=(mat[1]*mat[6]-mat[0]*mat[7]);
  res[6]=(mat[1]*mat[5]-mat[2]*mat[4]);
  res[7]=(mat[2]*mat[3]-mat[5]*mat[0]);
  res[8]=(mat[0]*mat[4]-mat[1]*mat[3]);
  }
}
void matPerspec(float *mat,float c){
  matId4(mat);
  mat[14]=-1/c;
}
void matViewport(int width, int height, int depth, float* mat){
  matId4(mat);
  mat[0]=width;
  mat[3]=width;
  mat[5]=height;
  mat[7]=height;
  mat[10]=depth;
  mat[11]=depth;
}
void matRotX(double angle,double* mat){
  matId4(mat);
  mat[5]=cos(angle);
  mat[6]=-sin(angle);
  mat[9]=sin(angle);
  mat[10]=cos(angle);
}
void matRotY(double angle,double* mat){
  matId4(mat);
  mat[0]=cos(angle);
  mat[2]=-sin(angle);
  mat[8]=sin(angle);
  mat[10]=cos(angle);
}
void matRotZ(double angle,double* mat){
  matId4(mat);
  mat[0]=cos(angle);
  mat[1]=-sin(angle);
  mat[4]=sin(angle);
  mat[5]=cos(angle);
}

void lookat(float cameraX,float cameraY, float cameraZ, float centreX, float centreY, float centreZ, float upX, float upY, float upZ,float * mat) {
  normalize(cameraX,cameraY,cameraZ,&cameraX,&cameraY,&cameraZ);
  normalize(centreX,centreY,centreZ,&centreX,&centreY,&centreZ);
  normalize(upX,upY,upZ,&upX,&upY,&upZ);
  matId4(mat);
  float vecZ[3],vecY[3],vecX[3];
  normalize2(cameraX-centreX,cameraY-centreY,cameraZ-centreZ,&vecZ[0],&vecZ[1],&vecZ[2]);
  produitVec(upX,upY,upZ,vecZ[0],vecZ[1],vecZ[2],&vecX[0],&vecX[1],&vecX[2]);
  normalize2(vecX[0],vecX[1],vecX[2],&vecX[0],&vecX[1],&vecX[2]);
  produitVec(vecZ[0],vecZ[1],vecZ[2],vecX[0],vecX[1],vecX[2],&vecY[0],&vecY[1],&vecY[2]);
  normalize2(vecY[0],vecY[1],vecY[2],&vecY[0],&vecY[1],&vecY[2]);
  float Tr[16],Minv[16];
  matId4(Tr);

    for (int i=0; i<3; i++) {
        Minv[i] = vecX[i];
        Minv[4+i] = vecY[i];
        Minv[8+i] = vecZ[i];

    }
        Tr[3] = -centreX;
	Tr[7] = -centreY;
	Tr[11] = -centreZ;
	multMatSimpl4(Minv,Tr,mat);
}

int main(int argc, char** argv) {
  int reso= 1000;
  TGAImage image(reso,reso, TGAImage::RGB);
  std::ifstream fichier2("obj/african_head.obj",std::ios::in);
  TGAImage texture;
  TGAImage normal;
  TGAImage eclairage;
  normal.read_tga_file("obj/african_head_nm.tga");
  texture.read_tga_file("obj/african_head_diffuse.tga");
  eclairage.read_tga_file("obj/african_head_spec.tga");
  texture.flip_vertically();
  normal.flip_vertically();
  eclairage.flip_vertically();
  int resow=texture.get_width();
  int resoh=texture.get_height();
  float cameraX=0;
  float cameraY=0;
  float cameraZ=5;
  float centreX=0,centreY=-2,centreZ=0;
  float upX=0,upY=1,upZ=0;
  float lumX=5,lumY=-50,lumZ= 50;
  normalize(lumX,lumY,lumZ,&lumX,&lumY,&lumZ);
  float transfo[16];
  float transfo2[16];
  double rot[16];
  double rot2[16];
  double rot3[16];
  double rotTotal[16];
  float camera[16];
  float axe[16];
  float lum[16];
  float orient[16];
  float look[16];
  float lock[16];
  double tmp[16];
  double alpha=0;//-acos(lumY);//y
  double beta=0;//-acos(lumX);//x
  double gama=0;//-acos(lumZ);//z

  cout << alpha <<"/"<< beta<<"/"<< gama<< endl;
  // return 0;
  /*  alpha=rad(alpha);
  beta=rad(beta);
  gama=rad(gama);
  */
  float res[16];
  cout << "debut des matrice" << endl;
  // rotation par rapport a y 
  matRotX(alpha,rot);
  // roation par rapport a x 
  matRotY(beta,rot2);
  // rotation par rapport a z
  matRotZ(gama,rot3);
  // matrice camera
  for(int i =0;i<16;i++)camera[i]=0;
  camera[3]=cameraX;
  camera[7]=cameraY;
  camera[11]=cameraZ;
  // matrice orientation lumiere
  for(int i=0;i<16;i++)lum[i]=0;
  lum[3]=lumX;
  lum[7]=lumY;
  lum[11]=lumZ;
  // matrive viewport
  matViewport(reso/2,reso/2,reso/8,transfo);
  //matrice rot total
  multMatSimpl4(rot,rot2,tmp);
  //printmat44(tmp);
  multMatSimpl4(tmp,rot3,rotTotal);
  //matrice axer camera 
  addMat44(rotTotal,camera,axe);
  addMat44(rotTotal,lum,orient);
  cout << "axe" << endl;
  printmat44(axe);
  cout << "lum" << endl;
  printmat44(orient);
  // matrice perspective
  float dcam=sqrt(cameraX*cameraX+cameraY*cameraY+cameraZ*cameraZ);
  matPerspec(transfo2,dcam);
  // look at 
  lookat(cameraX,cameraY,cameraZ,centreX,centreY,centreZ,upX,upY,upZ,look);
  cout << "fin debut des matrice" << endl ;

  printmat44(transfo);
  std::cout << std::endl;
  printmat44(transfo2);
  std::cout << std::endl;
  printmat44(rot);
  cout << endl;
  printmat44(rot2);
  cout << endl;
  printmat44(rot3);
  cout << "project"<< endl;
  float resi[16];
  multMatSimpl4(look,transfo2,lock);
  multMatSimpl4(orient,transfo2,lum);
  multMatSimpl4(axe,transfo2,resi);
  resi[15]=1;
  printmat44(resi);
  printmat44(orient);
  cout << endl;
  multMatSimpl4(transfo,lock,look);
  multMatSimpl4(transfo,resi,res);
  multMatSimpl4(transfo,lum,orient);
  printmat44(res);
  printmat44(orient);
  printmat44(look);
  std::cout << std::endl;
  int info[4];
  //return 0;
  std::string ligne;
  string id0;
  cout << "fin mult matrice" << endl<< endl;
  int nbpoint=0;
  int nbvectext=0;
  int nbvectnorm=0;
  int nbface=0;
  if(fichier2){
    int q=0;

    while(fichier2.good()){
      fichier2 >> id0;
      if(id0.size()>0){
	char id=id0.at(0);
	char d='#';
	if(id==d){
	  fichier2 >> info[q];
	  q++;
	  getline(fichier2,ligne);
	}
      }
    }
    int nbpoint=info[0];
    cout << "nbpoint : " << nbpoint << endl;
    int nbvectext=info[1];
    cout << "nbvctext : " << nbvectext << endl;
    int nbvectnorm=info[2];
    cout << "nbvectnorm : " << nbvectnorm << endl;
    int nbface=info[3];
    cout << "nbface : " << nbface << endl;
    cout <<"fin de premiere lecture "<< endl << endl;
    fichier2.close();
  }
  std::ifstream fichier("obj/african_head.obj",std::ios::in);
  if(fichier){
    printf("Lecture du fichier : \n");
     
    string p1,p2,p3;
    int x[2492],y[2492],z[2492];
    int vx[2492],vy[2492],vz[2492];
    int vnx[2492],vny[2492],vnz[2492];
    size_t str,str2,str3;
    int i=0;
    int j=0;
    int k=0;
    int l=0;
    float x1[1258],y1[1258],z1[1258];
    float vtx[1339],vty[1339],vtz[1339];
    float vnx1[1258],vny1[1258],vnz1[1258];

    int tailletotal=image.get_width()*image.get_height();
    float tabZ[tailletotal];
   
    for(int i=0;i<tailletotal;i++)tabZ[i]=-5000;
    while( fichier.good()){
      //  cout << i <<"/"<< j <<"/"<< k <<"/"<< l << endl;
      
      
      fichier >> id0;
      if(id0.size()>0){
	char id= id0.at(0);
	char f='f';
	char v='v';
	if(id==v){
	  if(id0.size()==1){
	    fichier >> x1[j];
	    fichier >> y1[j];
	    fichier >> z1[j];	 
	    j++;
	  }else{
	    if(id0.at(1)=='t'){
	      fichier >> vtx[k];
	      fichier >> vty[k];
	      fichier >> vtz[k];
	      k++;
	    }else{
	      if(id0.at(1)=='n'){
		fichier >> vnx1[l];
		fichier >> vny1[l];
		fichier >> vnz1[l];
		l++;		     
	      }else{
		getline(fichier,ligne);
	      }
	    }
	  }
	}
	if(id==f){
	  getline(fichier, ligne);
	  str=ligne.find(' ');
	  str2=ligne.find(' ', str+1);
	  p1=(ligne.substr(str,str2));
	  str3=ligne.find(' ', str2+1);
	  p2=(ligne.substr(str2,str3));
	  p3=(ligne.substr(str3,ligne.length()));
	  str=p1.find('/');		
	  x[i]=atoi((p1.substr(0,str)).c_str());
	  str2=p1.find('/',str+1);
	  vx[i]=atoi(p1.substr(str+1,str2).c_str());
	  vnx[i]=atoi(p1.substr(str2+1,p1.size()).c_str());
	  str=p2.find('/');
	  y[i]=atoi((p2.substr(0,str)).c_str());
	  str2=p2.find('/',str+1);
	  vy[i]=atoi(p2.substr(str+1,str2).c_str());
	  vny[i]=atoi(p2.substr(str2+1,p2.size()).c_str());
	  str=p3.find('/');
	  z[i]=atoi((p3.substr(0,str)).c_str());
	  str2=p3.find('/',str+1);
	  vz[i]=atoi(p3.substr(str+1,str2).c_str());
	  vnz[i]=atoi(p3.substr(str2+1,p3.size()).c_str());
        		
	  float X1,X2,X3;
	  float Y1,Y2,Y3;
	  float Z1,Z2,Z3;
	  
	  X1=x1[x[i]-1];
	  X2=x1[y[i]-1];
	  X3=x1[z[i]-1];
	  Y1=y1[x[i]-1];
	  Y2=y1[y[i]-1];
	  Y3=y1[z[i]-1];
	  Z1=z1[x[i]-1];
	  Z2=z1[y[i]-1];
	  Z3=z1[z[i]-1];
	  float X_1=0,X_2=0,X_3=0,Y_1=0,Y_2=0,Y_3=0,Z_1=0,Z_2=0,Z_3=0;

	  transformationVec(&X1,&Y1,&Z1,res,&X_1,&Y_1,&Z_1);
	  transformationVec(&X2,&Y2,&Z2,res,&X_2,&Y_2,&Z_2);
	  transformationVec(&X3,&Y3,&Z3,res,&X_3,&Y_3,&Z_3);
	  
	  X_1=max(min(X_1,(float)(reso)),0.f);
	  X_2=max(min(X_2,(float)(reso)),0.f);
	  X_3=max(min(X_3,(float)(reso)),0.f);
	  Y_1=max(min(Y_1,(float)(reso)),0.f);
	  Y_2=max(min(Y_2,(float)reso),0.f);
	  Y_3=max(min(Y_3,(float)reso),0.f);
	  Z_1=max(Z_1,0.f);
	  Z_2=max(Z_2,0.f);
	  Z_3=max(Z_3,0.f);
	  
	  float VXu,VXv,VYu,VYv,VZu,VZv;
	        
	  VXu=vtx[vx[i]-1]*resow;
	  VXv=vty[vx[i]-1]*resoh;
	  VYu=vtx[vy[i]-1]*resow;
	  VYv=vty[vy[i]-1]*resoh;
	  VZu=vtx[vz[i]-1]*resow;
	  VZv=vty[vz[i]-1]*resoh;

	  float Xn1,Xn2,Xn3,Yn1,Yn2,Yn3,Zn1,Zn2,Zn3;
	  resow=eclairage.get_width();
	  resoh=eclairage.get_height();
	  Xn1=vnx1[vnx[i]-1];
	  Xn2=vny1[vnx[i]-1];
	  Xn3=vnz1[vnx[i]-1];
	  Yn1=vnx1[vny[i]-1];
	  Yn2=vny1[vny[i]-1];
	  Yn3=vnz1[vnz[i]-1];
	  Zn1=vnx1[vnz[i]-1];
	  Zn2=vny1[vnz[i]-1];
	  Zn3=vnz1[vnz[i]-1];

	  triangleBarycentre(X_1,Y_1,Z_1,X_2,Y_2,Z_2,X_3,Y_3,Z_3,image,texture,normal,eclairage,VXu,VXv,VYu,VYv,VZu,VZv,tabZ,Xn1,Xn2,Xn3,Yn1,Yn2,Yn3,Zn1,Zn2,Zn3,lumX,lumY,lumZ);
	 
	  i++;
	}else{
	  getline(fichier,ligne);
	}
      }
    }

    fichier.close();
    printf("On a fini de lire le fichier\n");
  
  }else{
  }

  image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
  image.write_tga_file("output.tga");
  return 0;
}

