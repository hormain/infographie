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
  *xmin=min(x1,x2);
  *xmin=min(*xmin,x3);
  *ymin=min(y1,y2);
  *ymin=min(*ymin,y3);
  *xmax=max(x1,x2);
  *xmax=max(*xmax,x3);
  *ymax=max(y1,y2);
  *ymax=max(*ymax,y3);

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
void normalize(float x, float y, float z, float* X, float* Y, float *Z){
  float tmp=sqrt(x*x+y*y+z*z);
  *X=x/tmp;
  *Y=y/tmp;
  *Z=z/tmp;
}
void matId4(float * mat){
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
  //  cout << *x1 <<"/" << *y1 << endl;
  //cout << "transformation" << endl;
  /* for(int i=0;i<16;i++){
     cout<<mat[i]<<endl;
     }*/
  x_1=*x1*mat[0]+*y1*mat[1]+*z1*mat[2]+caj*mat[3];
  y_1=*x1*mat[4]+*y1*mat[5]+*z1*mat[6]+caj*mat[7];
  z_1=*x1*mat[8]+*y1*mat[9]+*z1*mat[10]+caj*mat[11];
  // cout << x_1 <<"/"<< y_1 <<"/"<< z_1 << endl;
  //projection 4D -> 3D
  float rapport1, rapport2, rapport3;
  //cout << m4_1 << "/" << m4_2 << "/" << m4_3 << "/" << m4_4 << endl;
  rapport1=*x1*mat[12]+*y1*mat[13]+*z1*mat[14]+caj*mat[15];
  // cout << rapport1 << endl;
   if(rapport1!=0){
    *X1=x_1/rapport1;
    *Y1=y_1/rapport1;
    *Z1=z_1/rapport1;
    }
   /*
     for(int i=0;i<4;i++){
     for(int j=0;j<4;j++){
     res[i]+=mat[i*4+j]*mat1[j];
     }
     }
    */
  // cout << *X1 << "/" << *Y1 <<"/" << *Z1 << endl;
}
void triangleBarycentre(float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3, TGAImage &image, TGAImage &text,TGAImage &normal,TGAImage &eclairage,float vtx1,float vtx2,float vty1,float vty2,float vtz1,float vtz2,float* tabZ, float* transfo,float Xn1,float Xn2,float Xn3, float Yn1, float Yn2, float Yn3, float Zn1, float Zn2, float Zn3){
  float xmin=0,xmax=0,ymin=0,ymax=0;
  float undeuxX=0,untroisX=0,undeuxY=0,untroisY=0,undeuxZ=0,untroisZ=0,pointunX=0,pointunY=0,pointunZ=0;
  float coordX=0,coordY=0,coordZ=0;
  int tailleX=image.get_width();
  //float tabZ[tailleX*image.get_height()];
  float lum=0,lum2=0;
  float planX=0,planY=0,planZ=0;
  float lumX=10,lumY=10,lumZ=10;
  float z;
  float x_1,y_1,z_1,x_2,y_2,z_2,x_3,y_3,z_3;
  normalize(lumX,lumY,lumZ,&lumX,&lumY,&lumZ);
  //boxe 
  box(x1,y1,x2,y2,x3,y3,&xmin,&xmax,&ymin,&ymax);
  
  coordonneVec(x1,y1,z1,x2,y2,z2,&undeuxX,&undeuxY,&undeuxZ);
  coordonneVec(x1,y1,z1,x3,y3,z3,&untroisX,&untroisY,&untroisZ);
   produitVec(undeuxX,undeuxY,undeuxZ,untroisX,untroisY,untroisZ,&planX,&planY,&planZ);
  //produitScal(planX,planY,planZ,lumX,lumY,lumZ,&lum);
  
   //  lum=abs(lum);
  // cout << xmin << '/' << xmax << '/' << ymin << '/' << ymax <<'/'<< lum << endl;
  for(int x=xmin;x<xmax;x++){
    for(int y=ymin;y<ymax;y++){
      float X=0,Y=0,Z=0;
      // cout << x << "/" << y << endl;

      coordonneVec2(x,y,x1,y1,&pointunX,&pointunY);
      //coodonne des vecteur des absices et ordonne des points du triangle
      //prmier vecteur composer de(undeuxX,untroisX,pointunX) et deuxieme vecteur composer de (undeuxY,untroisY,pointunY)
      produitVec(undeuxX,untroisX,pointunX,undeuxY,untroisY,pointunY,&coordX,&coordY,&coordZ);
      X=(1.-(coordX+coordY)/coordZ);
      Y=coordY/coordZ;
      Z=coordX/coordZ;
      z=X*z1+Y*z2+Z*z3;
      float textureX1 = vtx1*X; 
      float textureX2 = vtx2*X;
      float textureY1 = vty1*Z;
      float textureY2 = vty2*Z;
      float textureZ1 = vtz1*Y;
      float textureZ2 = vtz2*Y;
   float n[3];
      // gouraud 
   /*
      float vXn1=Xn1*X;
      float vYn1=Xn2*X;
      float vZn1=Xn3*X;
      float vXn2=Yn1*Z;
      float vYn2=Yn2*Z;
      float vZn2=Yn3*Z;
      float vXn3=Zn1*Y;
      float vYn3=Zn2*Y;
      float vZn3=Zn3*Y;
   
      n[0]=vXn1+vXn2+vXn3;
      n[1]=vYn1+vYn2+vYn3;
      n[2]=vZn1+vZn2+vZn3;
      normalize(n[0],n[1],n[2],&n[0],&n[1],&n[2]);
      produitScal(lumX,lumY,lumZ,n[0],n[1],n[2],&lum2);
   */
      // normal mapping
      TGAColor color2=normal.get((textureX1+textureY1+textureZ1),(textureX2+textureY2+textureZ2));
      n[0]=color2.r;
      n[1]=color2.g;
      n[2]=color2.b;
      normalize(n[0],n[1],n[2],&n[0],&n[1],&n[2]);
      produitScal(lumX,lumY,lumZ,n[0],n[1],n[2],&lum2);
      lum2=lum2*2-1;
      //  cout << z<< endl;
      float X_1,Y_1,Z_1;
      // cout << x << "/" << y << endl;
      
      if(X>=0 && Y>=0 && Z>=0){
	if(tabZ[tailleX*y+x]<z){
	   TGAColor color=text.get((textureX1+textureY1+textureZ1),(textureX2+textureY2+textureZ2));
	   
	   //  TGAColor color=TGAColor(255,255,255,255);
	  color.r=color.r*lum2;
	  
	  color.g=color.g*lum2;
	  
	  color.b=color.b*lum2;
	   
	  image.set(x,y,color);
	  tabZ[tailleX*y+x]=z;
	}
      }
    }
  }
}


void triangleSwip(int x1, int y1, int x2, int y2, int x3, int y3, TGAImage &image, TGAColor color){
  // cout << "on est rentrer ici" << endl;
  float a,b;
  int X1,X2;
  //tri des 3 points par rapport a leur y
  if(y1>y2){ std::swap(x1,x2);std::swap(y1,y2);}
  if(y1>y3){ std::swap(x1,x3);std::swap(y1,y3);}
  if(y2>y3){ std::swap(x2,x3);std::swap(y2,y3);}
  //premier triangle
  int hauteurmax=y3-y1;
  for(int y=y1;y<y2;y++){
    int hauteur=y2-y1+1;
    a=(float)(y-y1)/hauteurmax;
    /*if(hauteur!=0)*/b=(float)(y-y1)/hauteur;
    // else b=(float)(y-y1);
    X1=x1+(x3-x1)*a;
    X2=x1+(x2-x1)*b;
    // cout << X1 << endl;
    //if(X1>0){
    if(X1>X2){std::swap(X1,X2);}
      for(int i=X1;i<X2;i++){
	//cout << "on colorie ! " << endl;
	image.set(i,y,color);
      }
    //}
  }
  //deuxieme triangle
  for(int y=y2;y<y3;y++){
    int hauteur=y3-y2+1;
    a=(float)(y-y1)/hauteurmax;
    /*if(hauteur!=0)*/ b=(float)(y-y2)/hauteur;
    // else b=(float)(y-y2);
    X1=x1+(x3-x1)*a;
    X2=x2+(x3-x2)*b;
    //cout << X1 << endl;
    if(X1>0){
    if(X1>X2){std::swap(X1,X2);}
    for(int i=X1;i<X2;i++){
      image.set(i,y,color);
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


int main(int argc, char** argv) {
  int reso= 1000;
  TGAImage image(reso,reso, TGAImage::RGB);
  std::ifstream fichier("obj/african_head.obj",std::ios::in);
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
  float regardCamX=0;
  float regardCamY=0;
  float regardCamZ=0;
  float transfo[16];
  float transfo2[16];
  float transfo3[16];
  double rot[16];
  double rot2[16];
  double rot3[16];
  double rotTotal[16];
  float camera[16];
  float axe[16];
  double tmp[16];
  double alpha=0;//y
  double beta=0;//x
  double gama=0;//z
  alpha=rad(alpha);
  beta=rad(beta);
  gama=rad(gama);
  float res[16];
  cout << "debut des matrice" << endl;
  // rotation par rapport a y 
  rot[0]=cos(alpha);rot[1]=0;rot[2]=-sin(alpha);rot[3]=0;rot[4]=0;rot[5]=1;rot[6]=0;rot[7]=0;rot[8]=sin(alpha);rot[9]=0;rot[10]=cos(alpha);rot[11]=0;rot[12]=0;rot[13]=0;rot[14]=0;rot[15]=1;

  // roation par rapport a x 

  rot2[0]=1;rot2[1]=0;rot2[2]=0;rot2[3]=0;rot2[4]=0;rot2[5]=cos(beta);rot2[6]=-sin(beta);rot2[7]=0;rot2[8]=0;rot2[9]=sin(beta);rot2[10]=cos(beta);rot2[11]=0;rot2[12]=0;rot2[13]=0;rot2[14]=0;rot2[15]=1;

 // rotation par rapport a z

 rot3[0]=cos(gama);rot3[1]=-sin(gama);rot3[2]=0;rot3[3]=0;rot3[4]=sin(gama);rot3[5]=cos(gama);rot3[6]=0;rot3[7]=0;rot3[8]=0;rot3[9]=0;rot3[10]=1;rot3[11]=0;rot3[12]=0;rot3[13]=0;rot3[14]=0;rot3[15]=1;
 // matruce camera
 for(int i =0;i<16;i++)camera[i]=0;
 camera[3]=cameraX;
 camera[7]=cameraY;
 camera[11]=cameraZ;
 // matrive viewport
  transfo[0]=reso/2;transfo[1]=0;transfo[2]=0;transfo[3]=reso/2;transfo[4]=0;transfo[5]=reso/2;transfo[6]=0;transfo[7]=reso/2;transfo[8]=0;transfo[9]=0;transfo[10]=reso/8;transfo[11]=reso/8;transfo[12]=0;transfo[13]=0;transfo[14]=0;transfo[15]=1;
  //matrice rot total
  multMatSimpl4(rot,rot2,tmp);
  //printmat44(tmp);
 
  multMatSimpl4(tmp,rot3,rotTotal);
 
 
  //matrice axer camera 
  addMat44(rotTotal,camera,axe);
  cout << "axe" << endl;
  printmat44(axe);
  // matrice perspective
  transfo2[0]=1;transfo2[1]=0;transfo2[2]=0;transfo2[3]=0;transfo2[4]=0;transfo2[5]=1;transfo2[6]=0;transfo2[7]=0;transfo2[8]=0;transfo2[9]=0;transfo2[10]=1;transfo2[11]=0;transfo2[12]=0;transfo2[13]=0;transfo2[14]=-1/cameraZ;transfo2[15]=1;

 cout << "fin debut des matrice" << endl;

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
 multMatSimpl4(axe,transfo2,resi);
resi[15]=1;
 printmat44(resi);

 cout << endl;
 multMatSimpl4(transfo,resi,res);
 printmat44(res);

 std::cout << std::endl;

 //return 0;

 cout << "fin mult matrice" << endl;
  if(fichier){
    printf("Lecture du fichier : \n");
    std::string ligne;
    string p1,p2,p3;
    int x[2942],y[2942],z[2942];
    int vx[2942],vy[2942],vz[2942];
    int vnx[2942],vny[2942],vnz[2942];
    size_t str,str2,str3;
    int i=0;
    int j=0;
    int k=0;
    int l=0;
    int nbpoint =1258;
    float x1[nbpoint],y1[nbpoint],z1[nbpoint];
    float vtx[1339],vty[1339],vtz[1339];
    float vnx1[nbpoint],vny1[nbpoint],vnz1[nbpoint];
    string id0;
    int tailletotal=image.get_width()*image.get_height();
    float tabZ[tailletotal];
    for(int i=0;i<tailletotal;i++)tabZ[i]=-5000;
    while( fichier.good()){
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
	  /*
	  X_1=x1[x[i]-1]*reso/2+reso/2;
	  X_2=x1[y[i]-1]*reso/2+reso/2;
	  X_3=x1[z[i]-1]*reso/2+reso/2;
	  Y_1=y1[x[i]-1]*reso/2+reso/2;
	  Y_2=y1[y[i]-1]*reso/2+reso/2;
	  Y_3=y1[z[i]-1]*reso/2+reso/2;
	  Z_1=z1[x[i]-1]*reso/2+reso/2;
	  Z_2=z1[y[i]-1]*reso/2+reso/2;
	  Z_3=z1[z[i]-1]*reso/2+reso/2;
	  */
	 
	  // transformationFinal(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,reso/2,0,0,reso/2,0,reso/2,0,reso/2,0,0,reso/2,reso/2,0,0,0,1,&X1,&Y1,&Z1,&X2,&Y2,&Z2,&X3,&Y3,&Z3);

	  transformationVec(&X1,&Y1,&Z1,res,&X_1,&Y_1,&Z_1);
	  transformationVec(&X2,&Y2,&Z2,res,&X_2,&Y_2,&Z_2);
	  transformationVec(&X3,&Y3,&Z3,res,&X_3,&Y_3,&Z_3);
	  /*    float tmp[9];
	  tmp[0]=X_1;
	  tmp[1]=Y_1;
	  tmp[2]=Z_1;
	    printVec(tmp);
	  	  cout << endl;
	  tmp[0]=X_2;
	  tmp[1]=Y_2;
	  tmp[2]=Z_2;
	   printVec(tmp);
	  cout << endl;
	  tmp[0]=X_3;
	  tmp[1]=Y_3;
	  tmp[2]=Z_3;

	  printVec(tmp);

	  cout << endl;
	  */
	  X_1=max(X_1,0.f);
	  X_2=max(X_2,0.f);
	  X_3=max(X_3,0.f);
	  Y_1=max(Y_1,0.f);
	  Y_2=max(Y_2,0.f);
	  Y_3=max(Y_3,0.f);
	  Z_1=max(Z_1,0.f);
	  Z_2=max(Z_2,0.f);
	  Z_3=max(Z_3,0.f);
	 
	   /*   cout << X1 <<"/"<< X1*reso/2+reso/2<<"/"<< X_1 << endl;
	       /* cout << X2 <<"/"<< X2*reso/2+reso/2<<"/"<< X_2+reso/2 << endl;
	    cout << X3 <<"/"<< X3*reso/2+reso/2<<"/"<< X_3+reso/2 << endl;
	    cout << Y1 <<"/"<< Y1*reso/2+reso/2<<"/"<< Y_1+reso/2 << endl;
	    cout << Y2 <<"/"<< Y2*reso/2+reso/2<<"/"<< Y_2+reso/2 << endl;
	    cout << Y3 <<"/"<< Y3*reso/2+reso/2<<"/"<< Y_3+reso/2 << endl;
	    cout << Z1 <<"/"<< Z1*reso/2+reso/2<<"/"<< Z_1+reso/2 << endl;
	    cout << Z2 <<"/"<< Z2*reso/2+reso/2<<"/"<< Z_2+reso/2 << endl;
	    cout << Z3 <<"/"<< Z3*reso/2+reso/2<<"/"<< Z_3+reso/2 << endl;*/
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
	  Xn1=vnx1[vnx[i]-1]*resow;
	  Xn2=vny1[vnx[i]-1]*resoh;
	  Xn3=vnz1[vnx[i]-1];
	  Yn1=vnx1[vny[i]-1]*resow;
	  Yn2=vny1[vny[i]-1]*resoh;
	  Yn3=vnz1[vnz[i]-1];
	  Zn1=vnx1[vnz[i]-1]*resow;
	  Zn2=vny1[vnz[i]-1]*resoh;
	  Zn3=vnz1[vnz[i]-1];
	  //	line(X1,Y1,X2,Y2,image,white);
	 
	  //line(X2,Y2,X3,Y3,image,white);
	 
	  //line(X3,Y3,X1,Y1,image,white);
		
	  //	triangleSwip(X1,Y1,X2,Y2,X3,Y3,image,TGAColor(rand()%255, rand()%255,   rand()%255,   255));
	  triangleBarycentre(X_1,Y_1,Z_1,X_2,Y_2,Z_2,X_3,Y_3,Z_3,image,texture,normal,eclairage,VXu,VXv,VYu,VYv,VZu,VZv,tabZ,res,Xn1,Xn2,Xn3,Yn1,Yn2,Yn3,Zn1,Zn2,Zn3);
		  //triangleBarycentre(X_1,Y_1,Z_1,X_2,Y_2,Z_2,X_3,Y_3,Z_3,image,texture,VXu,VXv,VYu,VYv,VZu,VZv,tabZ);
	  i++;
	}else{
	  getline(fichier,ligne);
	}
      }
    }

    fichier.close();
    printf("On a fini de lire le fichier\n");
    /*	  
	  for(int i=0;i<nbpoint;i++){
	  image.set(x1[i]*reso/2+reso/2,y1[i]*reso/2+reso/2,red);
	  }
    */
    //cout << "on a fini de mettre les point " << endl;

  }else{
  }

  image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
  image.write_tga_file("output.tga");
  return 0;
}

