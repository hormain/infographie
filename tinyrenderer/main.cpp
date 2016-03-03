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

void multMatriceQuatreQuatre(float m1_1,float m1_2,float m1_3, float m1_4,float m2_1,float m2_2,float m2_3, float m2_4,float m3_1,float m3_2,float m3_3, float m3_4,float m4_1,float m4_2,float m4_3, float m4_4,float n1_1,float n1_2,float n1_3, float n1_4,float n2_1,float n2_2,float n2_3, float n2_4,float n3_1,float n3_2,float n3_3, float n3_4,float n4_1,float n4_2,float n4_3, float n4_4,float *M1_1, float *M1_2,float *M1_3, float *M1_4,float *M2_1,float *M2_2,float *M2_3, float *M2_4,float *M3_1,float *M3_2,float *M3_3, float *M3_4,float *M4_1,float *M4_2,float *M4_3, float *M4_4){
  *M1_1=m1_1*n1_1+m1_2*n2_1+m1_3*n3_1+m1_4*n4_1;
 *M1_2=m1_1*n1_2+m1_2*n2_2+m1_3*n3_2+m1_4*n4_2;
 *M1_3=m1_1*n1_3+m1_2*n2_3+m1_3*n3_3+m1_4*n4_3;
 *M1_4=m1_1*n1_4+m1_2*n2_4+m1_3*n3_4+m1_4*n4_4;

 *M2_1=m2_1*n1_1+m2_2*n2_1+m2_3*n3_1+m2_4*n4_1;
 *M2_2=m2_1*n1_2+m2_2*n2_2+m2_3*n3_2+m2_4*n4_2;
 *M2_3=m2_1*n1_3+m2_2*n2_3+m2_3*n3_3+m2_4*n4_3;
 *M2_4=m2_1*n1_4+m2_2*n2_4+m2_3*n3_4+m2_4*n4_4;

 *M3_1=m3_1*n1_1+m3_2*n2_1+m3_3*n3_1+m3_4*n4_1;
 *M3_2=m3_1*n1_2+m3_2*n2_2+m3_3*n3_2+m3_4*n4_2;
 *M3_3=m3_1*n1_3+m3_2*n2_3+m3_3*n3_3+m3_4*n4_3;
 *M3_4=m3_1*n1_4+m3_2*n2_4+m3_3*n3_4+m3_4*n4_4;

 *M4_1=m4_1*n1_1+m4_2*n2_1+m4_3*n3_1+m4_4*n4_1;
 *M4_2=m4_1*n1_2+m4_2*n2_2+m4_3*n3_2+m4_4*n4_2;
 *M4_3=m4_1*n1_3+m4_2*n2_3+m4_3*n3_3+m4_4*n4_3;
 *M4_4=m4_1*n1_4+m4_2*n2_4+m4_3*n3_4+m4_4*n4_4;

}

void transformationFinal(float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3,float m1_1,float m1_2,float m1_3, float m1_4,float m2_1,float m2_2,float m2_3, float m2_4,float m3_1,float m3_2,float m3_3, float m3_4,float m4_1,float m4_2,float m4_3, float m4_4,float *X1, float *Y1, float *Z1, float *X2, float *Y2, float *Z2, float *X3, float *Y3, float *Z3){
  //calcule de laa matrice centrale 3D
  float x_1, y_1,z_1, x_2, y_2, z_2, x_3, y_3, z_3;
  float caj=0;//coordonne que l'on rapjoute a chaque vecteur pour avoir la 4D
  x_1=x1*m1_1+y1*m1_2+z1*m1_3+caj*m1_4;
  y_1=x1*m2_1+y1*m2_2+z1*m2_3+caj*m2_4;
  z_1=x1*m3_1+y1*m3_2+z1*m3_3+caj*m3_4;
  x_2=x2*m1_1+y2*m1_2+z2*m1_3+caj*m1_4;
  y_2=x2*m2_1+y2*m2_2+z2*m2_3+caj*m2_4;
  z_2=x2*m3_1+y2*m3_2+z2*m3_3+caj*m3_4;
  x_3=x3*m1_1+y3*m1_2+z3*m1_3+caj*m1_4;
  y_3=x3*m2_1+y3*m2_2+z3*m2_3+caj*m2_4;
  z_3=x3*m3_1+y3*m3_2+z3*m3_3+caj*m3_4;
  //projection 4D -> 3D
  float rapport1, rapport2, rapport3;
  rapport1=x1*m4_1+y1*m4_2+z1*m4_3+caj*m4_4;
  rapport2=x2*m4_1+y2*m4_2+z2*m4_3+caj*m4_4;
  rapport3=x3*m4_1+y3*m4_2+z3*m4_3+caj*m4_4;
  *X1=x_1/rapport1;
  *Y1=y_1/rapport1;
  *Z1=z_1/rapport1;
  *X2=x_2/rapport2;
  *Y2=y_2/rapport2;
  *Z2=z_2/rapport2;
  *X3=x_3/rapport3;
  *Y3=y_3/rapport3;
  *Z3=z_3/rapport3;
}
void transformationVec(float *x1, float *y1, float *z1,float m1_1,float m1_2,float m1_3, float m1_4,float m2_1,float m2_2,float m2_3, float m2_4,float m3_1,float m3_2,float m3_3, float m3_4,float m4_1,float m4_2,float m4_3, float m4_4,float *X1, float *Y1, float *Z1){
  //calcule de la matrice centrale 3D
  float x_1,y_1,z_1;
  float caj=1.f;//coordonne que l'on rapjoute a chaque vecteur pour avoir la 4D
   cout << *x1 <<"/" << *y1 << endl;
  //cout << "transformation" << endl;
  x_1=*x1*m1_1+*y1*m1_2+*z1*m1_3+caj*m1_4;
  y_1=*x1*m2_1+*y1*m2_2+*z1*m2_3+caj*m2_4;
  z_1=*x1*m3_1+*y1*m3_2+*z1*m3_3+caj*m3_4;
  //cout << x_1 <<"/"<< y_1 <<"/"<< z_1 << endl;
  //projection 4D -> 3D
  float rapport1, rapport2, rapport3;
  //cout << m4_1 << "/" << m4_2 << "/" << m4_3 << "/" << m4_4 << endl;
  rapport1=*x1*m4_1+*y1*m4_2+*z1*m4_3+caj*m4_4;
  //  cout << rapport1 << endl;
  *X1=x_1/rapport1;
  *Y1=y_1/rapport1;
  *Z1=z_1/rapport1;
  // cout << *X1 << "/" << *Y1 <<"/" << *Z1 << endl;
}
void triangleBarycentre(float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3, TGAImage &image, TGAImage &text,float vtx1,float vtx2,float vty1,float vty2,float vtz1,float vtz2,float* tabZ, float* transfo){
  float xmin=0,xmax=0,ymin=0,ymax=0;
  float undeuxX=0,untroisX=0,undeuxY=0,untroisY=0,undeuxZ=0,untroisZ=0,pointunX=0,pointunY=0,pointunZ=0;
  float coordX=0,coordY=0,coordZ=0;
  int tailleX=image.get_width();
  //float tabZ[tailleX*image.get_height()];
  float lum=0;
  float planX=0,planY=0,planZ=0;
  float lumX=0,lumY=0,lumZ=1;
  float z;
  float x_1,y_1,z_1,x_2,y_2,z_2,x_3,y_3,z_3;
  // TGAColor color=TGAColor(rand()%255, rand()%255,   rand()%255,   255);
  // cout << x1 <<endl;
  //cout << x2 << endl;
  //boxe 
  box(x1,y1,x2,y2,x3,y3,&xmin,&xmax,&ymin,&ymax);
  //perspective
  transformationVec(&x1,&y1,&z1,transfo[0],transfo[1],transfo[2],transfo[3],transfo[4],transfo[5],transfo[6],transfo[7],transfo[8],transfo[9],transfo[10],transfo[11],transfo[12],transfo[13],transfo[14],transfo[15],&x_1,&y_1,&z_1);
 transformationVec(&x2,&y2,&z2,transfo[0],transfo[1],transfo[2],transfo[3],transfo[4],transfo[5],transfo[6],transfo[7],transfo[8],transfo[9],transfo[10],transfo[11],transfo[12],transfo[13],transfo[14],transfo[15],&x_2,&y_2,&z_2);
 transformationVec(&x3,&y3,&z3,transfo[0],transfo[1],transfo[2],transfo[3],transfo[4],transfo[5],transfo[6],transfo[7],transfo[8],transfo[9],transfo[10],transfo[11],transfo[12],transfo[13],transfo[14],transfo[15],&x_3,&y_3,&z_3);
  //repere former de vers 1 vers 2 et 1 vers 3
  coordonneVec(x_1,y_1,z_1,x_2,y_2,z_2,&undeuxX,&undeuxY,&undeuxZ);
  coordonneVec(x_1,y_1,z_1,x_3,y_3,z_3,&untroisX,&untroisY,&untroisZ);
  produitVec(undeuxX,undeuxY,undeuxZ,untroisX,untroisY,untroisZ,&planX,&planY,&planZ);
  produitScal(planX,planY,planZ,lumX,lumY,lumZ,&lum);
  
   lum=abs(lum);
  // cout << xmin << '/' << xmax << '/' << ymin << '/' << ymax <<'/'<< lum << endl;
  for(int x=xmin;x<xmax;x++){
    for(int y=ymin;y<ymax;y++){
 float X=0,Y=0,Z=0;
       cout << x << "/" << y << endl;

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
      //  cout << z<< endl;
      float X_1,Y_1,Z_1;
      cout << x << "/" << y << endl;
      
      if(X>=0 && Y>=0 && Z>=0){
	if(tabZ[tailleX*y+x]<z){
	  TGAColor color=text.get((textureX1+textureY1+textureZ1),(textureX2+textureY2+textureZ2));
	  //TGAColor color=text.get(x,y);
	  // cout << coordX << '/' << coordY <<'/'<< coordZ << endl;
	  color.r=color.r*lum;
	  color.g=color.g*lum;
	  color.b=color.b*lum;	
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


int main(int argc, char** argv) {
  int reso= 1000;
  TGAImage image(reso,reso, TGAImage::RGB);
  std::ifstream fichier("obj/african_head.obj",std::ios::in);
  TGAImage texture;
  texture.read_tga_file("obj/african_head_diffuse.tga");
  texture.flip_vertically();
  int resow=texture.get_width();
  int resoh=texture.get_height();
  float camera=10;
  float transfo[16];
  multMatriceQuatreQuatre(reso/2,0,0,reso/2,0,reso/2,0,reso/2,0,0,reso/2,reso/2,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,0,0,/*-1./camera*/0,1,&transfo[0],&transfo[1],&transfo[2],&transfo[3],&transfo[4],&transfo[5],&transfo[6],&transfo[7],&transfo[8],&transfo[9],&transfo[10],&transfo[11],&transfo[12],&transfo[13],&transfo[14],&transfo[15]);
   /*multMatriceQuatreQuatre(reso/2,0,0,reso/2,0,reso/2,0,reso/2,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,&transfo[0],&transfo[1],&transfo[2],&transfo[3],&transfo[4],&transfo[5],&transfo[6],&transfo[7],&transfo[8],&transfo[9],&transfo[10],&transfo[11],&transfo[12],&transfo[13],&transfo[14],&transfo[15]);*/
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
	    //cout << ligne << endl;
	    /*
	      str=ligne.find(' ');
	      str2=ligne.find(' ',str+1);
	      str3=ligne.find(' ',str2+1);
	      x1[j]=ligne.substr(str,str2-str);
	      y1[j]=ligne.substr(str2,str3-str2);
	      z1[j]=ligne.substr(str3);
	    */
	    fichier >> x1[j];
	    fichier >> y1[j];
	    fichier >> z1[j];	 
	    cout << x1[j] << endl;
	    cout << y1[j] << endl;
	    cout << z1[j] << endl;
	    //cout << j << endl;
	    j++;
	  }else{
	    if(id0.at(1)=='t'){
	      fichier >> vtx[k];
	      fichier >> vty[k];
	      fichier >> vtz[k];
	      //cout << vtx[k] << endl;
	      //cout << vty[k] << endl;
	      //cout << vtz[k] << endl;
	      k++;
	    }else{
	      if(id0.at(1)=='n'){
		fichier >> vnx1[l];
		fichier >> vny1[l];
		fichier >> vnz1[l];
		//cout << vnx1[l] << endl;
		//cout << vny1[l] << endl;
		//cout << vnz1[l] << endl;
		l++;		     
	      }else{
		getline(fichier,ligne);
	      }
	    }
	  }
	}
	if(id==f){
	  getline(fichier, ligne);
	  //char trash;
	  //int gore;
	  /*
	    ligne >> x[i];
	    ligne >> trash;
	    atoi(ligne.c_str()) >> gore;
	    ligne >> trash;
	    atoi(ligne.c) >> gore ;
	    ligne >> y[i];
	    ligne >> trash;
	    atoi(ligne) >> gore;
	    ligne >> trash;
	    atoi(ligne) >> gore ;
	    ligne >> z[i];
	    ligne >> trash;
	    atoi(ligne) >> gore;
	    ligne >> trash;
	    atoi(ligne) >> gore ;
	  */
	  //	ligne >> y[i] >> z 
	  //	cout << ligne << endl;
	  str=ligne.find(' ');
	  //cout << str << endl;
	  str2=ligne.find(' ', str+1);
	  //cout << str2 << endl ;
	  p1=(ligne.substr(str,str2));
	  //cout << x << endl ;
	  str3=ligne.find(' ', str2+1);
	  p2=(ligne.substr(str2,str3));
	  //cout << y << endl ;
	  p3=(ligne.substr(str3,ligne.length()));
	  //cout << z << endl;
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
	  //  cout << vnx[i] <<'/'<< vny[i] <<'/'<< vnz[i] << endl;
	  /*
	    cout << x[i] << endl ;
	    cout << y[i] << endl ;
	    cout << z[i] << endl ;
	  */		
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
	  
	  /*
	  X1=x1[x[i]-1]*reso/2+reso/2;
	  X2=x1[y[i]-1]*reso/2+reso/2;
	  X3=x1[z[i]-1]*reso/2+reso/2;
	  Y1=y1[x[i]-1]*reso/2+reso/2;
	  Y2=y1[y[i]-1]*reso/2+reso/2;
	  Y3=y1[z[i]-1]*reso/2+reso/2;
	  Z1=z1[x[i]-1]*reso/2+reso/2;
	  Z2=z1[y[i]-1]*reso/2+reso/2;
	  Z3=z1[z[i]-1]*reso/2+reso/2;
	  */
	  float X_1=0,X_2=0,X_3=0,Y_1=0,Y_2=0,Y_3=0,Z_1=0,Z_2=0,Z_3=0;
	  // transformationFinal(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,reso/2,0,0,reso/2,0,reso/2,0,reso/2,0,0,reso/2,reso/2,0,0,0,1,&X1,&Y1,&Z1,&X2,&Y2,&Z2,&X3,&Y3,&Z3);

	  /* transformationVec(X1,Y1,Z1,transfo[0],transfo[1],transfo[2],transfo[3],transfo[4],transfo[5],transfo[6],transfo[7],transfo[8],transfo[9],transfo[10],transfo[11],transfo[12],transfo[13],transfo[14],transfo[15],&X_1,&Y_1,&Z_1);
	  transformationVec(X2,Y2,Z2,transfo[0],transfo[1],transfo[2],transfo[3],transfo[4],transfo[5],transfo[6],transfo[7],transfo[8],transfo[9],transfo[10],transfo[11],transfo[12],transfo[13],transfo[14],transfo[15],&X_2,&Y_2,&Z_2);
	  transformationVec(X3,Y3,Z3,transfo[0],transfo[1],transfo[2],transfo[3],transfo[4],transfo[5],transfo[6],transfo[7],transfo[8],transfo[9],transfo[10],transfo[11],transfo[12],transfo[13],transfo[14],transfo[15],&X_3,&Y_3,&Z_3);*/
	  //  cout << transfo[0] <<"/" <<  transfo[1] <<"/" <<  transfo[2] << "/" << transfo[3] << "/" << transfo[4] << "/" << transfo[5] << "/" << transfo[6] << "/" << transfo[7] << "/" << transfo[8] << "/" << transfo[9] << "/" << transfo[10] << "/" << transfo[11] << "/" << transfo[12] << "/" << transfo[13] << "/" << transfo[14] << "/" << transfo[15] << endl;
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

	  //	line(X1,Y1,X2,Y2,image,white);
	 
	  //line(X2,Y2,X3,Y3,image,white);
	 
	  //line(X3,Y3,X1,Y1,image,white);
		
	  //	triangleSwip(X1,Y1,X2,Y2,X3,Y3,image,TGAColor(rand()%255, rand()%255,   rand()%255,   255));
	  triangleBarycentre(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,image,texture,VXu,VXv,VYu,VYv,VZu,VZv,tabZ,transfo);
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

