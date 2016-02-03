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
  

int appartientAuTriangle(int x, int y, float a1, float a2, float a3, float b1, float b2, float b3){
  int res=0;
  int machin1=a1*x+b1;
  int machin2=a2*x+b2;
  int machin3=a3*x+b3;
  if(y-machin1<0)
    if(y-machin2<0)
      if(y-machin3<0)res=1;
  return res;
}
void triangleBarycentre(int x1, int y1, int x2, int y2, int x3, int y3, TGAImage &image, TGAColor color){
  int xmin,xmax,ymin,ymax;
  float a1,a2,a3;
  float b1,b2,b3;
  //boxe 
  xmin=fmin(x1,x2);
  xmin=fmin(xmin,x3);
  ymin=fmin(y1,y2);
  ymin=fmin(ymin,y3);
  xmax=fmax(x1,x2);
  xmax=fmax(xmax,x3);
  ymax=fmax(y1,y2);
  ymax=fmax(ymax,y3);
  // a se debrouiller de facon que les coef directeur auront des vecteur direct positif vers l'interieur du triangle
  a1=(y2-y1)/(x2-x1);
  a2=(y3-y2)/(x3-x2);
  a3=(y1-y3)/(x1-x3);
  b1=y1-a1*x1;
  b2=y2-a2*x2;
  b3=y3-a3*x3;
  for(int i = xmin ; i<xmax ; i++){
    for(int j = ymin ; j<ymax ; j++){
      if(appartientAuTriangle(i,j,a1,a2,a3,b1,b2,b3)){
	image.set(i,j,white);
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
  int hauteurmax=y1-y2;
  for(int y=y1;y<y2;y++){
    int hauteur=y2-y1+1;
    a=(float)(y-y1)/hauteurmax;
    /*if(hauteur!=0)*/b=(float)(y-y1)/hauteur;
    // else b=(float)(y-y1);
    X1=x1+(x3-x1)*a;
    X2=x1+(x2-x1)*b;
    cout << X1 << endl;
    if(X1>0){
    if(X1>X2){std::swap(X1,X2);}
      for(int i=X1;i<X2;i++){
	//cout << "on colorie ! " << endl;
	image.set(i,y,color);
      }
    }
  }
  //deuxieme triangle
  for(int y=y2;y<y3;y++){
    int hauteur=y3-y2+1;
    a=(float)(y-y1)/hauteurmax;
    /*if(hauteur!=0)*/ b=(float)(y-y2)/hauteur;
    // else b=(float)(y-y2);
    X1=x1+(x3-x1)*a;
    X2=x2+(x3-x2)*b;
    cout << X1 << endl;
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
	if(fichier){
	  printf("Lecture du fichier : \n");
	  std::string ligne;
	  string p1,p2,p3;
	  int x[2942],y[2942],z[2942];
	  size_t str,str2,str3;
	  int i=0;
	  int j=0;
	  int nbpoint =1258;
	  float x1[nbpoint],y1[nbpoint],z1[nbpoint];
	  string id0;
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
		  cout << j << endl;
		  j++;
		}else{
		  getline(fichier,ligne);
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
		str=p2.find('/');
		y[i]=atoi((p2.substr(0,str)).c_str());
		str=p3.find('/');
		z[i]=atoi((p3.substr(0,str)).c_str());
		/*
		cout << x[i] << endl ;
		cout << y[i] << endl ;
		cout << z[i] << endl ;
		*/
		
		float X1,X2,X3;
		float Y1,Y2,Y3;
		X1=x1[x[i]-1]*reso/2+reso/2;
		X2=x1[y[i]-1]*reso/2+reso/2;
		X3=x1[z[i]-1]*reso/2+reso/2;
		Y1=y1[x[i]-1]*reso/2+reso/2;
		Y2=y1[y[i]-1]*reso/2+reso/2;
		Y3=y1[z[i]-1]*reso/2+reso/2;

		
		line(X1,Y1,X2,Y2,image,white);
	 
		line(X2,Y2,X3,Y3,image,white);
	 
		line(X3,Y3,X1,Y1,image,white);
		
		triangleSwip(X1,Y1,X2,Y2,X3,Y3,image,white);
		i++;
	      }else{
		getline(fichier,ligne);
	      }
	      }
	  }

	  fichier.close();
	  printf("On a fini de lire le fichier\n");
	  
	for(int i=0;i<nbpoint;i++){
	  image.set(x1[i]*reso/2+reso/2,y1[i]*reso/2+reso/2,red);
	}
	cout << "on a fini de mettre les point " << endl;

	}else{
	}

	image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
	image.write_tga_file("output.tga");
	return 0;
}

