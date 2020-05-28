#include "math.h"
#include <iostream>
#include "string.h"
#include "MCCA_functions.H"







double MCCA_dot_product(double A[], double B[]) 
{
	return (A[0]*B[0]+A[1]*B[1]+A[2]*B[2]); 
}


void MCCA_polygon_cross_product(double a[],double b[], double c[])
{

	c[0] = a[1]*b[2] - b[1]*a[2];
    c[1] = a[2]*b[0] - b[2]*a[0];
    c[2] = a[0]*b[1] - b[0]*a[1];
}

void MCCA_get_polygon_normal(double v1[],double v2[],double v3[],double norm[])
{
	// produce the unit vector normal to the polygon
	// local variables
    double a[3],b[3];
    double mag;
    double tiny=1.0e-14;

    a[0] = v2[0] - v1[0];a[1] = v2[1] - v1[1];a[2] = v2[2] - v1[2];
    b[0] = v3[0] - v2[0];b[1] = v3[1] - v2[1];b[2] = v3[2] - v2[2];
    MCCA_polygon_cross_product(a,b,norm);
    mag  = std::max(pow((MCCA_dot_product(norm,norm)),0.5), tiny);
    norm[0] = norm[0]/mag;
    norm[1] = norm[1]/mag;
    norm[2] = norm[2]/mag;

}


void MCCA_get_polygon_local_coord(double x[],double y[],double z[],int &n,double xprime[],double yprime[],double zprime[])
{
	//implicit none
	//convert a polygon in (x,y,z) to a polygon in (x',y',0)
	//input
	//x(n)    = x-coordinates of polygon vertices
	//y(n)    = y-coordinates of polygon vertices
	//z(n)    = z-coordinates of polygon vertices
	//n       = number of vertices

	//output
	//xprime(n) = x-coordinates of polygon vertices
	//yprime(n) = y-coordinates of polygon vertices
	//zprime(n) = z-coordinates of polygon vertices

	//local variables
    int i;
    double xc[n],yc[n],zc[n];
    double iprime[3], jprime[3], kprime[3], pv[3], mag;

    
    for(int i=0;i<n;i++)
    {
    	xc[i]=x[i]; yc[i]=y[i]; zc[i]=z[i];
       	//subtract the displacement from all points
        xc[i]=xc[i]-x[0];
        yc[i]=yc[i]-y[0];
        zc[i]=zc[i]-z[0];
    }




    //select point 1 to point 2 as the x direction
    //unit vector in the xprime direction
    iprime[0]=xc[1]-xc[0]; iprime[1]=yc[1]-yc[0]; iprime[2]=zc[1]-zc[0];
    mag = std::max(pow(MCCA_dot_product(iprime,iprime),0.5), 1.0e-14);
    iprime[0]=iprime[0]/mag; iprime[1]=iprime[1]/mag; iprime[2]=iprime[2]/mag;
    
    //vector from point 2 to point 3
    pv[0]=xc[2]-xc[1]; pv[1]=yc[2]-yc[1]; pv[2]=zc[2]-zc[1];

    //unit vector in the zprime direction
    MCCA_polygon_cross_product(iprime,pv,kprime);
    mag = std::max(pow(MCCA_dot_product(kprime,kprime),0.5), 1.0e-14);
    kprime[0] = kprime[0]/mag; kprime[1] = kprime[1]/mag; kprime[2] = kprime[2]/mag;


    //unit vector jprime in the yprime direction
    MCCA_polygon_cross_product(kprime,iprime,jprime);
    mag = std::max(pow(MCCA_dot_product(jprime,jprime),0.5), 1.0e-14);
    jprime[0] = jprime[0]/mag; jprime[1]=jprime[1]/mag; jprime[2]=jprime[2]/mag;

    //for each point find the projections of xprime, yprime, and zprime
    //all zprime values should be zero
    for(int i=0; i<n; i++)
    {
    	pv[0]=xc[i]; pv[1]=yc[i]; pv[2]=zc[i];
        xprime[i] = MCCA_dot_product(iprime,pv);
        yprime[i] = MCCA_dot_product(jprime,pv);
        zprime[i] = MCCA_dot_product(kprime,pv);
    }

}



double MCCA_irregular_polygon_area(double x[],double y[],double z[],int &n)
{
	//given the coordinates of the vertices, this routine applies 
	//the surveyors formula for the area and centroid of an irregular 2d polygon
	//the vertices may be ordered clockwise or counterclockwise. if they are 
	//ordered clockwise, the area will be negative but correct in absolute value.

	//it is assumed the vertices are ordered in a counterclockwise manner.
	//x(n)    = x-coordinates of vertices
	//y(n)    = y-coordinates of vertices
	//z(n)    = z-coordinates of vertices
	//n       = number of vertices
	//area    = area of irregular polygon


	//local variables
    double xp[n],yp[n],zp[n],det,sum1,sum2,sum3;
    double xcenter, ycenter;
    double sixth = 1.0/6.0;
    double area;
    double tiny=1.0e-15;

	//put the polygon in its own x-y coordinates; all zp's should be zero

    MCCA_get_polygon_local_coord(x,y,z,n,xp,yp,zp);

	//apply the surveyor's formula
    sum1 = 0.0;
    sum2 = 0.0;
    sum3 = 0.0;
    for(int i=0; i<n-1; i++)
    {
    	det  = xp[i]*yp[i+1] - xp[i+1]*yp[i];
        sum1 = sum1 + det;
        sum2 = sum2 + det*(xp[i] + xp[i+1]);
        sum3 = sum3 + det*(yp[i] + yp[i+1]);
    }


    //close the polygon
    det  = xp[n-1]*yp[0] - xp[0]*yp[n-1];
    sum1 = sum1 + det;
    sum2 = sum2 + det*(xp[n-1] + xp[0]);
    sum3 = sum3 + det*(yp[n-1] + yp[0]);



    //area and center
    area    = 0.5*sum1;
    det     = 1.0/(area+tiny);
    xcenter = sixth*det*sum2;
    ycenter = sixth*det*sum3;
    area    = abs(area);

    return area;
}



double MCCA_pyramid_height(double x[],double y[],double z[],int &n,double &px,double &py,double &pz)
{
	//get the height of the pyramid from the base to the top point
	//input
	//x(n)    = x-coordinates of polygon vertices
	//y(n)    = y-coordinates of polygon vertices
	//z(n)    = z-coordinates of polygon vertices
	//n       = number of vertices
	//px      = x-coordinate of point
	//py      = x-coordinate of point
	//px      = x-coordinate of point
	//output
	//dist    = distance from base to point
	//declare the pass
    //local variables
    double v1[3],v2[3],v3[3],norm[3],vt[3];
    double dist;


    //get unit vector normal to polygon 
    v1[0]=x[0]; v1[1]=y[0]; v1[2]=z[0];
    v2[0]=x[1]; v2[1]=y[1]; v2[2]=z[1];
    v3[0]=x[2]; v3[1]=y[2]; v3[2]=z[2];
    MCCA_get_polygon_normal(v1,v2,v3,norm);

    //the height is the projection of vt in the normal direction.
    //absolute value to be independent of orientation 
    vt[0] = px-x[1]; vt[1]=py-y[1]; vt[2]=pz-z[1];
    dist = abs(MCCA_dot_product(norm,vt));
    return dist;
}


double MCCA_pyramid_volume(double x[],double y[],double z[],int &n,double &px,double &py,double &pz)
{
	//volume of a pyramid
	//input
	//x(n)    = x-coordinates of polygon vertices
	//y(n)    = y-coordinates of polygon vertices
    //z(n)    = z-coordinates of polygon vertices
    //n       = number of vertices
    //px      = x-coordinate of point
    //py      = x-coordinate of point
    //px      = x-coordinate of point
    //output
    //volume = volume of pyramid
    //declare the pass
    //local variables
    double area,height,volume;
    double third = 1.0/3.0;

    area = MCCA_irregular_polygon_area(x,y,z,n);
    height = MCCA_pyramid_height(x,y,z,n,px,py,pz);
    volume = area * height * third;

    return volume;
}


void MCCA_add_vertex_to_face(int iface,int &nfaces,int (&iface_mask)[7],int (&face)[7],int (&nface_vertices)[7],int nvert,int (&ivertices)[7][7])
{
//adds a face for the plane - cube intersection

//go
      if (iface_mask[iface-1] == 0){
       iface_mask[iface-1] = 1;
       nfaces = nfaces + 1;
       face[nfaces-1] = iface;
      }
      nface_vertices[iface-1] = nface_vertices[iface-1] + 1;
      ivertices[nface_vertices[iface-1]-1][iface-1] = nvert;
      
}



void MCCA_check_edge(int iedge,int jedge,const double (&vx)[8],const double (&vy)[8],const double (&vz)[8], double a,double b,double c,double d,double &lambda,double &xx,double &yy,double &zz)
{
 

// computes the distance lambda along an edge between two vertices where an intersection occurs.
// lambda will be between 0 and 1 in the unit cube for a valid intersection

// input: 
// iedge = edge number
// jedge = edge numver
// vx, vy, vz = coordinates of unit cube
// a, b, c, d = coefficients of the hessian normal plane ax + b*y * c*z + d = 0

// output:
// lambda = distance along edges for an intersection
// xx, yy, zz = coordinates of intersection point

// local variables
double eij[3];
const double tiny = 1.0/pow(10,12);
const double supertiny = 1.0/pow(10,300);

eij[0] = vx[jedge] - vx[iedge] ; eij[1] = vy[jedge] - vy[iedge] ;  eij[2] = vz[jedge] - vz[iedge];

lambda = (d - (a*vx[iedge] + b*vy[iedge] + c*vz[iedge]) ) / (a*eij[0] + b*eij[1] + c*eij[2]+ supertiny);

if (abs(lambda) < tiny){ lambda = 0.0;}
if (abs(lambda - 1.0) < tiny){ lambda = 1.0;}


// coordinates of intersection 
xx = 0.0 ; yy = 0.0; zz = 0.0;

if (lambda >= 0.0 && lambda <= 1.0){
xx=vx[iedge]+lambda*eij[0];
yy=vy[iedge]+lambda*eij[1];
zz=vz[iedge]+lambda*eij[2];
}
}



void MCCA_check_if_vertex_exists(double (&x_intersect)[6],double (&y_intersect)[6],double (&z_intersect)[6],int &n,double &xx,double &yy,double &zz,int &idum)
{
//checks if a vertex already exists
//declare the pass
//local variables
const double tiny = 1.0/pow(10,12);

// loop over all existing intersection points
idum = 0;
for(int i =0; i<n ; i++){
	if (abs(x_intersect[i]-xx) <= tiny && abs(y_intersect[i]-yy) <= tiny && abs(z_intersect[i]-zz) <= tiny)
          {idum = 1; break;}
}

}



void MCCA_volume_fraction(double dl, double theta, double phi, double &volume, int &npoints, double (&x_final)[6], double (&y_final)[6], double (&z_final)[6])
{

// computes the volume of plane cube interface
//
// task 1 - find the plane-cube intersection points.
//          if this is all that was wanted, this routine would be a *lot* shorter.
//
// task 2 - form the list of face vertices.
//          this is the bulk of the routine. 
//
// task 3 - form the volume associated with each face.
//          relatively trivial once the face vertices are known
//
// it is useful to know the orderings in cube_path.pdf to understand this routine

// input:
// dl    = length of ray i.e., the radius from the origin
// theta = spherical coordinate angle theta of ray dl;  theta from 0 to 2*pi
// phi   = spherical coordinate angle phi of ray dl; phi from 0 to pi.

// output
// volume = volume of box cut by the intersection of the ray's normal plane and the cube
// npoints = number of intersection points; only values of 3, 4, 5, and 6 are possible.

// for convenience only
const double pi=3.1415926535897932384, a2rad=pi/180.0, rad2a = 180.0/pi;
int dummyinput;

// declare the pass
//int npoints;
//double dl, theta, phi, volume;

//dl = 0.5;
//theta = pi/4.0;
//phi = pi/4.0;

// ipath stores on which pathway an intersection was found
// ipath = 1 = light gray +x +z +y path
// ipath = 2 = gray +y +x +z path 
// ipath = 3 = black +z +y +x path
// ipath = 4 = dotted light gray, parallel to y-axis, path
// ipath = 5 = dotted gray, parallel to z-axis, path
// ipath = 6 = dotted black, parallel to x-axis, path

int ipath[12];

// for the plane-cube intersection points
// npoints is the number of intersections with the unit cube
// x_intersect, y_intersect, z_intersect store the coordinates of the intersection points 

double x_intersect[6], y_intersect[6], z_intersect[6];

// for the vertices
// nvert is the number of vertices, which includes both cube vertices (8) and intersection points (6 at maximum)
// x_vertex, y_vertex, z_vertex store the coordinates of the vertices

int nvert;
double x_vertex[14], y_vertex[14], z_vertex[14];

// for the faces
// nfaces counts the number of faces with intersection, a maximum of 7 - six for the cube plus one for the plane
// face stores the face number
// nface_vertices stores how many vertices are with each face
// iverticies stores the vertices of eah face
// iface_mask is a convenient was of turning faces on or off

int nfaces, face[7], nface_vertices[7], ivertices[7][7], iface_mask[7];

// for the edges
// kedge stores the edge number

int iedge, jedge, kedge[6];

// other variables
int m, n, idum, jdum, kdum;
double xn, yn, zn, mag, a, b, c, d, xx, yy, zz, lambda, xdum[8], ydum[8], zdum[8], area, dv, dl_max;

//coordinates of the unit cube
const double vx[8] = {0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0};
const double vy[8] = {0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0};
const double vz[8] = {0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0};



//initialize
npoints = 0;
volume  = -1.0; 

//check if there is an intersection at all
dl_max  = (cos(theta) + sin(theta))*sin(phi) + cos(phi);
if(dl > dl_max){ 
    x_final[0] = x_intersect[0]; x_final[1] = x_intersect[1]; x_final[2] = x_intersect[2]; x_final[3] = x_intersect[3]; x_final[4] = x_intersect[4]; x_final[5] = x_intersect[5];
    y_final[0] = y_intersect[0]; y_final[1] = y_intersect[1]; y_final[2] = y_intersect[2]; y_final[3] = y_intersect[3]; y_final[4] = y_intersect[4]; y_final[5] = y_intersect[5];
    z_final[0] = z_intersect[0]; z_final[1] = z_intersect[1]; z_final[2] = z_intersect[2]; z_final[3] = z_intersect[3]; z_final[4] = z_intersect[4]; z_final[5] = z_intersect[5];
    return;}

//coordinates of the normal vector of length dl
xn =  dl*cos(theta)*sin(phi);
yn =  dl*sin(theta)*sin(phi);
zn =  dl*cos(phi);

//hessian normal form of the cutting plane, a*x + b*y + c*z = d
mag = sqrt(xn*xn + yn*yn + zn*zn);
a  = xn/mag ;  b = yn/mag ; c = zn/mag ; d = dl;

//initialize all face and vertex data
//ivertices = 0;
double doubleZero = 0.0;
int   intZero = 0;

npoints = 0; nfaces = 0;
memset(x_vertex, doubleZero, sizeof(x_vertex));
memset(y_vertex, doubleZero, sizeof(y_vertex));
memset(z_vertex, doubleZero, sizeof(z_vertex));
memset(ipath, intZero, sizeof(ipath));
memset(iface_mask, intZero, sizeof(iface_mask));
memset(face, intZero, sizeof(face));
memset(nface_vertices, intZero, sizeof(nface_vertices));
memset(ivertices, intZero, sizeof(ivertices));

//starting vertex v1 contributes to three faces - 1, 4, and 5
nvert = 1;
x_vertex[nvert] = vx[0]  ; y_vertex[nvert] = vy[0]  ; z_vertex[nvert] = vz[0]; 
ipath[nvert] = 0;


MCCA_add_vertex_to_face(1,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
MCCA_add_vertex_to_face(4,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
MCCA_add_vertex_to_face(5,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);

//face 7, the intersecting plane, is slightly special - don't add this vertex, but set it up
nfaces = nfaces + 1 ; face[nfaces-1] = 7 ; iface_mask[6] = 1;


//-------------------------------------------------------------------------------------------------------------------     

// starting from vertex 1, first independent path is +x +z +y (light gray solid lines in cube_paths.pdf)
// an intersection can only occur once along this path, hence the triple if block

// check along +x 

iedge = 0 ; jedge = 1;
MCCA_check_edge(iedge,jedge,vx,vy,vz,a,b,c,d,lambda,xx,yy,zz);

//there is a valid intersection 

if (lambda >= 0.0 && lambda <= 1.0)
{

	//see if this vertex already exists

	MCCA_check_if_vertex_exists(x_intersect,y_intersect,z_intersect,npoints,xx,yy,zz,idum);


	//if not, add this vertex to faces 1,5,7
	npoints = npoints + 1;       
	kedge[npoints-1] = 1;
	x_intersect[npoints-1]=xx;y_intersect[npoints-1]=yy;z_intersect[npoints-1] = zz;

	nvert = nvert + 1;
	x_vertex[nvert-1] = x_intersect[npoints-1];y_vertex[nvert-1]=y_intersect[npoints-1];z_vertex[nvert-1]=z_intersect[npoints-1];
	ipath[nvert-1] = 1;
	MCCA_add_vertex_to_face(1,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
	MCCA_add_vertex_to_face(5,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
	MCCA_add_vertex_to_face(7,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);





	//if needed, add end point vertices to faces 4 and 2
	if(lambda == 0)
	{
		MCCA_add_vertex_to_face(4,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
		}
		else if(lambda ==1)
		{
		MCCA_add_vertex_to_face(2,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
		}
}

// no intersection along +x, now check along +z
// vertex v2, contributes to faces 1, 2 and 5
// less commentary as the pattern emerges

else
{
	nvert = nvert + 1;
	x_vertex[nvert-1] = vx[jedge] ; y_vertex[nvert-1] = vy[jedge] ; z_vertex[nvert-1] = vz[jedge];
	ipath[nvert-1] = 1;

	MCCA_add_vertex_to_face(1,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
	MCCA_add_vertex_to_face(2,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
	MCCA_add_vertex_to_face(5,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);

	iedge=1; jedge=4;
	MCCA_check_edge(iedge,jedge,vx,vy,vz,a,b,c,d,lambda,xx,yy,zz);

	if (lambda >= 0.0 && lambda <= 1.0)
	{
		MCCA_check_if_vertex_exists(x_intersect,y_intersect,z_intersect,npoints,xx,yy,zz,idum);
		if (idum == 0){
			npoints = npoints + 1;
			kedge[npoints] = 2;
			x_intersect[npoints-1] = xx ; y_intersect[npoints-1] = yy ; z_intersect[npoints-1] = zz;

			nvert = nvert + 1;
            x_vertex[nvert-1] = x_intersect[npoints-1] ; y_vertex[nvert-1] = y_intersect[npoints-1] ; z_vertex[nvert-1] = z_intersect[npoints-1];
            ipath[nvert-1] = 1;
            MCCA_add_vertex_to_face(2,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
            MCCA_add_vertex_to_face(5,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
            MCCA_add_vertex_to_face(7,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
            if(lambda==0){
            	MCCA_add_vertex_to_face(1,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
            }else if(lambda==1){
            	MCCA_add_vertex_to_face(3,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
            }
        }
	}
	//no intersection along +z either, now check +y
	//vertex v5, contributes to faces 2, 3 and 5 
	else{

		nvert = nvert + 1;
		x_vertex[nvert-1] = vx[jedge] ; y_vertex[nvert-1] = vy[jedge] ; z_vertex[nvert-1] = vz[jedge];
		ipath[nvert-1] = 1;
		MCCA_add_vertex_to_face(2,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
		MCCA_add_vertex_to_face(3,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
		MCCA_add_vertex_to_face(5,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);

		iedge = 4; jedge = 7;
		MCCA_check_edge(iedge,jedge,vx,vy,vz,a,b,c,d,lambda,xx,yy,zz);
		if (lambda >= 0.0 && lambda <= 1.0){
			MCCA_check_if_vertex_exists(x_intersect,y_intersect,z_intersect,npoints,xx,yy,zz,idum);
			if (idum == 0){
				npoints = npoints + 1;
				kedge[npoints-1] = 3;
				x_intersect[npoints-1] = xx ; y_intersect[npoints-1] = yy ; z_intersect[npoints-1] = zz;

				nvert = nvert + 1;
				x_vertex[nvert-1] = x_intersect[npoints-1] ; y_vertex[nvert-1] = y_intersect[npoints-1] ; z_vertex[nvert-1] = z_intersect[npoints-1];
				ipath[nvert-1] = 1;
				MCCA_add_vertex_to_face(2,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
				MCCA_add_vertex_to_face(3,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
				MCCA_add_vertex_to_face(7,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);

				if (lambda == 0.0){
					MCCA_add_vertex_to_face(5,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
				} else if (lambda == 1.0) {
					MCCA_add_vertex_to_face(6,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
				}

			}


        //end of +y lambda between 0 and 1 if block
		}
    //end of +z +y lambda between 0 and 1 if block
	}
//end of +x lambda between 0 and 1 if block
}


// line intersection +x is parallel to y 
// light gray dashed line between v2 and v6, in this direction, in figure cube_paths.pdf

iedge = 1 ; jedge = 5;
MCCA_check_edge(iedge,jedge,vx,vy,vz,a,b,c,d,lambda,xx,yy,zz);

if (lambda >= 0.0 && lambda <= 1.0){
	MCCA_check_if_vertex_exists(x_intersect,y_intersect,z_intersect,npoints,xx,yy,zz,idum);
	// add this vertex to faces 1, 2, and 7
	if(idum==0){
		npoints = npoints + 1;       
        kedge[npoints-1] = 4;
        x_intersect[npoints-1] = xx ; y_intersect[npoints-1] = yy ; z_intersect[npoints-1] = zz;

        nvert = nvert + 1;
        x_vertex[nvert-1] = x_intersect[npoints-1]; y_vertex[nvert-1] = y_intersect[npoints-1]; z_vertex[nvert-1] = z_intersect[npoints-1];
        ipath[nvert-1] = 4;
        MCCA_add_vertex_to_face(1,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
        MCCA_add_vertex_to_face(2,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
        MCCA_add_vertex_to_face(7,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
        //if needed, add end point vertices to faces 5 and 6
        if(lambda==0.0){
        	MCCA_add_vertex_to_face(5,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
        }else if(lambda == 1.0){
        	MCCA_add_vertex_to_face(6,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
        }

	} 

}





//--------------------------------------------------------------------------------------
//starting from vertex 1, second independent path is +y +x +z (gray solid lines in cube_paths.pdf)
//an intersection can only occur once along this path, hence the triple if block
//check +y

iedge = 0; jedge = 2;
MCCA_check_edge(iedge,jedge,vx,vy,vz,a,b,c,d,lambda,xx,yy,zz);
if (lambda >= 0.0 && lambda <= 1.0){
	MCCA_check_if_vertex_exists(x_intersect,y_intersect,z_intersect,npoints,xx,yy,zz,idum);
	//add this vertex to faces 1, 4, and 7
	if(idum == 0){
		npoints = npoints + 1;       
        kedge[npoints-1] = 5;
        x_intersect[npoints-1] = xx; y_intersect[npoints-1] = yy; z_intersect[npoints-1] = zz;

        nvert = nvert + 1;
        x_vertex[nvert-1] = x_intersect[npoints-1] ; y_vertex[nvert-1] = y_intersect[npoints-1] ; z_vertex[nvert-1] = z_intersect[npoints-1];
        ipath[nvert-1] = 2;
        MCCA_add_vertex_to_face(1,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
        MCCA_add_vertex_to_face(4,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
        MCCA_add_vertex_to_face(7,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
        //if needed, add end point vertices to faces 5 and 6
        if(lambda == 0.0){
        	MCCA_add_vertex_to_face(5,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
        }else if(lambda == 1.0){
        	MCCA_add_vertex_to_face(6,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
        }

	}
} 
//no intersection along +y, now check +x
//vertex v3 which contributes to faces 1, 4, and 6
else {
	nvert = nvert + 1; 
    x_vertex[nvert-1] = vx[jedge] ; y_vertex[nvert-1] = vy[jedge] ; z_vertex[nvert-1] = vz[jedge];
    ipath[nvert-1] = 2;
    MCCA_add_vertex_to_face(1,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
    MCCA_add_vertex_to_face(4,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
    MCCA_add_vertex_to_face(6,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);


    iedge = 2; jedge = 5;
    MCCA_check_edge(iedge,jedge,vx,vy,vz,a,b,c,d,lambda,xx,yy,zz);

    if (lambda >= 0.0 && lambda <= 1.0){
    	MCCA_check_if_vertex_exists(x_intersect,y_intersect,z_intersect,npoints,xx,yy,zz,idum);

    	//add this vertex to faces 1, 6, and 7
    	if(idum == 0){
    		npoints = npoints + 1;
    		kedge[npoints-1] = 6;
            x_intersect[npoints-1] = xx ; y_intersect[npoints-1] = yy ; z_intersect[npoints-1] = zz;

            nvert = nvert + 1;
            x_vertex[nvert-1] = x_intersect[npoints-1] ; y_vertex[nvert-1] = y_intersect[npoints-1] ; z_vertex[nvert-1] = z_intersect[npoints-1];
            ipath[nvert-1] = 2;
            MCCA_add_vertex_to_face(1,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
            MCCA_add_vertex_to_face(6,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
            MCCA_add_vertex_to_face(7,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
            //if needed, add end point vertices to faces 4 and 2
            if(lambda == 0.0){
            	MCCA_add_vertex_to_face(4,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
            }else if(lambda == 1.0){
            	MCCA_add_vertex_to_face(2,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
            }
    	}

    } 
    //no intersection along +x, now check +z
    //meaning we pick up vertex v6 which contributes to faces 1, 2, and  6
    else {
    	nvert = nvert + 1;
    	x_vertex[nvert-1] = vx[jedge] ; y_vertex[nvert-1] = vy[jedge] ; z_vertex[nvert-1] = vz[jedge];
        ipath[nvert-1] = 2;
        MCCA_add_vertex_to_face(1,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
        MCCA_add_vertex_to_face(2,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
        MCCA_add_vertex_to_face(6,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);

        iedge = 5 ; jedge = 7;
        MCCA_check_edge(iedge,jedge,vx,vy,vz,a,b,c,d,lambda,xx,yy,zz);

        if (lambda >= 0.0 && lambda <= 1.0){
        	MCCA_check_if_vertex_exists(x_intersect,y_intersect,z_intersect,npoints,xx,yy,zz,idum);
        	//add this vertex to faces 2, 6, and 7
        	if(idum == 0){
        		npoints = npoints + 1;       
                kedge[npoints-1] = 7;
                x_intersect[npoints-1] = xx ; y_intersect[npoints-1] = yy ; z_intersect[npoints-1] = zz;
                nvert = nvert + 1;
                x_vertex[nvert-1] = x_intersect[npoints-1] ; y_vertex[nvert-1] = y_intersect[npoints-1] ; z_vertex[nvert-1] = z_intersect[npoints-1];
                ipath[nvert-1] = 2;
                MCCA_add_vertex_to_face(2,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
                MCCA_add_vertex_to_face(6,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
                MCCA_add_vertex_to_face(7,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);

                //if needed, add end point vertices to faces 1 and 3
                if(lambda == 0.0){
                	MCCA_add_vertex_to_face(1,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
                }else if(lambda == 1.0){
                	MCCA_add_vertex_to_face(3,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
                }


        	}

        //end of +z lambda between 0 and 1 if block
        }

    //end of +x +z lambda between 0 and 1 if block
    }


//end of +y +x +z lambda between 0 and 1 if block
}

//line intersections +y parallel to +z
//gray dashed line between v3 and v7 in figure cube_paths.pdf
iedge = 2; jedge = 6;
MCCA_check_edge(iedge,jedge,vx,vy,vz,a,b,c,d,lambda,xx,yy,zz);
if (lambda >= 0.0 && lambda <= 1.0){
	MCCA_check_if_vertex_exists(x_intersect,y_intersect,z_intersect,npoints,xx,yy,zz,idum);

	//add this vertex to faces 2, 6, and 7
	if(idum == 0){
		npoints = npoints + 1;       
        kedge[npoints-1] = 8;
        x_intersect[npoints-1] = xx ; y_intersect[npoints-1] = yy ; z_intersect[npoints-1] = zz;

        nvert = nvert + 1;
        x_vertex[nvert-1] = x_intersect[npoints-1] ; y_vertex[nvert-1] = y_intersect[npoints-1] ; z_vertex[nvert-1] = z_intersect[npoints-1];
        ipath[nvert-1] = 5;
        MCCA_add_vertex_to_face(4,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
        MCCA_add_vertex_to_face(6,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
        MCCA_add_vertex_to_face(7,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);

        //if needed, add end point vertices to faces 1 and 3
        if(lambda == 0.0){
        	MCCA_add_vertex_to_face(1,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
        } else if(lambda == 1.0){
        	MCCA_add_vertex_to_face(3,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
        }
	}
}

//--------------------------------------------------------------------------------------
//starting from vertex 1, third independent path is +z +y +x (black solid lines in cube_paths.pdf)
//an intersection can only occur once along this path, hence the triple if block
//check +z

iedge = 0; jedge = 3;
MCCA_check_edge(iedge,jedge,vx,vy,vz,a,b,c,d,lambda,xx,yy,zz);

if (lambda >= 0.0 && lambda <= 1.0){
	MCCA_check_if_vertex_exists(x_intersect,y_intersect,z_intersect,npoints,xx,yy,zz,idum);

	//add this vertex to faces 4, 5, and 7
	if(idum == 0){
		npoints = npoints + 1;       
        kedge[npoints-1] = 9;
        x_intersect[npoints-1] = xx ; y_intersect[npoints-1] = yy ; z_intersect[npoints-1] = zz;

        nvert = nvert + 1;
        x_vertex[nvert-1] = x_intersect[npoints-1] ; y_vertex[nvert-1] = y_intersect[npoints-1] ; z_vertex[nvert-1] = z_intersect[npoints-1];
        ipath[nvert-1] = 3;
        MCCA_add_vertex_to_face(4,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
        MCCA_add_vertex_to_face(5,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
        MCCA_add_vertex_to_face(7,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);

        //if needed, add end point vertices to faces 1 and 3
        if(lambda == 0.0){
        	MCCA_add_vertex_to_face(1,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
        } else if (lambda == 1.0){
        	MCCA_add_vertex_to_face(3,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
        }
	}
}
//no intersection along +z, now check +y
//meaning we pick up vertex v4 which contributes to faces 3, 4, and 5
else {
	nvert = nvert + 1;
    x_vertex[nvert-1] = vx[jedge] ; y_vertex[nvert-1] = vy[jedge] ; z_vertex[nvert-1] = vz[jedge];
    ipath[nvert-1] = 3;
    MCCA_add_vertex_to_face(3,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
    MCCA_add_vertex_to_face(4,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
    MCCA_add_vertex_to_face(5,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);

    iedge = 3; jedge = 6;
    MCCA_check_edge(iedge,jedge,vx,vy,vz,a,b,c,d,lambda,xx,yy,zz);

    if (lambda >= 0.0 && lambda <= 1.0){
    	MCCA_check_if_vertex_exists(x_intersect,y_intersect,z_intersect,npoints,xx,yy,zz,idum);
    	//add this vertex to faces 3, 4, and 7
    	if(idum == 0){
    		npoints = npoints + 1;       
            kedge[npoints-1] = 10;
            x_intersect[npoints-1] = xx ; y_intersect[npoints-1] = yy ; z_intersect[npoints-1] = zz;

            nvert = nvert + 1;
            x_vertex[nvert-1] = x_intersect[npoints-1] ; y_vertex[nvert-1] = y_intersect[npoints-1] ; z_vertex[nvert-1] = z_intersect[npoints-1];
            ipath[nvert-1] = 3;
            MCCA_add_vertex_to_face(3,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
            MCCA_add_vertex_to_face(4,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
            MCCA_add_vertex_to_face(7,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
            //if needed, add end point vertices to faces 5 and 6
            if(lambda == 0.0){
            	MCCA_add_vertex_to_face(5,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
            } else if(lambda == 1.0){
            	MCCA_add_vertex_to_face(6,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
            }

    	}
    }
    //no intersection along +y, now check +x
    //we pick up vertex v7 which contributes to faces 3, 4, and 6
    else {
    	nvert = nvert + 1;
        x_vertex[nvert-1] = vx[jedge] ; y_vertex[nvert-1] = vy[jedge] ; z_vertex[nvert-1] = vz[jedge];
        ipath[nvert-1] = 3;
        MCCA_add_vertex_to_face(3,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
        MCCA_add_vertex_to_face(4,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
        MCCA_add_vertex_to_face(6,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);

        iedge = 6 ; jedge = 7;
        MCCA_check_edge(iedge,jedge,vx,vy,vz,a,b,c,d,lambda,xx,yy,zz);

        if (lambda >= 0.0 && lambda <= 1.0){
        	MCCA_check_if_vertex_exists(x_intersect,y_intersect,z_intersect,npoints,xx,yy,zz,idum);

        	//add this vertex to faces 3, 6, and 7
        	if(idum == 0){
        		npoints = npoints + 1;       
                kedge[npoints-1] = 11;
                x_intersect[npoints-1] = xx ; y_intersect[npoints-1] = yy ; z_intersect[npoints-1] = zz;

                nvert = nvert + 1;
                x_vertex[nvert-1] = x_intersect[npoints-1] ; y_vertex[nvert-1] = y_intersect[npoints-1] ; z_vertex[nvert-1] = z_intersect[npoints-1];
                ipath[nvert-1] = 3;
                MCCA_add_vertex_to_face(3,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
                MCCA_add_vertex_to_face(6,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
                MCCA_add_vertex_to_face(7,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);

                //if needed, add end point vertices to faces 4 and 2
                if(lambda == 0.0){
                	MCCA_add_vertex_to_face(4,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
                } else if (lambda == 1.0){
                	MCCA_add_vertex_to_face(2,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
                }

        	}
        //end of +x lambda between 0 and 1 if block
        }
    //end of +y +x lambda between 0 and 1 if block
    }
// end of +z +y +x lambda between 0 and 1 if block
}


//line intersection +z parallel to x
//black dashed line between v4 and v5 in figure cube_paths.pdf

iedge = 3; jedge = 4;
MCCA_check_edge(iedge,jedge,vx,vy,vz,a,b,c,d,lambda,xx,yy,zz);
if (lambda >= 0.0 && lambda <= 1.0){
	MCCA_check_if_vertex_exists(x_intersect,y_intersect,z_intersect,npoints,xx,yy,zz,idum);

	//add this vertex to faces 3, 5, and 7
	if(idum == 0){
		npoints = npoints + 1;       
        kedge[npoints-1] = 12;
        x_intersect[npoints-1] = xx ; y_intersect[npoints-1] = yy ; z_intersect[npoints-1] = zz;

        nvert = nvert + 1;
        x_vertex[nvert-1] = x_intersect[npoints-1] ; y_vertex[nvert-1] = y_intersect[npoints-1] ; z_vertex[nvert-1] = z_intersect[npoints-1];
        ipath[nvert-1] = 6;
        MCCA_add_vertex_to_face(3,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
        MCCA_add_vertex_to_face(5,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
        MCCA_add_vertex_to_face(7,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
        //if needed, add end point vertices to faces 4 and 2
        if(lambda == 0.0){
        	MCCA_add_vertex_to_face(4,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
        } else if(lambda == 1.0){
        	MCCA_add_vertex_to_face(2,nfaces,iface_mask,face,nface_vertices,nvert,ivertices);
        }
	}

}



//whew. now we have a list of potential faces and vertices
//-------------------------------------------------------------------------------------------------------------------     
// now address various pathologies
// the utility of face_mask is thus revealed.
//
// in the degenerate case, the plane intersects one of the faces.
// remove the double counting here.
// construct a unique integer from the vertex numbers of the intersecting plane









idum = 0;
for(int i=0; i<nface_vertices[6]; i++){
	idum = idum + pow(10,i)*ivertices[i][6];
}

//construct a unique integer from the vertex numbers of the other faces
for(int m= 0; m<nfaces; m++){
	jdum = 0;
	for(int n=0; n<nface_vertices[face[m]-1]; n++){
		jdum = jdum + pow(10,n)*ivertices[n][face[m]-1];
	}


	//zero the face_mask if the two integers match
	if(jdum == idum && (face[m]-1) != 6){
		iface_mask[face[m]-1] = 0;
	}

///CHECK THESE BEFORE CONTINUING
	// another case to address is no intersection at all and only cube vertices.
	// remove the faces where the number of vertices is less than 3
	if (iface_mask[face[m]-1] == 1 && nface_vertices[face[m]-1] < 3){
       iface_mask[face[m]-1] = 0;
     }
    }



//----------------------------------------------------------------------------------------------------------------------


//compute the volume 
//if we got any intersection


if(npoints > 0)
{
	volume = 0.0;

	//loop over valid faces
	for(int m=3; m<nfaces; m++){

		if(iface_mask[face[m]-1] == 1){
			//the vertices of the intersecting plane, by construction, are in clockwise order.
            //the cube face intersections do not have a consistent winding; 
            //we put them in clockwise order by reversing paths 

            //face 1
            //reverse order of vertexes along path 2, dark gray in cube_paths.pdf

            if ((face[m]) == 1){
            	idum = 0;
            	for (int n = 0; n<nface_vertices[face[m]-1]; n++){
            		kdum = ivertices[n][face[m]-1];
            		if (ipath[kdum-1] == 0 || ipath[kdum-1] == 1 || ipath[kdum-1] == 4){
            			idum = idum + 1;
            			xdum[idum-1] = x_vertex[kdum-1] ; ydum[idum-1] = y_vertex[kdum-1] ; zdum[idum-1] = z_vertex[kdum-1];
            		}
            	}

            	for (int n = nface_vertices[face[m]-1]-1; n>-1; n--){
            		kdum = ivertices[n][face[m]-1];
            		if (ipath[kdum-1] == 2){
            			idum = idum + 1;
            			xdum[idum-1] = x_vertex[kdum-1] ; ydum[idum-1] = y_vertex[kdum-1] ; zdum[idum-1] = z_vertex[kdum-1];
            		}
            	}
            //face 2
            //reverse order of vertexes along paths 2 and 4, dark gray and dotted light gray in cube_paths.pdf
            }else if(face[m]==2){
            	idum = 0;
            	for (int n = 0; n<nface_vertices[face[m]-1]; n++){
            		kdum = ivertices[n][face[m]-1];
            		if (ipath[kdum-1] == 1){
            			idum = idum + 1;
            			xdum[idum-1] = x_vertex[kdum-1] ; ydum[idum-1] = y_vertex[kdum-1] ; zdum[idum-1] = z_vertex[kdum-1];
            		}
            	}

            	for (int n = nface_vertices[face[m]-1]-1; n>-1; n--){
            		kdum = ivertices[n][face[m]-1];
            		if (ipath[kdum-1] == 2 || ipath[kdum-1] == 4){
            			idum = idum + 1;
            			xdum[idum-1] = x_vertex[kdum-1] ; ydum[idum-1] = y_vertex[kdum-1] ; zdum[idum-1] = z_vertex[kdum-1];
            		}
            	}
            //face 3
            //reverse order of vertexes along path 3 , black in cube_paths.pdf
            }else if(face[m]==3){
            	idum = 0;
            	for (int n = 0; n<nface_vertices[face[m]-1]; n++){
            		kdum = ivertices[n][face[m]-1];
            		if (ipath[kdum-1] == 1 || ipath[kdum-1] == 6){
            			idum = idum + 1;
            			xdum[idum-1] = x_vertex[kdum-1] ; ydum[idum-1] = y_vertex[kdum-1] ; zdum[idum-1] = z_vertex[kdum-1];
            		}
            	}

            	for (int n = nface_vertices[face[m]-1]-1; n>-1; n--){
            		kdum = ivertices[n][face[m]-1];
            		if (ipath[kdum-1] == 3){
            			idum = idum + 1;
            			xdum[idum-1] = x_vertex[kdum-1] ; ydum[idum-1] = y_vertex[kdum-1] ; zdum[idum-1] = z_vertex[kdum-1];
            		}
            	}
            //face 4
            //reverse order of vertexes along path 3 , black in cube_paths.pdf 
            }else if(face[m]==4){
            	idum = 0;
            	for (int n = 0; n<nface_vertices[face[m]-1]; n++){
            		kdum = ivertices[n][face[m]-1];
            		if (ipath[kdum-1] == 0 || ipath[kdum-1] == 2 || ipath[kdum-1] == 5){
            			idum = idum + 1;
            			xdum[idum-1] = x_vertex[kdum-1] ; ydum[idum-1] = y_vertex[kdum-1] ; zdum[idum-1] = z_vertex[kdum-1];
            		}
            	}

            	for (int n = nface_vertices[face[m]-1]-1; n>-1; n--){
            		kdum = ivertices[n][face[m]-1];
            		if (ipath[kdum-1] == 3){
            			idum = idum + 1;
            			xdum[idum-1] = x_vertex[kdum-1] ; ydum[idum-1] = y_vertex[kdum-1] ; zdum[idum-1] = z_vertex[kdum-1];
            		}
            	}
            //face 5
            //reverse order of vertexes along paths 3 and 6, black and dark dashed in cube_paths.pdf
            }else if(face[m]==5){
            	idum = 0;
            	for (int n = 0; n<nface_vertices[face[m]-1]; n++){
            		kdum = ivertices[n][face[m]-1];
            		if (ipath[kdum-1] == 0 || ipath[kdum-1] == 1){
            			idum = idum + 1;
            			xdum[idum-1] = x_vertex[kdum-1] ; ydum[idum-1] = y_vertex[kdum-1] ; zdum[idum-1] = z_vertex[kdum-1];
            		}
            	}

            	for (int n = nface_vertices[face[m]-1]-1; n>-1; n--){
            		kdum = ivertices[n][face[m]-1];
            		if (ipath[kdum-1] == 3 || ipath[kdum-1] == 6){
            			idum = idum + 1;
            			xdum[idum-1] = x_vertex[kdum-1] ; ydum[idum-1] = y_vertex[kdum-1] ; zdum[idum-1] = z_vertex[kdum-1];
            		}
            	}


            //face 6
            //reverse order of vertexes along paths 3 and 5, black and gray in cube_paths.pdf 
            }else if(face[m]==6){
            	idum = 0;
            	for (int n = 0; n<nface_vertices[face[m]-1]; n++){
            		kdum = ivertices[n][face[m]-1];
            		if (ipath[kdum-1] == 2){
            			idum = idum + 1;
            			xdum[idum-1] = x_vertex[kdum-1] ; ydum[idum-1] = y_vertex[kdum-1] ; zdum[idum-1] = z_vertex[kdum-1];
            		}
            	}

            	for (int n = nface_vertices[face[m]-1]-1; n>-1; n--){
            		kdum = ivertices[n][face[m]-1];
            		if (ipath[kdum-1] == 3 || ipath[kdum-1] == 5){
            			idum = idum + 1;
            			xdum[idum-1] = x_vertex[kdum-1] ; ydum[idum-1] = y_vertex[kdum-1] ; zdum[idum-1] = z_vertex[kdum-1];
            		}
            	}

            //face 7 is the intersecting plane and the final face
            //they are in clockwise order by construction
            }else if(face[m]==7){
            	for (int n = 0; n<nface_vertices[face[m]-1]; n++){
            		kdum = ivertices[n][face[m]-1];
           			xdum[n] = x_vertex[kdum-1] ; ydum[n] = y_vertex[kdum-1] ; zdum[n] = z_vertex[kdum-1];
            	}

            }
            //volume of this polyhedra
            //choose vertex 1 as the reference height of the pyramids.
            //faces 1, 4, and 5 then contribute zero volume because their heights are zero.
            //only faces 2, 3, 6 and 7 contribute.
            dv = MCCA_pyramid_volume(xdum, ydum, zdum, nface_vertices[face[m]-1], x_vertex[0], y_vertex[0], z_vertex[0]);
            
            volume = volume + dv;
        }
    }
}


x_final[0] = x_intersect[0]; x_final[1] = x_intersect[1]; x_final[2] = x_intersect[2]; x_final[3] = x_intersect[3]; x_final[4] = x_intersect[4]; x_final[5] = x_intersect[5];
y_final[0] = y_intersect[0]; y_final[1] = y_intersect[1]; y_final[2] = y_intersect[2]; y_final[3] = y_intersect[3]; y_final[4] = y_intersect[4]; y_final[5] = y_intersect[5];
z_final[0] = z_intersect[0]; z_final[1] = z_intersect[1]; z_final[2] = z_intersect[2]; z_final[3] = z_intersect[3]; z_final[4] = z_intersect[4]; z_final[5] = z_intersect[5];



}


void MCCA_Iterative_Area(double &alpha, double &theta, double &phi, double &area)
{


double x_final[6];
double y_final[6];
double z_final[6];
double target_volume = alpha;
double dl_max = (cos(theta)+sin(theta))*sin(phi) + cos(phi);
double b = dl_max;
double a = 0;
double p = (a+b)/2;

int npoints;
double volume;
double error = 1.0;
double etol = 1e-3;






while (error>etol){
    MCCA_volume_fraction(p,theta,phi,volume,npoints,x_final,y_final,z_final);
    error = abs(volume - target_volume);
    if(abs(error>etol)){
            if(volume < target_volume){
                a = p;
                p = (a+b)/2;
            }
            else{
                b = p;
                p = (a+b)/2;
            }
    }
}


if(npoints > 0){
    double x_[npoints];
    double y_[npoints];
    double z_[npoints];

    for(int i = 0; i<npoints; i++){
        x_[i]=x_final[i];
        y_[i]=y_final[i];
        z_[i]=z_final[i];
    }
    area = MCCA_irregular_polygon_area(x_,y_,z_,npoints);
}

}

