#ifndef MCCA_FUNCTIONS_H_INCLUDED
#define MCCA_FUNCTIONS_H_INCLUDED

double MCCA_dot_product(double A[], double B[]); 
void MCCA_polygon_cross_product(double a[],double b[], double c[]);
void MCCA_get_polygon_normal(double v1[],double v2[],double v3[],double norm[]);
void MCCA_get_polygon_local_coord(double x[],double y[],double z[],int &n,double xprime[],double yprime[],double zprime[]);
double MCCA_irregular_polygon_area(double x[],double y[],double z[],int &n);
double MCCA_pyramid_height(double x[],double y[],double z[],int &n,double &px,double &py,double &pz);
double MCCA_pyramid_volume(double x[],double y[],double z[],int &n,double &px,double &py,double &pz);
void MCCA_add_vertex_to_face(int iface,int &nfaces,int (&iface_mask)[7],int (&face)[7],int (&nface_vertices)[7],int nvert,int (&ivertices)[7][7]);
void MCCA_check_edge(int iedge,int jedge,const double (&vx)[8],const double (&vy)[8],const double (&vz)[8], double a,double b,double c,double d,double &lambda,double &xx,double &yy,double &zz);
void MCCA_check_if_vertex_exists(double (&x_intersect)[6],double (&y_intersect)[6],double (&z_intersect)[6],int &n,double &xx,double &yy,double &zz,int &idum);
void MCCA_volume_fraction(double dl, double theta, double phi, double &volume, int &npoints, double (&x_final)[6], double (&y_final)[6], double (&z_final)[6]);
void MCCA_Iterative_Area(double &alpha, double &theta, double &phi, double &area);





#endif