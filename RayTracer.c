/*
  CSC418 - RayTracer code - Winter 2017 - Assignment 3&4

  Written Dec. 9 2010 - Jan 20, 2011 by F. J. Estrada
  Freely distributable for adacemic purposes only.

  Uses Tom F. El-Maraghi's code for computing inverse
  matrices. You will need to compile together with
  svdDynamic.c

  You need to understand the code provided in
  this file, the corresponding header file, and the
  utils.c and utils.h files. Do not worry about
  svdDynamic.c, we need it only to compute
  inverse matrices.

  You only need to modify or add code in sections
  clearly marked "TO DO"
*/

// the ppm images used here are from http://www.cs.cornell.edu/courses/cs664/2003fa/images/

#include "utils.h"
#include <math.h>

// A couple of global structures and data: An object list, a light list, and the
// maximum recursion depth
struct object3D *object_list;
struct pointLS *light_list;
int MAX_DEPTH;

// the original scene of assignment
void buildScene(void)
{
  // Sets up all objects in the scene. This involves creating each object,
  // defining the transformations needed to shape and position it as
  // desired, specifying the reflectance properties (albedos and colours)
  // and setting up textures where needed.
  // Light sources must be defined, positioned, and their colour defined.
  // All objects must be inserted in the object_list. All light sources
  // must be inserted in the light_list.
  //
  // To create hierarchical objects:
  //   Copy the transform matrix from the parent node to the child, and
  //   apply any required transformations afterwards.
  //
  // NOTE: After setting up the transformations for each object, don't
  //       forget to set up the inverse transform matrix!

  struct object3D *o;
  struct pointLS *l;
  struct point3D p;

  ///////////////////////////////////////
  // TO DO: For Assignment 3 you have to use
  //        the simple scene provided
  //        here, but for Assignment 4 you
  //        *MUST* define your own scene.
  //        Part of your mark will depend
  //        on how nice a scene you
  //        create. Use the simple scene
  //        provided as a sample of how to
  //        define and position objects.
  ///////////////////////////////////////

  // Simple scene for Assignment 3:
  // Insert a couple of objects. A plane and two spheres
  // with some transformations.

  // Let's add a plane
  // Note the parameters: ra, rd, rs, rg, R, G, B, alpha, r_index, and shinyness)
  o=newPlane(.05,.75,.05,.05,.55,.8,.75,1,1,2, 0.5);	// Note the plane is highly-reflective (rs=rg=.75) so we
  					// should see some reflections if all is done properly.
  					// Colour is close to cyan, and currently the plane is
  					// completely opaque (alpha=1). The refraction index is
  					// meaningless since alpha=1
  // for scene signature
  // o=newPlane(1.0, 0.0, 0.0, 0.0,.55,.8,.75,1,1,2);
  // for diffuse and ambient
  // o=newPlane(.05,.75, 0.0, 0.0,.55,.8,.75,1,1,2);

  Scale(o,6,6,1);				// Do a few transforms...
  //RotateZ(o,PI/1.20);
  RotateX(o,-PI/1.85);
  Translate(o,0,-3,10);
  invert(&o->T[0][0],&o->Tinv[0][0]);		// Very important! compute
  					// and store the inverse
  					// transform for this object!
  loadTexture(o, "./textures/plane_wood2.ppm");
  insertObject(o,&object_list);			// Insert into object list

  // Let's add a couple spheres
  o=newSphere(.05,.95,.35,.35,1,.25,.25,1,1,6, 0.5);
  // for scene signature
  // o=newSphere(1.0, 0.0, 0.0, 0.0,1,.25,.25,1,1,50);
  // for diffuse and ambient
  //o=newSphere(.05,.95, 0.0, 0.0,1,.25,.25,1,1,50);
  // for refraction
  //o=newSphere(0,.7,0.9,.9,0.1,.1,.1,0.1,1.33,96);
  //Scale(o,.75,.5,1.5);
  //RotateY(o,PI/2);
  //Translate(o,-1.45,1.1,3.5);
  Translate(o,0,-1.5,5);
  invert(&o->T[0][0],&o->Tinv[0][0]);
  loadTexture(o, "./textures/sphere_cherryblossom.ppm");
  insertObject(o,&object_list);

  o=newSphere(.05,.95,.95,.15,.75,.95,.55,1,1, 6, 0.5);
  //o=newSphere(.05,.95,.95,.75,.75,.95,.55,1,1,6);
  // for signature
  // o=newSphere(1.0, 0.0, 0.0, 0.0,.75,.95,.55,1,1,50);
  // for ambient and signature
  // o=newSphere(.05,.95, 0.0, 0.0,.75,.95,.55,1,1,50);
  // for refraction
  // o=newSphere(0,.7,0.9,.95,    0.95,.95,.95,    0.1,1.33,96, 0.5);
  //Scale(o,0.5,0.5,0.5);
  //RotateZ(o,PI/1.5);
  //Translate(o,1.75,1.25,5.0);
  Translate(o,-2.5,-1,5);
  invert(&o->T[0][0],&o->Tinv[0][0]);
  loadTexture(o, "./textures/sphere_pattern2.ppm");
  insertObject(o,&object_list);

  o=newSphere(0.4,.7,0.9,.95,    0.35,.35,.95,    0.8,1.2,96, 0.5);
  Scale(o,0.5,0.5,0.5);
  RotateZ(o,PI/1.5);
  //Translate(o,1.75,1.25,5.0);
  Translate(o,2.5,-2,5);
  invert(&o->T[0][0],&o->Tinv[0][0]);
  loadTexture(o, "./textures/sphere_pattern3.ppm");
  insertObject(o,&object_list);

  // Insert a single point light source.
  p.px=0;
  p.py=15.5;
  p.pz=-5.5;
  p.pw=1;
  l=newPLS(&p,.95,.95,.95);
  insertPLS(l,&light_list);


 // End of simple scene for Assignment 3
 // Keep in mind that you can define new types of objects such as cylinders and parametric surfaces,
 // or, you can create code to handle arbitrary triangles and then define objects as surface meshes.
 //
 // Remember: A lot of the quality of your scene will depend on how much care you have put into defining
 //           the relflectance properties of your objects, and the number and type of light sources
 //           in the scene.
}

// build new scene
void buildSceneA(void)
{
 // Sets up all objects in the scene. This involves creating each object,
 // defining the transformations needed to shape and position it as
 // desired, specifying the reflectance properties (albedos and colours)
 // and setting up textures where needed.
 // Light sources must be defined, positioned, and their colour defined.
 // All objects must be inserted in the object_list. All light sources
 // must be inserted in the light_list.
 //
 // To create hierarchical objects:
 //   Copy the transform matrix from the parent node to the child, and
 //   apply any required transformations afterwards.
 //
 // NOTE: After setting up the transformations for each object, don't
 //       forget to set up the inverse transform matrix!

 struct object3D *o;
 struct pointLS *l;
 struct point3D p;

 ///////////////////////////////////////
 // TO DO: For Assignment 3 you have to use
 //        the simple scene provided
 //        here, but for Assignment 4 you
 //        *MUST* define your own scene.
 //        Part of your mark will depend
 //        on how nice a scene you
 //        create. Use the simple scene
 //        provided as a sample of how to
 //        define and position objects.
 ///////////////////////////////////////

 // Let's add a plane placed on the bottom
 // Note the parameters: ra, rd, rs, rg, R, G, B, alpha, r_index, and shinyness)
 o=newPlane(.05,.75,.05,.05,.55,.8,.75,1,1,2, 0.5);  // Note the plane is highly-reflective (rs=rg=.75) so we
 //o=newPlane(.25,.25,.05,.05,.25,.25,.25,1,1,2);           
            // should see some reflections if all is done properly.
            // Colour is close to cyan, and currently the plane is
            // completely opaque (alpha=1). The refraction index is
            // meaningless since alpha=1
 Scale(o,4.8,3.2,1);        // Do a few transforms...
 //RotateZ(o,PI/1.20);
 RotateX(o, PI/2);
 Translate(o,0,-4.8, 8);
 invert(&o->T[0][0],&o->Tinv[0][0]);    // Very important! compute
 //loadTexture(o, "wmp.ppm");
 insertObject(o,&object_list);      // Insert into object list

 // add another plane placed on the left side
 o=newPlane(.05,.75,.05,.05,.55,.8,.75,1,1,2, 0.5);  // Note the plane is highly-reflective (rs=rg=.75) so we
 Scale(o, 3.2, 4.8,1);        // Do a few transforms...
 RotateX(o, PI);
 RotateY(o, PI/2);
 Translate(o,-4.8, 0.0, 8);
 invert(&o->T[0][0],&o->Tinv[0][0]);    // Very important! compute
 //loadTexture(o, "mandrill.ppm");
 insertObject(o,&object_list);      // Insert into object list

 // add another plane placed on the right side
 o=newPlane(.05,.75,.05,.05,.55,.8,.75,1,1,2, 0.5);  // Note the plane is highly-reflective (rs=rg=.75) so we
 Scale(o, 3.2, 4.8,1);        // Do a few transforms...
 RotateX(o, PI);
 RotateY(o, -PI/2);
 Translate(o, 4.8, 0.0, 8);
 invert(&o->T[0][0],&o->Tinv[0][0]);    // Very important! compute
 //loadTexture(o, "mandrill.ppm");
 insertObject(o,&object_list);      // Insert into object list

 // Let's add another plane with texture colours on the far side
 // Note the parameters: ra, rd, rs, rg, R, G, B, alpha, r_index, and shinyness)
 o=newPlane(.05,.75,.05,.05,.55,.8,.75, 0.0,1,2, 0.5);  // Note the plane is highly-reflective (rs=rg=.75)
 Scale(o,12,12,1);        // Do a few transforms...
 //RotateZ(o,PI/1.20);
 RotateX(o,PI);
 Translate(o,0,0,28);
 invert(&o->T[0][0],&o->Tinv[0][0]);    // Very important! compute
 // add texture
 //loadTexture(o, "mcfaddin_2.ppm");
 insertObject(o,&object_list);      // Insert into object list

 o=newPlane(.05,.75,.05,.05,.55,.8,.75, 0.0,1,2, 0.5);  // Note the plane is highly-reflective (rs=rg=.75)
 Scale(o,8,8,1);        // Do a few transforms...
 RotateX(o,PI);
 Translate(o,0,2,27.9);
 invert(&o->T[0][0],&o->Tinv[0][0]);    // Very important! compute
 // add texture
 loadTexture(o, "mcfaddin_2.ppm");
 insertObject(o,&object_list);      // Insert into object list

 // Let's add a couple spheres
 o=newSphere(.05,.95,.35,.35,1,.25,.25,1,1,50, 0.5);
 Scale(o,.4,.25,0.8);
 RotateY(o,PI/2);
 Translate(o,-1.3,1.0,3.5);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 loadTexture(o, "mandrill.ppm");
 insertObject(o,&object_list);

 o=newSphere(.05,.95,.95,.75,.75,.95,.55,1,1,50, 0.5);
 Scale(o,.25,1.0,0.5);
 RotateZ(o,PI/1.5);
 Translate(o,1.6,1.1,5.0);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 loadTexture(o, "mandrill.ppm");
 insertObject(o,&object_list);

 o=newSphere(.55,.95,.95,.05,.10,.10,.95,1,1,50, 0.5);
 Scale(o, 0.5, 0.5, 0.5);
 RotateX(o,PI/2);
 Translate(o, 4.3, -4.3, 6.0);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 //loadTexture(o, "mcfaddin_2.ppm");
 loadTexture(o, "Earth-is-flat.ppm");
 insertObject(o,&object_list);

 o=newSphere(.95,.95,.95,.01,.95,.95,.95,1,1,50, 0.5);
 Scale(o, 1.0, 1.0, 1.0);
 RotateX(o,PI/2);
 RotateY(o,PI); 
 Translate(o, 0.0, -3.0, 5.0);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 loadTexture(o, "Earth-is-flat.ppm");
 insertObject(o,&object_list);


  // add an area light source
  addAreaLight(0.5, 0.5, 0.0, 1.0, 0.0,\
                  0.5, 15.5, -5.5, 1, 1,\
                  0.95, 0.95, 0.95, &object_list, &light_list);

 // End of simple scene for Assignment 3
 // Keep in mind that you can define new types of objects such as cylinders and parametric surfaces,
 // or, you can create code to handle arbitrary triangles and then define objects as surface meshes.
 //
 // Remember: A lot of the quality of your scene will depend on how much care you have put into defining
 //           the relflectance properties of your objects, and the number and type of light sources
 //           in the scene.
}


void rtShade(struct object3D *obj, struct point3D *p, struct point3D *n, struct ray3D *ray, int depth, double a, double b, struct colourRGB *col)
{
  // This function implements the shading model as described in lecture. It takes
  // - A pointer to the first object intersected by the ray (to get the colour properties)
  // - The coordinates of the intersection point (in world coordinates)
  // - The normal at the point
  // - The ray (needed to determine the reflection direction to use for the global component, as well as for
  //   the Phong specular component)
  // - The current racursion depth
  // - The (a,b) texture coordinates (meaningless unless texture is enabled)
  //
  // Returns:
  // - The colour for this ray (using the col pointer)

  struct colourRGB tmp_col;	// Accumulator for colour components
  double R,G,B;			// Colour for the object in R G and B

  // This will hold the colour as we process all the components of
  // the Phong illumination model
  tmp_col.R=0;
  tmp_col.G=0;
  tmp_col.B=0;

  if (obj->texImg==NULL)		// Not textured, use object colour
  {
  R=obj->col.R;
  G=obj->col.G;
  B=obj->col.B;
  }
  else
  {
    // Get object colour from the texture given the texture coordinates (a,b), and the texturing function
    // for the object. Note that we will use textures also for Photon Mapping.
    obj->textureMap(obj->texImg,a,b,&R,&G,&B);
  }

  //////////////////////////////////////////////////////////////
  // TO DO: Implement this function. Refer to the notes for
  // details about the shading model.
  //////////////////////////////////////////////////////////////

  // Be sure to update 'col' with the final colour computed here!

  struct pointLS *cur_light = light_list;
  struct point3D *cur_shadow_ray_p0 = p;
  struct point3D *cur_shadow_ray_d = newPoint(0.0, 0.0, 0.0, 0.0);
  struct ray3D *cur_shadow_ray = newRay(cur_shadow_ray_p0, cur_shadow_ray_d);
  struct point3D *first_hit_p = newPoint(0.0, 0.0, 0.0, 1.0);
  struct point3D *first_hit_n = newPoint(0.0, 0.0, 0.0, 0.0); 
  struct colourRGB *phong_col = (struct colourRGB *) malloc(sizeof(struct colourRGB));
  double lambda;
  struct object3D *findFirstHit_obj;

  if (p != NULL){
    while (cur_light != NULL){
      struct point3D *light_ray = newPoint(cur_light->p0.px, cur_light->p0.py, cur_light->p0.pz, cur_light->p0.pw);
      subVectors(&cur_shadow_ray->p0, light_ray);
      memcpy(&cur_shadow_ray->d, light_ray, sizeof(struct point3D));
      lambda = -1;
      findFirstHit(cur_shadow_ray, &lambda, obj, &findFirstHit_obj, first_hit_p, first_hit_n, &a, &b);
      if (lambda > 0 && lambda < 1){
      tmp_col.R += obj->alb.ra * cur_light->col.R * R;
      tmp_col.G += obj->alb.ra * cur_light->col.G * G;
      tmp_col.B += obj->alb.ra * cur_light->col.B * B;
      } else {
      phongIllumination(cur_light, ray, cur_shadow_ray, obj, p, n, phong_col);
      tmp_col.R += R * phong_col->R;
      tmp_col.G += G * phong_col->G;
      tmp_col.B += B * phong_col->B;
      }
      cur_light = cur_light->next;
      free(light_ray);
    }
  
    col->R = tmp_col.R;
    col->G = tmp_col.G;
    col->B = tmp_col.B;

    if (depth < MAX_DEPTH){

      // // Single Reflection ray.
      // struct point3D *reflect_ray_p0 = newPoint(p->px, p->py, p->pz, 1);
      // struct point3D *offset_n = newPoint(n->px/pow(2, 20), n->py/pow(2, 20), n->pz/pow(2, 20), 0);
      // addVectors(offset_n, reflect_ray_p0);
      // struct point3D *reversed_ray_d = newPoint(-ray->d.px, -ray->d.py, -ray->d.pz, 0);
      // double vn = dot(reversed_ray_d, n);
      // struct point3D *reflect_ray_d = newPoint(2*vn*n->px - reversed_ray_d->px, 2*vn*n->py - reversed_ray_d->py, 2*vn*n->pz - reversed_ray_d->pz, 0.0);
      // normalize(reflect_ray_d);
      // struct ray3D *reflect_ray = newRay(reflect_ray_p0, reflect_ray_d);
      // rayTrace(reflect_ray, depth+1, col, obj);
      // free(reflect_ray_p0);
      // free(reflect_ray_d);
      // free(reflect_ray);
      // free(reversed_ray_d);
      // free(offset_n);
      // // end of Single Reflection ray

      // Glossy Reflection rays for A4
      double theta, phi;
      double x, y, z;
      struct point3D *u, *v;
      struct ray3D *glossy_ray;
      struct point3D *glossy_ray_p0;
      struct point3D *glossy_ray_d;
      double A2W[4][4]; // Local to World conversion matrix 
      struct colourRGB *glossy_col = (struct colourRGB *) malloc(sizeof(struct colourRGB));
      int numGlossyRay = 5;
      int i;

      struct point3D *reflect_ray_p0 = newPoint(p->px, p->py, p->pz, 1);
      struct point3D *offset_n = newPoint(n->px/pow(2, 20), n->py/pow(2, 20), n->pz/pow(2, 20), 0);
      addVectors(offset_n, reflect_ray_p0);
      struct point3D *reversed_ray_d = newPoint(-ray->d.px, -ray->d.py, -ray->d.pz, 0);
      double vn = dot(reversed_ray_d, n);
      struct point3D *reflect_ray_d = newPoint(2*vn*n->px - reversed_ray_d->px, 2*vn*n->py - reversed_ray_d->py, 2*vn*n->pz - reversed_ray_d->pz, 0.0);
      normalize(reflect_ray_d);

      u=cross(reflect_ray_d, n);
      normalize(u);
      v=cross(u, reflect_ray_d);
      normalize(v);
      // construct Area2World matrix
      A2W[0][0]=u->px;
      A2W[1][0]=u->py;
      A2W[2][0]=u->pz;
      A2W[3][0]=0;

      A2W[0][1]=v->px;
      A2W[1][1]=v->py;
      A2W[2][1]=v->pz;
      A2W[3][1]=0;

      A2W[0][2]=reflect_ray_d->px;
      A2W[1][2]=reflect_ray_d->py;
      A2W[2][2]=reflect_ray_d->pz;
      A2W[3][2]=0;

      A2W[0][3]=reflect_ray_p0->px;
      A2W[1][3]=reflect_ray_p0->py;
      A2W[2][3]=reflect_ray_p0->pz;
      A2W[3][3]=1;

      for (i=0; i<numGlossyRay; i++) {
        glossy_col->R = 0;
        glossy_col->G = 0;
        glossy_col->B = 0;
        // Choose uniformly sampled random direction to send the ray
        // theta = (PI/2)*random()*obj->roughness;
        theta = (PI/2)*random();
        phi = 2*PI*random();
        x = sin(theta)*cos(phi);
        y = sin(theta)*sin(phi);
        z = cos(theta);

        // ray direction in local coordinates
        glossy_ray_d = newPoint(x, y, z, 0);  
        // convert to global coordinates
        matVecMult(A2W, glossy_ray_d);  
        // check if the glossy_ray_d vector is below the surface
        if (dot(glossy_ray_d, n) < 0) {
          glossy_ray_d = newPoint(-x, -y, z, 0);
        }
        glossy_ray = newRay(reflect_ray_p0, reflect_ray_d);

        rayTrace(glossy_ray, depth+1, glossy_col, obj);

        col->R +=  glossy_col->R/numGlossyRay;
        col->G +=  glossy_col->G/numGlossyRay;
        col->B +=  glossy_col->B/numGlossyRay;
      }

      free(u);
      free(v);
      free(glossy_ray);
      free(glossy_ray_p0);
      free(glossy_ray_d);
      free(glossy_col);
      
      free(reflect_ray_p0);
      free(reflect_ray_d);
      free(reversed_ray_d);
      free(offset_n);
      // end of Glossy Reflection rays

   
      // Refraction ray for A4
      if (obj->alpha < 1){
        struct colourRGB refract_col;
        refract_col.R = 0;
        refract_col.G = 0;
        refract_col.B = 0;
        struct point3D *t_p0 = newPoint(p->px, p->py, p->pz, 1);
        struct point3D *offset_neg_n = newPoint(-n->px/pow(2, 20), -n->py/pow(2, 20), -n->pz/pow(2, 20), 0);
        addVectors(offset_neg_n, t_p0);
        struct point3D *t_d = newPoint(1, 1, 1, 0);
        getRefractionVector(obj->r_index, 1, &(ray->d), n, t_d);
        struct ray3D *t_ray = newRay(t_p0, t_d);

        double t_lambda = -1;
        struct point3D *out_p = newPoint(1,1,1,1);
        struct point3D *out_n = newPoint(1,1,1,0);
        obj->intersect(obj, t_ray, &t_lambda, out_p, out_n, &a, &b);
    	
      	if (t_lambda < 0){
      		printf("Refraction ray lambda < 0\n");
      		exit(0);
      	}
  	
        // Add offset to out refraction ray p0
        out_p->px = out_p->px + out_n->px/pow(2, 20);
        out_p->py = out_p->py + out_n->py/pow(2, 20);
        out_p->pz = out_p->pz + out_n->pz/pow(2, 20);

        struct point3D *neg_out_n = newPoint(-out_n->px, -out_n->py, -out_n->pz, 0);

        struct point3D *out_d1 = newPoint(1, 1, 1, 0);
        getRefractionVector(1, obj->r_index, t_d, neg_out_n, out_d1);
    	
        if (dot(out_n, out_d1) >= 0){
        	struct ray3D *out_ray1 = newRay(out_p, out_d1);

        	rayTrace(out_ray1, depth+1, &refract_col, obj);
        	
        	col->R = col->R + refract_col.R * (1-obj->alpha);
        	col->G = col->G + refract_col.G * (1-obj->alpha);
        	col->B = col->B + refract_col.B * (1-obj->alpha);
        	free(out_ray1);
        }
      	free(t_p0);
      	free(offset_neg_n);
      	free(t_d);
      	free(t_ray);
      	free(out_p);
      	free(out_n);
      	free(neg_out_n);
      	free(out_d1);
      }
    }
  }

  free(cur_shadow_ray_d);
  free(cur_shadow_ray);
  free(first_hit_p);
  free(first_hit_n);
  free(phong_col);

  return;
}


void getRefractionVector(double nt, double n, struct point3D *d, struct point3D *normal, struct point3D *t){
	
	normalize(d);
	normalize(normal);
	
	//printf("\nnt = %f, n = %f, d = (%f, %f, %f), normal = (%f, %f, %f)\n", nt, n, d->px, d->py, d->pz, normal->px, normal->py, normal->pz);
	
	t->px = d->px;
	t->py = d->py;
	t->pz = d->pz;
	
	double d_dot_n = dot(t, normal);
	struct point3D *n_times_d_dot_n = newPoint(normal->px * d_dot_n, normal->py * d_dot_n, normal->pz * d_dot_n, 0);
	
	subVectors(n_times_d_dot_n, t);
	
	t->px = t->px * n / nt;
	t->py = t->py * n / nt;
	t->pz = t->pz *n / nt;
	
	double mult = sqrt(1 - n * n * (1 - d_dot_n * d_dot_n)/(nt * nt));
	struct point3D *vec2 = newPoint(normal->px * mult, normal->py * mult, normal->pz * mult, 0);
	
	subVectors(vec2, t);
	
	normalize(t);
	
	free(n_times_d_dot_n);
	free(vec2);
	
	//printf("t = (%f, %f, %f)\n----------------------------------\n", t->px, t->py, t->pz);
	
	return;
}

void phongIllumination(struct pointLS *light, struct ray3D *ray, struct ray3D *light_ray, struct object3D *obj, struct point3D *p, struct point3D *n, struct colourRGB *col)
{
	/* Return the phone illumination value col
	 * light: light source
	 * ray: ray from view to intersection point
	 * light_ray: ray from intersection point to light source
	 * obj: the intersection object
	 * p: the intersection point
	 * n: the normal at intersection point
	 * col: the RETURN colour
	 */

	struct point3D *L = newPoint(light_ray->d.px, light_ray->d.py, light_ray->d.pz, 0.0);
	struct point3D *N = newPoint(n->px, n->py, n->pz, 0.0);
	struct point3D *V = newPoint(-ray->d.px, -ray->d.py, -ray->d.pz, 0.0);
	normalize(L);
	normalize(N);
	normalize(V);
	
	double ln = dot(L, N);
	
	struct point3D *R = newPoint(2*ln*N->px - L->px, 2*ln*N->py - L->py, 2*ln*N->pz - L->pz, 0);
	normalize(R);

	double nl = dot(N, L);
	if (nl < 0){
		nl = 0;
	}
	
	double vr = pow(dot(V, R), obj->shinyness);
	if (vr < 0){
		vr = 0;
	}
	double ambient = obj->alb.ra;
	double diffuse = obj->alb.rd * nl;
	double specular = obj->alb.rs * vr;
	col->R = (ambient + diffuse + specular) * light->col.R;
	col->G = (ambient + diffuse + specular) * light->col.G;
	col->B = (ambient + diffuse + specular) * light->col.B;
	free(L);
	free(N);
	free(V);
	free(R);
	return;
}

void findFirstHit(struct ray3D *ray, double *lambda, struct object3D *Os, struct object3D **obj, struct point3D *p, struct point3D *n, double *a, double *b)
{
 // Find the closest intersection between the ray and any objects in the scene.
 // It returns:
 //   - The lambda at the intersection (or < 0 if no intersection)
 //   - The pointer to the object at the intersection (so we can evaluate the colour in the shading function)
 //   - The location of the intersection point (in p)
 //   - The normal at the intersection point (in n)
 //
 // Os is the 'source' object for the ray we are processing, can be NULL, and is used to ensure we don't 
 // return a self-intersection due to numerical errors for recursive raytrace calls.
 //

 /////////////////////////////////////////////////////////////
 // TO DO: Implement this function. See the notes for
 // reference of what to do in here
 /////////////////////////////////////////////////////////////

  double lambda1;    // Lambda at intersection
  struct object3D *cur_obj;  // Pointer to walk through object_list
  struct point3D p1;  // Intersection point
  struct point3D n1;  // Normal at intersection

  cur_obj = object_list;
   
  if (cur_obj==NULL) {
    return;
  }

  while (cur_obj!=NULL) {
  	if (cur_obj != Os){
  		cur_obj->intersect(cur_obj, ray, &lambda1, &p1, &n1, a, b);
  		if (lambda1 >= 0 && (*lambda == -1 || lambda1 < *lambda)) {
  		*lambda = lambda1;
  		*obj = cur_obj;
  		*p = p1;
  		*n = n1;
  		}
  	}
    cur_obj=cur_obj->next;
  }
  return;
}

void rayTrace(struct ray3D *ray, int depth, struct colourRGB *col, struct object3D *Os)
{
 // Ray-Tracing function. It finds the closest intersection between
 // the ray and any scene objects, calls the shading function to
 // determine the colour at this intersection, and returns the
 // colour.
 //
 // Os is needed for recursive calls to ensure that findFirstHit will
 // not simply return a self-intersection due to numerical
 // errors. For the top level call, Os should be NULL. And thereafter
 // it will correspond to the object from which the recursive
 // ray originates.
 //
  
 double lambda;		// Lambda at intersection
 double a,b;		// Texture coordinates
 struct object3D *obj;	// Pointer to object at intersection
 struct point3D p;	// Intersection point
 struct point3D n;	// Normal at intersection
 struct colourRGB I;	// Colour returned by shading function

 if (depth>MAX_DEPTH)	// Max recursion depth reached. Return invalid colour.
 {
  // col->R=-1;
  // col->G=-1;
  // col->B=-1;
  return;
 }

 ///////////////////////////////////////////////////////
 // TO DO: Complete this function. Refer to the notes
 // if you are unsure what to do here.
 ///////////////////////////////////////////////////////

 //find the closest intersection
  lambda = -1;

  findFirstHit(ray, &lambda, NULL, &obj, &p, &n, &a, &b);

  if (lambda < 0) {
    return;
  }

  // evaluate shading mode and get the colour
  if (obj!=Os) {
    rtShade(obj, &p, &n, ray, depth, a, b, &I);

    if (Os!=NULL) {
      col->R += Os->alb.rg*I.R;
      col->G += Os->alb.rg*I.G;
      col->B += Os->alb.rg*I.B;
    } else {
      col->R += I.R;
      col->G += I.G;
      col->B += I.B;
    } 
  }

  return;
}

int main(int argc, char *argv[])
{
 // Main function for the raytracer. Parses input parameters,
 // sets up the initial blank image, and calls the functions
 // that set up the scene and do the raytracing.
 struct image *im;	// Will hold the raytraced image
 struct view *cam;	// Camera and view for this scene
 int sx;		// Size of the raytraced image
 int antialiasing;	// Flag to determine whether antialiaing is enabled or disabled
 char output_name[1024];	// Name of the output file for the raytraced .ppm image
 struct point3D e;		// Camera view parameters 'e', 'g', and 'up'
 struct point3D g;
 struct point3D up;
 double du, dv;			// Increase along u and v directions for pixel coordinates
 struct point3D pc,d;		// Point structures to keep the coordinates of a pixel and
				// the direction or a ray
 struct ray3D *ray;		// Structure to keep the ray from e to a pixel
 struct colourRGB col, sub_col;		// Return colour for raytraced pixels
 struct colourRGB background;   // Background colour
 int i,j;			// Counters for pixel coordinates
 unsigned char *rgbIm;

 if (argc<5)
 {
  fprintf(stderr,"RayTracer: Can not parse input parameters\n");
  fprintf(stderr,"USAGE: RayTracer size rec_depth antialias output_name\n");
  fprintf(stderr,"   size = Image size (both along x and y)\n");
  fprintf(stderr,"   rec_depth = Recursion depth\n");
  fprintf(stderr,"   antialias = A single digit, 0 disables antialiasing. Anything else enables antialiasing\n");
  fprintf(stderr,"   output_name = Name of the output file, e.g. MyRender.ppm\n");
  exit(0);
 }
 sx=atoi(argv[1]);
 MAX_DEPTH=atoi(argv[2]);
 if (atoi(argv[3])==0) antialiasing=0; else antialiasing=1;
 strcpy(&output_name[0],argv[4]);

 fprintf(stderr,"Rendering image at %d x %d\n",sx,sx);
 fprintf(stderr,"Recursion depth = %d\n",MAX_DEPTH);
 if (!antialiasing) fprintf(stderr,"Antialising is off\n");
 else fprintf(stderr,"Antialising is on\n");
 fprintf(stderr,"Output file name: %s\n",output_name);

 object_list=NULL;
 light_list=NULL;

 // Allocate memory for the new image
 im=newImage(sx, sx);
 if (!im)
 {
  fprintf(stderr,"Unable to allocate memory for raytraced image\n");
  exit(0);
 }
 else rgbIm=(unsigned char *)im->rgbdata;

 ///////////////////////////////////////////////////
 // TO DO: You will need to implement several of the
 //        functions below. For Assignment 3, you can use
 //        the simple scene already provided. But
 //        for Assignment 4 you need to create your own
 //        *interesting* scene.
 ///////////////////////////////////////////////////
 buildScene();		// Create a scene. This defines all the
			// objects in the world of the raytracer
 //buildSceneA();

 //////////////////////////////////////////
 // TO DO: For Assignment 3 you can use the setup
 //        already provided here. For Assignment 4
 //        you may want to move the camera
 //        and change the view parameters
 //        to suit your scene.
 //////////////////////////////////////////

 // Mind the homogeneous coordinate w of all vectors below. DO NOT
 // forget to set it to 1, or you'll rgbImget junk out of the
 // geometric transformations later on.

 // Camera center is at (0,0,-1)
 e.px=0;
 e.py=0;
 e.pz=-1;
 e.pw=1;

 // To define the gaze vector, we choose a point 'pc' in the scene that
 // the camera is looking at, and do the vector subtraction pc-e.
 // Here we set up the camera to be looking at the origin, so g=(0,0,0)-(0,0,-1)
 g.px=0;
 g.py=0;
 g.pz=1;
 g.pw=0;

 // Define the 'up' vector to be the Y axis
 up.px=0;
 up.py=1;
 up.pz=0;
 up.pw=0;
  

 // Set up view with given the above vectors, a 4x4 window,
 // and a focal length of -1 (why? where is the image plane?)
 // Note that the top-left corner of the window is at (-2, 2)
 // in camera coordinates.
 cam=setupView(&e, &g, &up, -3, -2, 2, 4);

 if (cam==NULL)
 {
  fprintf(stderr,"Unable to set up the view and camera parameters. Our of memory!\n");
  cleanup(object_list,light_list);
  deleteImage(im);
  exit(0);
 }

 // Set up background colour here
 background.R=0.01;
 background.G=0.01;
 background.B=0.01;

 // Do the raytracing
 //////////////////////////////////////////////////////
 // TO DO: You will need code here to do the raytracing
 //        for each pixel in the image. Refer to the
 //        lecture notes, in particular, to the
 //        raytracing pseudocode, for details on what
 //        to do here. Make sure you undersand the
 //        overall procedure of raytracing for a single
 //        pixel.
 //////////////////////////////////////////////////////
 du=cam->wsize/(sx-1);		// du and dv. In the notes in terms of wl and wr, wt and wb,
 dv=-cam->wsize/(sx-1);		// here we use wl, wt, and wsize. du=dv since the image is
				// and dv is negative since y increases downward in pixel
				// coordinates and upward in camera coordinates.

 fprintf(stderr,"View parameters:\n");
 fprintf(stderr,"Left=%f, Top=%f, Width=%f, f=%f\n",cam->wl,cam->wt,cam->wsize,cam->f);
 fprintf(stderr,"Camera to world conversion matrix (make sure it makes sense!):\n");
 printmatrix(cam->C2W);
 fprintf(stderr,"World to camera conversion matrix\n");
 printmatrix(cam->W2C);
 fprintf(stderr,"\n");

 fprintf(stderr,"Rendering row: ");
 for (j=0;j<sx;j++)		// For each of the pixels in the image
 {
  fprintf(stderr,"%d/%d, ",j,sx);
  for (i=0;i<sx;i++)
  {
	
	// antialiasing: the pixel into n by n subpixels.
	if (antialiasing == 1){
		int n = 3;
		double sub_du = du/n;
		double sub_dv = dv/n;
		
		col.R=0.0;
		col.G=0.0;
		col.B=0.0;
		
		// Generate n^2 rays to get the sum color.
		for (int iter_v = 0; iter_v < n; iter_v++){
			for (int iter_h = 0; iter_h < n; iter_h++){ 
				pc.px = cam->wl + i*du + iter_h*sub_du + rand()/(RAND_MAX/sub_du);
				pc.py = cam->wt + j*dv + iter_v*sub_dv + rand()/(RAND_MAX/sub_dv);
				pc.pz = cam->f;
				pc.pw = 1.0;
				
				d.px = pc.px-cam->e.px;
				d.py = pc.py-cam->e.py;
				d.pz = pc.pz-cam->e.pz;
				d.pw = 0.0;
				
				matVecMult(cam->C2W, &pc);
				matVecMult(cam->C2W, &d);
				
				ray=newRay(&pc, &d);
				
				sub_col.R=0.0;
				sub_col.G=0.0;
				sub_col.B=0.0;
				
				rayTrace(ray, 0, &sub_col, NULL);
				free(ray);
				
				col.R += sub_col.R;
				col.G += sub_col.G;
				col.B += sub_col.B;
			}
		}
		
		// pixel color is the average of all sub pixel colors.
		col.R = col.R / (n*n);
		col.G = col.G / (n*n);
		col.B = col.B / (n*n);
	} else {
	
    ///////////////////////////////////////////////////////////////////
    // TO DO - complete the code that should be in this loop to do the
    //         raytracing!
    ///////////////////////////////////////////////////////////////////

    // the coordinates of a pixel in view coordinator
	 
		pc.px = cam->wl+(i+0.5)*du;
		pc.py = cam->wt+(j+0.5)*dv;
		pc.pz = cam->f;
		pc.pw = 1.0;

		// the direction of a ray in view coordinator
		d.px = pc.px-cam->e.px;
		d.py = pc.py-cam->e.py;
		d.pz = pc.pz-cam->e.pz;
		d.pw = 0.0;
		// convert to world-space
		matVecMult(cam->C2W, &pc);
		matVecMult(cam->C2W, &d);
		//construct viewing ray
		ray=newRay(&pc, &d);

		// initialize pixel colour
		col.R=0.0;
		col.G=0.0;
		col.B=0.0;

		// call rayTrace
		rayTrace(ray, 0, &col, NULL);
		free(ray);
	}
    
    if (col.R <= 0) {
      rgbIm[(j*sx + i)*3] = background.R * 255;
      rgbIm[(j*sx + i)*3 + 1] = background.G * 255;
      rgbIm[(j*sx + i)*3 + 2] = background.B * 255;
    } 
    else {
	  if (col.R > 1){
		rgbIm[(j*sx + i)*3] = 255;
	  }
	  else{
		rgbIm[(j*sx + i)*3] = col.R * 255;
	  }
	  if (col.G > 1){
		rgbIm[(j*sx + i)*3+1] = 255;
	  }
	  else{
		rgbIm[(j*sx + i)*3+1] = col.G * 255;
	  }
	  if (col.B > 1){
		rgbIm[(j*sx + i)*3+2] = 255;
	  }
	  else{
		rgbIm[(j*sx + i)*3+2] = col.B * 255;
	  }
    }
  } // end for i
 } // end for j

 fprintf(stderr,"\nDone!\n");

 // Output rendered image
 imageOutput(im,output_name);

 // Exit section. Clean up and return.
 cleanup(object_list,light_list);		// Object and light lists
 deleteImage(im);				// Rendered image
 free(cam);					// camera view
 exit(0);
}
