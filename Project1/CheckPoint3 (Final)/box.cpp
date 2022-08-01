#include <limits>
#include "box.h"

// Return whether the ray intersects this box.
bool Box::Intersection(const Ray& ray) const
{
    //https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection
    //https://gamedev.stackexchange.com/questions/18436/most-efficient-aabb-vs-ray-collision-algorithms

    //Basically, the intersection would be at Ox + tDx = Bx for x
    //If we solve for x, we would get t = (Bx-Ox)/Dx

    //int sign[3];  //Store signs so instead of potentially dividing by 0, we are multiplying by the inverse (multiplying is also probably faster than dividing?)
    double txmin, txmax, tymin, tymax, tzmin, tzmax;

    vec3 inverseDir; //getting the inverse of all the Dx, Dy, Dz
    inverseDir[0] = 1/ray.direction[0];
    inverseDir[1] = 1/ray.direction[1];
    inverseDir[2] = 1/ray.direction[2];

    //Setting txmin and txmax based on the closest x side of to the endpoint                                                //Basically just this for each axis
    txmin = inverseDir[0] >= 0 ? (lo[0] - ray.endpoint[0])*inverseDir[0] : (hi[0] - ray.endpoint[0])*inverseDir[0];         //if(inverseDir[0] >= 0){
    txmax = inverseDir[0] >= 0 ? (hi[0] - ray.endpoint[0])*inverseDir[0] : (lo[0] - ray.endpoint[0])*inverseDir[0];         //    txmin = (lo[0] - ray.endpoint[0])*inverseDir[0];
    tymin = inverseDir[1] >= 0 ? (lo[1] - ray.endpoint[1])*inverseDir[1] : (hi[1] - ray.endpoint[1])*inverseDir[1];         //    txmax = (hi[0] - ray.endpoint[0])*inverseDir[0];
    tymax = inverseDir[1] >= 0 ? (hi[1] - ray.endpoint[1])*inverseDir[1] : (lo[1] - ray.endpoint[1])*inverseDir[1];         //}else{
                                                                                                                            //    txmin = (hi[0] - ray.endpoint[0])*inverseDir[0];
    //For this case:     or the same case on the bottom right corner                                                        //    txmax = (lo[0] - ray.endpoint[0])*inverseDir[0];
    //     /                                                                                                                //}
    //    /
    //   /   _____ tymax
    //  /   |
    // /    | txmin
    if((txmin > tymax) || (tymin > txmax)) return false;

    //Resetting to make sure that tx is the min and max of outer boundaries of t
    //Uhhh, tx from now on no longer represents the x axis
    if(tymin > txmin) txmin = tymin;
    if(tymax < txmax) txmax = tymax;

    tzmin = inverseDir[2] >= 0 ? (lo[2] - ray.endpoint[2])*inverseDir[2] : (hi[2] - ray.endpoint[2])*inverseDir[2];
    tzmax = inverseDir[2] >= 0 ? (hi[2] - ray.endpoint[2])*inverseDir[2] : (lo[2] - ray.endpoint[2])*inverseDir[2];

    if((txmin > tzmax) || (tzmin > txmax)) return false;
    if(tymin > txmin) txmin = tymin; //Dunno why we need this, prob don't
    if(tymax < txmax) txmax = tymax;

    return true;
}

// Compute the smallest box that contains both *this and bb.
Box Box::Union(const Box& bb) const
{
    Box box;
    for(int i = 0; i < 3; i++){
        box.lo[i] = std::min(this->lo[i], bb.lo[i]);
        box.hi[i] = std::max(this->hi[i], bb.hi[i]);
    }
    return box;
}

// Enlarge this box (if necessary) so that pt also lies inside it.
void Box::Include_Point(const vec3& pt)
{
    for(int i = 0; i < 3; i++){
        if(pt[i] < this->lo[i])
            this->lo[i] = pt[i];
        if(pt[i] > this->hi[i])
            this->hi[i] = pt[i];
    }
}

// Create a box to which points can be correctly added using Include_Point.
void Box::Make_Empty()
{
    lo.fill(std::numeric_limits<double>::infinity());
    hi=-lo;
}

// Create a box that contains everything.
void Box::Make_Full()
{
    hi.fill(std::numeric_limits<double>::infinity());
    lo=-hi;
}
