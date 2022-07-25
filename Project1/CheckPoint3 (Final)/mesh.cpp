#include "mesh.h"
#include <iostream>
#include <fstream>
#include <string>
#include <limits>

// Consider a triangle to intersect a ray if the ray intersects the plane of the
// triangle with barycentric weights in [-weight_tolerance, 1+weight_tolerance]
static const double weight_tolerance = 1e-4;

// Read in a mesh from an obj file.  Populates the bounding box and registers
// one part per triangle (by setting number_parts).
void Mesh::Read_Obj(const char* file)
{
    std::ifstream fin(file);
    if(!fin)
    {
        exit(EXIT_FAILURE);
    }
    std::string line;
    ivec3 e;
    vec3 v;
    box.Make_Empty();
    while(fin)
    {
        getline(fin,line);

        if(sscanf(line.c_str(), "v %lg %lg %lg", &v[0], &v[1], &v[2]) == 3)
        {
            vertices.push_back(v);
            box.Include_Point(v);
        }

        if(sscanf(line.c_str(), "f %d %d %d", &e[0], &e[1], &e[2]) == 3)
        {
            for(int i=0;i<3;i++) e[i]--;
            triangles.push_back(e);
        }
    }
    number_parts=triangles.size();
}

// Check for an intersection against the ray.  See the base class for details.
Hit Mesh::Intersection(const Ray& ray, int part) const
{
    Hit result = {nullptr, 0, 0};

    // Do not return intersections where dist<small_t.
    // If part>=0 only test for intersections against the specified part.
    // If part<0 intersect against all parts.

    double currentClosestDist = std::numeric_limits<double>::max();
    double dist = 0;

    if(part >= 0){
        vec3 n = Normal(ray.endpoint, part);
        int isPerpendicular = dot(ray.direction, n);
        if(isPerpendicular == 0)
            return result;
        dist = dot(vertices[triangles[part][0]] - ray.endpoint, n)/(isPerpendicular);
        if(dist > small_t && Intersect_Triangle(ray, part, dist))
            return {this, dist, part};
    }else{
        int i = 0;
        for(auto triangle : triangles){
            vec3 n = Normal(ray.endpoint, i);
            double isPerpendicular = dot(ray.direction, n);
            if(!(isPerpendicular == 0)){
                dist = dot(vertices[triangle[0]] - ray.endpoint, n)/(isPerpendicular);
                if(dist > small_t && dist < currentClosestDist && Intersect_Triangle(ray, i, dist)){
                    currentClosestDist = dist;
                    result = {this, dist, i};
                }
            }
            i++;
        }
    }

    //********************Attempt 3************************
    // if(part >= 0){
    //     vec3 n = Normal(ray.direction, part);
    //     int isPerpendicular = dot(ray.direction, n);
    //     if(isPerpendicular == 0)
    //         return result;
    //     //double dist = (dot(n, ray.endpoint) + dot(vertices[triangles[part][0]], n)) / dot(ray.direction, n);
    //     double dist = dot(vertices[triangles[part][0]] - ray.endpoint, n)/(isPerpendicular);
    //     if(dist < small_t)
    //         return result;
    //     if(Intersect_Triangle(ray, part, dist))
    //         return {this, dist, part};
    // }else{
    //     //double currentClosestDist = std::numeric_limits<double>::max();
    //     int i = 0;
    //     for(auto triangle : triangles){
    //         vec3 n = Normal(ray.direction, i);
    //         double isPerpendicular = dot(ray.direction, n);
    //         if(isPerpendicular == 0){
    //             i++;
    //             continue;
    //         }
    //         // double dist = (dot(n, ray.endpoint) + dot(vertices[triangle[0]], n)) / dot(ray.direction, n);
    //         double dist = dot(vertices[triangle[0]] - ray.endpoint, n)/(isPerpendicular);
    //         if(dist < small_t){
    //             i++;
    //             continue;
    //         }
    //         if(Intersect_Triangle(ray, i, dist)){
    //             //currentClosestDist = dist;
    //             result = {this, dist, i};
    //         }
    //         i++;
    //     }
    // }

    //********************Attempt 2************************
    // if(part >= 0){
    //     vec3 normal = Normal(vertices[triangles[part][0]], part);
    //     double orthogonalTest = dot(ray.direction, normal);
    //     if(orthogonalTest != 0){
    //         double dist = dot(vertices[triangles[part][0]] - ray.endpoint, normal) / orthogonalTest;
    //         if(dist > small_t && Intersect_Triangle(ray, part, dist)){
    //             return {this, dist, part};
    //         }
    //     }
    // }else{
    //     double currentClosestDist = std::numeric_limits<double>::max();
    //     for(size_t i = 0; i < triangles.size(); i++){
    //         vec3 normal = Normal(vertices[triangles[i][0]], (int)i);
    //         double orthogonalTest = dot(ray.direction, normal);
    //         if(orthogonalTest != 0){
    //             double dist = dot(vertices[triangles[i][0]] - ray.endpoint, normal) / orthogonalTest;
    //             if(dist > small_t && Intersect_Triangle(ray, i, dist) && dist < currentClosestDist){
    //                 currentClosestDist = dist;
    //                 result = {this, dist, (int)i};
    //             }
    //         }
    //     }
    // }

    //********************Attempt 1************************
    // if(part < 0){
    //     double currentClosestDist = std::numeric_limits<double>::max();
    //     for(int i = 0; i < triangles.size(); i++){
    //         vec3 pos = vertices[triangles[i][0]];
    //         double t = dot(pos - ray.endpoint, Normal(ray.direction, i)) / dot(ray.direction, Normal(ray.direction, i));
    //         if(t > 0 && Intersect_Triangle(ray, i, t) && t >= small_t && t <= currentClosestDist){
    //             currentClosestDist = t;
    //             result = {this, t, i};
    //         }
    //     }
    // }else{
    //     vec3 pos = vertices[triangles[part][0]];
    //     double t = dot(pos - ray.endpoint, Normal(ray.direction, part)) / dot(ray.direction, Normal(ray.direction, part));
    //     if(t > 0 && Intersect_Triangle(ray, part, t)){
    //         return {this, t, part};
    //     }
    // }

    //********************Attempt 0************************
    // double orthogonalTest = dot(ray.direction, Normal(ray.direction, part));

    // if(orthogonalTest){
    //     vec3 pos = vertices[triangles[part][0]];
    //     double t = dot(pos - ray.endpoint, Normal(ray.direction, part)) / orthogonalTest;
    //     if(t > 0){
    //         if(part < 0){
    //             double currentClosestDist = std::numeric_limits<double>::max();
    //             for(int i = 0; i < triangles.size(); i++){
    //                 pos = vertices[triangles[i][0]];
    //                 t = dot(pos - ray.endpoint, Normal(ray.direction, i)) / dot(ray.direction, Normal(ray.direction, i));
    //                 if(Intersect_Triangle(ray, i, t) && t >= small_t && t <= currentClosestDist){
    //                     result = {this, t, i};
    //                 }
    //             }
    //         }else{
    //             if(Intersect_Triangle(ray, part, t)){
    //                 return {this, t, part};
    //             }
    //         }
    //     }
    // }
    return result; //textbook (if ray intersects with plane)
    //  e + td = a + Beta(b-a) + Gamma(c-a) for some t, Beta, Gamma
}

// Compute the normal direction for the triangle with index part.
vec3 Mesh::Normal(const vec3& point, int part) const
{
    assert(part>=0);

    //2 arrays (vectors)
    //vertices: Ax, Ay, Az, Bx, By, Bz, etc 
    //  => <double, double, double><double, double, double> (the vertices)
    //Triangles: (Ia, Ib, Ic)(Ia, Id, Ic), etc 
    //  => <int, int, int><int, int, int> (index to the vertices)
    
    vec3 vertex1 = vertices[triangles[part][1]] - vertices[triangles[part][0]];
    vec3 vertex2 = vertices[triangles[part][2]] - vertices[triangles[part][0]];
    return cross(vertex1, vertex2).normalized();
}

// This is a helper routine whose purpose is to simplify the implementation
// of the Intersection routine.  It should test for an intersection between
// the ray and the triangle with index tri.  If an intersection exists,
// record the distance and return true.  Otherwise, return false.
// This intersection should be computed by determining the intersection of
// the ray and the plane of the triangle.  From this, determine (1) where
// along the ray the intersection point occurs (dist) and (2) the barycentric
// coordinates within the triangle where the intersection occurs.  The
// triangle intersects the ray if dist>small_t and the barycentric weights are
// larger than -weight_tolerance.  The use of small_t avoid the self-shadowing
// bug, and the use of weight_tolerance prevents rays from passing in between
// two triangles.
bool Mesh::Intersect_Triangle(const Ray& ray, int tri, double& dist) const
{
    //https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/barycentric-coordinates#:~:text=To%20compute%20the%20position%20of,(barycentric%20coordinates%20are%20normalized).
    
    //triangle area is (||B-A|| * ||C-A||sin(theta))/2
    //                = (||B-A|| X ||C-A||)/2                   cause ||u X v|| = ||u||||v||sin(theta) "the magnitude of the cross product can be interpreted as the area of the parallelogram."
    
    //barycentric coordinates are therefore then 
    //                = (||B-A|| X ||P-A||)/2  /  (||B-A|| X ||C-A||)/2 
    //                = (||B-A|| X ||P-A||)  /  (||B-A|| X ||C-A||)     Getting rid of the /2 cause its both in the numerator and denominator

    //Dot product is dot(A, B) = ||A||||B||cos(theta)
    //if A and B are parallel, cos(theta) is 1

    //Let's replace (AB X AP) with D and (AB X AC) with N

    //Let's pretend the numerator is dot(D, N)
    //Since ABC and ABP are coplanar, then dot(D, E) = cos(0)||D||||N|| = ||D||||N||

    //But since we dot product the numerator, we have to do the same to the denominator
    //                = dot((AB X AP), N) / dot(N, N)

    //we are then left with 
    //                = dot((AB X AP), N) / dot(N, N) where N = (AB X AC)


    ivec3 curTriangle = triangles[tri];
    vec3 AB = vertices[curTriangle[1]] - vertices[curTriangle[0]];
    vec3 AC = vertices[curTriangle[2]] - vertices[curTriangle[0]];
    vec3 CA = vertices[curTriangle[0]] - vertices[curTriangle[2]];
    vec3 n = cross(AB, AC);
    double area = dot(n, n);
    vec3 p = ray.endpoint + (ray.direction * dist);

    double alpha = dot(cross(AB, p - vertices[curTriangle[0]]), n)/area;
    double gamma = dot(cross(CA, p - vertices[curTriangle[2]]), n)/area;
    double beta = 1.0 - alpha - gamma;
    
    //********************Attempt 0************************
    // ivec3 curTriangle = triangles[tri];
    // vec3 p = ray.endpoint + (ray.direction * dist);

    //check if parallel to plane
    // if(!dot(ray.endpoint, Normal(ray.direction, tri)))
    //     return false;

    // double areaTri = (cross(vertices[curTriangle[1]] - vertices[curTriangle[0]], vertices[curTriangle[2]] - vertices[curTriangle[0]]).magnitude());
    // double alpha = cross(p - vertices[curTriangle[0]], p - vertices[curTriangle[1]]).magnitude() / areaTri;
    // double beta = cross(p - vertices[curTriangle[1]], p - vertices[curTriangle[2]]).magnitude() / areaTri;
    // double gamma = 1 - alpha - beta; //(cross(vertices[curTriangle[2]] - p, vertices[curTriangle[0]] - p).magnitude() / 2) / areaTri;

    return alpha > -weight_tol && beta > -weight_tol && gamma > -weight_tol;
}

// Compute the bounding box.  Return the bounding box of only the triangle whose
// index is part.
Box Mesh::Bounding_Box(int part) const
{
    Box b;
    TODO;
    return b;
}
