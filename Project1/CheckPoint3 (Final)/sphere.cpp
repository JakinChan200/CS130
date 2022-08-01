#include "sphere.h"
#include "ray.h"

// Determine if the ray intersects with the sphere
Hit Sphere::Intersection(const Ray& ray, int part) const
{
    vec3 distanceCenter = ray.endpoint-center;
    double discriminant = pow(dot(ray.direction, distanceCenter), 2) - dot(ray.direction, ray.direction) * (dot(distanceCenter, distanceCenter) - pow(radius, 2));

    //If there is a solution (granted t is >= small_t)
    if(discriminant >= 0){
        double solution1 = -(dot(ray.direction, distanceCenter) + sqrt(discriminant)) / dot(ray.direction, ray.direction);
        double solution2 = -(dot(ray.direction, distanceCenter) - sqrt(discriminant)) / dot(ray.direction, ray.direction);

        if(solution1 <= solution2 && solution1 >= small_t){ //-b/2a also works for when discriminant == 0 too (note for later?)
            return {this, solution1, part};
        }else if(solution2 >= small_t){
            return {this, solution2, part};
        }
    }
    return {nullptr, 0, part};
}

vec3 Sphere::Normal(const vec3& point, int part) const
{
    vec3 normal;
    normal = -(center-point).normalized();
    //compute the normal direction

    // Return the normal.  For objects with multiple parts (meshes), you
    // will need to use part to determine which piece was intersected.
    // It will be set to the part structure entry returned from the
    // intersection routine.
    return normal;
}

Box Sphere::Bounding_Box(int part) const
{
    Box box;
    box.lo = center - radius;
    box.hi = center + radius;
    return box;
}

