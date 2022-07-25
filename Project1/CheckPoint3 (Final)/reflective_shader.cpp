#include "reflective_shader.h"
#include "ray.h"
#include "render_world.h"

vec3 Reflective_Shader::
Shade_Surface(const Ray& ray,const vec3& intersection_point,
    const vec3& normal,int recursion_depth) const
{
    vec3 color;
     // determine the color
    vec3 v = (ray.endpoint-intersection_point).normalized();
    Ray reflection(intersection_point, -(v - (2.0*dot(v, normal) * normal)).normalized());

    color = ((1-reflectivity) * shader->Shade_Surface(reflection, intersection_point, normal, recursion_depth+1));
    if(recursion_depth < world.recursion_depth_limit)
        color += (reflectivity * world.Cast_Ray(reflection, recursion_depth+1));

    return color;
}
